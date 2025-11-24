# CONTRASTS
library(lme4)
library(MASS)

get_bf_contrast <- function(N=100, 
                            t.points=c(0,1,5,10), 
                            hypothesis = "W[c0](t1_t3)", 
                            betas=c(0, 0.4, 0.5, 1, 0, 1, 3, 5),
                            var.u0=.1,
                            var.u1=.1,
                            var.e=.2,
                            cov=.01,
                            seed=NULL,
                            fraction=1
){
  
  if(!is.null(seed)) {set.seed(seed)}
  
  ################################################################################
  ## Data Generation
  n.conditions <- length(betas)/length(t.points)
  
  Sigma_u <- matrix(c(var.u0, cov,
                      cov, var.u1), 2, 2)
  
  n.time <- length(t.points)
  n.cells <- n.conditions * n.time
  
  # Build full design
  cond <- 0:(n.conditions - 1)
  dat <- expand.grid(id = 1:N,  cond = cond, t = t.points)
  dat <- dat[order(dat$id, dat$cond, dat$t), ]
  
  # Create dummy variables, ordered by condition
  dummy_names <- paste0("D", 1:n.cells)
  dat[dummy_names] <- 0
  
  counter <- 1
  for (c in cond) {
    for (tp in seq_along(t.points)) {
      dat[dummy_names[counter]] <- as.numeric(dat$t == t.points[tp] & dat$cond == c)
      counter <- counter + 1
    }
  }
  
  # Random effects
  u <- MASS::mvrnorm(N, mu = c(0, 0), Sigma = Sigma_u)
  dat$u0 <- u[dat$id, 1]
  dat$u1 <- u[dat$id, 2]
  
  # Outcome (exactly as in your original version)
  y_terms <- paste0("betas[", seq_along(betas), "] * ", dummy_names, collapse = " + ")
  dat$y <- with(dat, eval(parse(text = paste0(
    y_terms, " + u0 + u1 * t + rnorm(nrow(dat), 0, sqrt(var.e))"
  ))))
  
  ################################################################################
  ## Fit MLM
  formula_str <- paste("y ~ -1 +", paste(dummy_names, collapse = " + "), "+ (1 + t|id)")
  mod <- lme4::lmer(as.formula(formula_str), data = dat, REML = F)
  
  beta_hat <- lme4::fixef(mod)
  var_beta_hat <- vcov(mod)
  
  # Name betas for condition × time
  beta_names <- dummy_names
  cond_labels <- paste0("cond", rep(cond, each = n.time))
  time_labels <- paste0("t", rep(seq_len(n.time), times = n.conditions))
  names(beta_hat) <- paste0(cond_labels, "_", time_labels)
  

  ################################################################################
  ## Extract Conditions and Timepoints from Names
  
  cond_labels <- paste0("c", cond)
  time_labels <- paste0("t", seq_along(t.points))
  beta_names <- unlist(lapply(cond_labels, function(c)
    paste0(c, "_", time_labels)
  ))
  names(beta_hat) <- beta_names
  
  
  ################################################################################
  ## Parse Hypothesis to Identify Needed Contrasts (robust token extraction)
  
  # Extract all contrast identifiers from hypothesis
  contrast_patterns <- gregexpr("[WB]\\[[^]]+\\]\\([^)]*\\)", hypothesis, perl = TRUE)
  contrast_names <- unique(regmatches(hypothesis, contrast_patterns)[[1]])
  
  contrasts_list <- list()
  
  for (cn in contrast_names) {
    type <- substr(cn, 1, 1)                        # 'W' or 'B'
    cond_str <- sub(".*\\[([^]]+)\\].*", "\\1", cn) # inside [...]
    conds <- unlist(strsplit(cond_str, "_"))
    
    time_str <- sub(".*\\(([^)]*)\\).*", "\\1", cn) # inside (...)
    time_parts <- unlist(strsplit(time_str, "_"))   # split groups by '_'
    
    # mean_sets <- lapply(time_parts, function(x) {
    #   ts <- regmatches(x, gregexpr("t[0-9]+", x, perl = TRUE))[[1]]
    #   if (length(ts) == 0) {
    #     stop("No time tokens found in part: '", x, "'. Use t# or mean(t#t#).")
    #   }
    #   ts
    # })
    # mean_sets is a list where each element is a vector of "tX" tokens
    
    # Build contrast vector
    v <- rep(0, length(beta_hat))
    
    if (type == "W") {
      # Within-condition: require two groups (e.g. mean(...) _ t3)
      if (length(mean_sets) < 2) {
        stop("Within-condition contrast requires two time sets, e.g. W[c0](mean(t1t2)_t3)")
      } else if (length(conds)>1) {
        stop("Within-condition contrast requires only one condition, e.g. W[c0](...)")
      }
      cond <- conds
      tset1 <- time_parts[1]
      tset2 <- time_parts[2]
      for (t in tset1) {v[match(paste0(cond, "_", t), names(beta_hat))] <- -1 / length(tset1)}
      for (t in tset2) {v[match(paste0(cond, "_", t), names(beta_hat))] <-  1 / length(tset2)}
      
    } else if (type == "B") {
      # Between-condition: compare same time-set across two conditions
      if (length(conds) < 2) stop("Between contrast needs two conditions, e.g. B[c0_c1](mean(t1t2))")
      cond1 <- conds[1]
      cond2 <- conds[2]
      tset <- mean_sets   
      for (t in tset) {v[match(paste0(cond1, "_", t), names(beta_hat))] <- -1 / length(tset)}
      for (t in tset) {v[match(paste0(cond2, "_", t), names(beta_hat))] <-  1 / length(tset)}
    }
    
    contrasts_list[[cn]] <- v
  }
  
  ################################################################################
  ## Bayes Factor Computation
  
  
  relevant_c <- unlist(strsplit(hypothesis, "\\s*[<>=]\\s*"))
  relevant_c <- relevant_c[relevant_c %in% rownames(c)]
  
  if (length(relevant_c) == 0) stop("No matching contrasts found in hypothesis.")
  
  est <- drop(c %*% beta_hat)[relevant_c]
  var_est_all <- c %*% var_beta_hat %*% t(c)
  var_est <- var_est_all[relevant_c, relevant_c]
  
  var_matrices <- lapply(1:length(est), function(i) {
    matrix(var_est[i, i], nrow = 1, ncol = 1,
           dimnames = list(relevant_c[i], relevant_c[i]))
  })
  
  bf_res <- bain::bain(
    x = est,
    Sigma = var_matrices,
    hypothesis = hypothesis,
    n = rep(N, length(est)),
    group_parameters = 1,
    joint_parameters = 0,
    fraction = fraction
  )
  
  bf_c <- bf_res[["fit"]][["BF.c"]][1]
  
  ################################################################################
  ################################################################################
  
  hypothesis <- "B[c0_c1](t1_t2) > W[c0](t2_t3) & B[c0_c1](mean(t1t2)) < 0"
  
  # Find all contrast-like expressions
  pattern <- "[WB]\\[[^]]+\\]\\([^)]*\\)"
  contrasts <- unique(regmatches(hypothesis, gregexpr(pattern, hypothesis, perl = TRUE))[[1]])
  
  # Assign short names x1, x2, ...
  short_names <- paste0("x", seq_along(contrasts))
  
  # Mapping table
  mapping <- setNames(short_names, contrasts)
  
  # Replace in hypothesis
  hypothesis_short <- hypothesis
  for(i in seq_along(contrasts)) {
    hypothesis_short <- gsub(contrasts[i], short_names[i], hypothesis_short, fixed = TRUE)
  }
  
  mapping
  hypothesis_short
  
  
  
  
  
  
  
  return(bf_c)
}

