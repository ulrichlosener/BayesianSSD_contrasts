# CONTRASTS
library(lme4)
library(MASS)
library(matrixcalc) # for is.singular
library(restriktor)

get_bf_contrast <- function(N=100, 
                            t.points=c(0,2,6), 
                            hypothesis="maintenance_diff>0",
                            betas=c(0, 0.1, 0.5, 0, .5, 1.5),
                            var.u0=.1,
                            var.u1=.1,
                            var.e=.2,
                            cov=.01,
                            seed=NULL,
                            fraction=1,
                            attrition=c(1, .8, .6),
                            method="bfc"
                    ){

  if(!is.null(seed)) {set.seed(seed)}
  
  ################################################################################
  ## Data Generation
  
  treat <- 0:1
  n <- length(t.points)
  
  # Create full design
  dat <- data.frame(id=rep(1:N, each=n), t=rep(t.points, N), j=rep(c(1,2,3), N), treat=rep(treat, each=N*n/2))
  dat$treat <- factor(dat$treat)
  X <- model.matrix(~ t:treat - 1, data = dat)
  
  # Dummies for each time × condition combination
  dat$D1 <- as.numeric(dat$t == t.points[1] & dat$treat == 0)
  dat$D2 <- as.numeric(dat$t == t.points[2] & dat$treat == 0)
  dat$D3 <- as.numeric(dat$t == t.points[3] & dat$treat == 0)
  dat$D4 <- as.numeric(dat$t == t.points[1] & dat$treat == 1)
  dat$D5 <- as.numeric(dat$t == t.points[2] & dat$treat == 1)
  dat$D6 <- as.numeric(dat$t == t.points[3] & dat$treat == 1)
  
  # Random effects variance-covariance matrix
  sigma_u <- matrix(c(var.u0, cov,
                      cov, var.u1), 2, 2)
  
  # if sigma_u is not positive definite, use the nearest PD matrix
  if(!matrixcalc::is.positive.definite(sigma_u)){
    sigma_u <- Matrix::nearPD(sigma_u)$mat
  }
  
  # draw random effects
  u <- MASS::mvrnorm(N, mu = c(0, 0), Sigma = sigma_u)
  dat$u0 <- u[dat$id, 1]
  dat$u1 <- u[dat$id, 2]
  
  # compute outcome y
  dat$y <- with(dat,
                betas[1]*D1 + betas[2]*D2 + betas[3]*D3 +
                betas[4]*D4 + betas[5]*D5 + betas[6]*D6 +
                u0 + u1*t + rnorm(nrow(dat), 0, sqrt(var.e))
    )

  # introduce attrition
  if(isFALSE(attrition)){
    attr <- matrix(rep(1, 6), nrow = 2, byrow = TRUE)
  } else if (length(attrition) == length(t.points)) {
      attr <- matrix(rep(attrition, 2), nrow = 2, byrow = TRUE)
  } else if (length(attrition) == 2 * length(t.points)) {
      attr <- matrix(attrition, nrow = 2, byrow = TRUE)
  } else {
      stop("length of 'attrition' vector must be the same or double the length of 't.points'.")
  }

  # dropout probability per timepoint (FACTOR-SAFE)
  dat$drop <- 1 - attr[cbind(as.integer(dat$treat), dat$j)]
  
  # introduce NAs
  dat$y <- ifelse(
    rbinom(nrow(dat), 1, dat$drop) == 1,
    NA,
    dat$y
  )
  
  # enforce monotonous dropout (no intermittent missingness)
  first_na <- tapply(
    dat$j[is.na(dat$y)],
    dat$id[is.na(dat$y)],
    min
  )
  
  drop_time <- first_na[as.character(dat$id)]
  dat$y[!is.na(drop_time) & dat$j >= drop_time] <- NA
  
  ################################################################################
  ## Fit MLM
  
  simplified <- FALSE # indication whether model simplification had to be done due to high attrition
  
  suppressMessages({ # suppress messages about singular fit
    mod <- tryCatch({
      lme4::lmer(y ~ -1 + D1 + D2 + D3 + D4 + D5 + D6 + (1 + t|id), 
                 data = dat,
                 REML = F,
                 control = lme4::lmerControl(calc.derivs = F))
    }, error = function(e) {
      if (grepl("Downdated VtV is not positive definite", e$message) ||
          grepl("number of observations.*number of random effects", e$message)) {
        simplified <<- TRUE # indicate when simplification happens
        tryCatch({
          lme4::lmer(formula = y ~ -1 + D1 + D2 + D3 + D4 + D5 + D6 + (1 | id) + (0 + t | id), # independent random effects
                     data = dat,
                     REML = F,
                     control = lme4::lmerControl(calc.derivs = F))
        })
      }
    })
  })

  # extract estimates and variance-covariance matrix
  beta_hat <- lme4::fixef(mod)
  var_beta_hat <- vcov(mod)
  
  if(length(beta_hat) < 6) {
    return(list(
      bf_c = NA,
      PMPc = NA,
      bf12 = NA,
      goric12 = NA,
      goric_c = NA,
      goric_weight = NA,
      simplified = simplified
    ))
  }
  
  ################################################################################
  ## Contrasts
  
  # create all contrasts of interest and name them
  c <- rbind(
    # effects within conditions
    "initial_0" = c(-1, 1, 0, 0, 0, 0),
    "initial_1" = c(0, 0, 0, -1, 1, 0),
    "maintenance_0" = c(0, -1, 1, 0, 0, 0),
    "maintenance_1" = c(0, 0, 0, 0, -1, 1),
    "overall_0"   = c(-1, 0, 1, 0, 0, 0),
    "overall_1"   = c(0, 0, 0, -1, 0, 1),
    # comparisons of means between conditions
    "avg_diff"  = c(-1/3, -1/3, -1/3, 1/3, 1/3, 1/3),
    "t1_diff"   = c(-1, 0, 0, 1, 0, 0),
    "t2_diff"   = c(0, -1, 0, 0, 1, 0),
    "t3_diff"   = c(0, 0, -1, 0, 0, 1),
    # comparisons of effects between conditions
    "initial_diff" = c(0, 0, 0, -1, 1, 0) - c(-1, 1, 0, 0, 0, 0), # initial_1 - initial_0
    "maintenance_diff" = c(0, 0, 0, 0, -1, 1) - c(0, -1, 1, 0, 0, 0), # maintenance_1 - maintenance_0
    "overall_diff"   = c(0, 0, 0, -1, 0, 1) - c(-1, 0, 1, 0, 0, 0)  # overall_1   - overall_0
  )
  
  # re-structure hypothesis object
  hyp <- paste(unlist(hypothesis), collapse = " ; ")

  # extract contrasts present in hypothesis
  rel_c <- unique(unlist(strsplit(hyp, "\\s*[<>=;]\\s*")))
  
  # keep only elements containing letters (no numbers)
  rel_c <- gsub("\\b[-+]?[0-9.]+\\s*\\*\\s*", "", rel_c)
  relevant_c <- unique(trimws(rel_c[grepl("[A-Za-z]", rel_c)]))
  
  relevant_c <- intersect(relevant_c, rownames(c))
  
  if (length(relevant_c) == 0) {
    return(list(
      bf_c = NA,
      PMPc = NA,
      bf12 = NA,
      goric12 = NA,
      goric_c = NA,
      goric_weight = NA,
      simplified = simplified
    ))
  }
  
  # compute estimates and their variance
  est <- drop(c %*% beta_hat)[relevant_c]
  var_est_all <- c %*% var_beta_hat %*% t(c)
  var_est <- as.matrix(var_est_all[relevant_c, relevant_c])
  
  # format Sigma for bain
  var_matrices <- lapply(1:length(est), function(i) {
    matrix(var_est[i, i], nrow = 1, ncol = 1,
           dimnames = list(relevant_c[i], relevant_c[i]))
  })
  
  ##############################################################################
  ## bain
  
  n_eff <- get_neff_contrast(model = mod)
  
  bf_res <- bain::bain(x = est, Sigma = var_matrices, 
                 hypothesis = hyp, n = rep(N, length(est)), 
                 group_parameters = 1, joint_parameters = 0, fraction=fraction)

  ##############################################################################
  ## GORIC(A)
  
  if(tolower(method) == "gorica" || tolower(method) == "both"){
    # re-structure hypothesis
    hyp_goric <- unlist(hypothesis, use.names = FALSE)
    
    # Create names H1, H2, etc
    names(hyp_goric) <- paste0("H", seq_along(hyp_goric))
    
    # Return as a named list
    hyp_goric <- as.list(hyp_goric)
    
    # average vcov matrix with its transpose to make it perfectly symmetric
    var_est_sym <- (var_est + t(var_est)) / 2
    
    # in case vcov matrix is not pd, add small ridge term to diagonal elements
    suppressWarnings({
      if(!matrixcalc::is.positive.definite(var_est_sym)) {
        var_est_pd <- var_est + diag(1e-10, nrow(var_est))
        goric_res <- restriktor::goric(object = est, VCOV = as.matrix(var_est_pd), hypotheses = hyp_goric, comparison = "complement")
      } else {
        goric_res <- restriktor::goric(object = est, VCOV = var_est, hypotheses = hyp_goric, comparison = "complement")
      }
    })
  } else {
    goric_res <- list()
  }

  
  # if there are exactly two hypotheses, return BF_12
  if(length(hypothesis) == 1){
    result <- list(bf_c = bf_res[["fit"]][["BF.c"]][1],
                   PMPc = bf_res[["fit"]][["PMPc"]][1],
                   bf12 = NA,
                   goric12 = NA,
                   goric_c = goric_res[["ratio.gw"]][1,2],
                   goric_weight = goric_res[["result"]][["gorica.weights_without_unc"]][1],
                   simplified = simplified) 

  } else if(length(hypothesis) == 2){
    result <- list(bf_c  = bf_res[["fit"]][["BF.c"]][1],
                   PMPc  = bf_res[["fit"]][["PMPc"]][1],
                   bf12  = bf_res[["BFmatrix"]][1,2],
                   goric12 = goric_res[["ratio.gw"]][1,2],
                   goric_c = NA,
                   goric_weight = goric_res[["result"]][["gorica.weights_without_unc"]][1],
                   simplified = simplified) 
    
  } else {
    result <- list(bf_c = bf_res[["fit"]][["BF.c"]][1],
                   PMPc = bf_res[["fit"]][["PMPc"]][1],
                   bf12 = NA,
                   goric12 = NA,
                   goric_c = NA,
                   goric_weight = goric_res[["result"]][["gorica.weights_without_unc"]][1],
                   simplified = simplified) 
  }
  
  return(result)
}


