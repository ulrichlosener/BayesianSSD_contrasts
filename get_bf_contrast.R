# CONTRASTS
library(lme4)
library(MASS)
library(matrixcalc) # for is.singular

# !So far only works for 2 conditions and 3 timepoints!

get_bf_contrast <- function(N=100, 
                           t.points=c(0,1,5), 
                           hypothesis="immediate_diff < retention_diff < overall_diff",
                           betas=c(0, 0.4, 0.5, 0, 1, 3),
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
  
  # Random effects covariance matrix
  Sigma_u <- matrix(c(var.u0, cov,
                      cov, var.u1), 2, 2)
  
  cond <- 0:1
  n <- length(t.points)
  
  # make dummies: 6 means for 3 timepoints and 2 conditions
  
  # Create full design
  dat <- data.frame(id=rep(1:N, each=n), t=rep(t.points, N), cond=rep(cond, each=N*n/2))
  
  # Dummies for each time × condition combination
  dat$D1 <- as.numeric(dat$t == t.points[1] & dat$cond == 0)
  dat$D2 <- as.numeric(dat$t == t.points[2] & dat$cond == 0)
  dat$D3 <- as.numeric(dat$t == t.points[3] & dat$cond == 0)
  dat$D4 <- as.numeric(dat$t == t.points[1] & dat$cond == 1)
  dat$D5 <- as.numeric(dat$t == t.points[2] & dat$cond == 1)
  dat$D6 <- as.numeric(dat$t == t.points[3] & dat$cond == 1)
  
  # Random effects
  u <- MASS::mvrnorm(N, mu = c(0, 0), Sigma = Sigma_u)
  dat$u0 <- u[dat$id, 1]
  dat$u1 <- u[dat$id, 2]
  
  # Outcome
  dat$y <- with(dat,
                betas[1]*D1 + betas[2]*D2 + betas[3]*D3 +
                betas[4]*D4 + betas[5]*D5 + betas[6]*D6 +
                u0 + u1*t + rnorm(nrow(dat), 0, sqrt(var.e))
  )
  
  ################################################################################
  ## Fit MLM
  
  mod <- lme4::lmer(y ~ -1 + D1 + D2 + D3 + D4 + D5 + D6 + (1 + t|id), data = dat)
  beta_hat <- lme4::fixef(mod)
  var_beta_hat <- vcov(mod)
  
  ################################################################################
  ## Contrasts
  
  # create all contrasts of interest and name them
  c <- rbind(
    # effects within conditions
    "immediate_0" = c(-1, 1, 0, 0, 0, 0),
    "immediate_1" = c(0, 0, 0, -1, 1, 0),
    "retention_0" = c(0, -1, 1, 0, 0, 0),
    "retention_1" = c(0, 0, 0, 0, -1, 1),
    "overall_0"   = c(-1, 0, 1, 0, 0, 0),
    "overall_1"   = c(0, 0, 0, -1, 0, 1),
    # comparisons of means between conditions
    "avg_diff"  = c(-1/3, -1/3, -1/3, 1/3, 1/3, 1/3),
    "t1_diff"   = c(-1, 0, 0, 1, 0, 0),
    "t2_diff"   = c(0, -1, 0, 0, 1, 0),
    "t3_diff"   = c(0, 0, -1, 0, 0, 1),
    # comparisons of effects between conditions
    "immediate_diff" = c(0, 0, 0, -1, 1, 0) - c(-1, 1, 0, 0, 0, 0), # immediate_1 - immediate_0
    "retention_diff" = c(0, 0, 0, 0, -1, 1) - c(0, -1, 1, 0, 0, 0), # retention_1 - retention_0
    "overall_diff"   = c(0, 0, 0, -1, 0, 1) - c(-1, 0, 1, 0, 0, 0)  # overall_1   - overall_0
  )
  
  # extract contrasts present in hypothesis
  rel_c <- unlist(strsplit(hypothesis, "\\s*[<>=]\\s*"))
  # keep only elements containing letters (no numbers)
  relevant_c <- rel_c[grepl("[A-Za-z]", rel_c)]  
  
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
  
  bf_res <- bain::bain(x = est, Sigma = var_matrices, 
                 hypothesis = hypothesis, n = rep(N, length(est)), 
                 group_parameters = 1, joint_parameters = 0, fraction=fraction)
  
  # extract BFc
  bf_c <- bf_res[["fit"]][["BF.c"]][1]
  
  ##############################################################################
  ## GORIC(A)
  
  # average vcov matrix with its transpose to make it perfectly symmetric
  var_est_sym <- (var_est + t(var_est)) / 2
  
  # in case vcov matrix is not pd, add small ridge term to diagonal elements
  suppressMessages({
    if(!is.positive.definite(var_est_sym)) {
      var_est_pd <- var_est + diag(1e-10, nrow(var_est))
      goric_res <- goric(object = est, VCOV = as.matrix(var_est_pd), hypotheses = list(H1 = hypothesis), comparison = "complement")
    } else {
      goric_res <- goric(object = est, VCOV = var_est, hypotheses = list(H1 = hypothesis), comparison = "complement")
    }
  })
  
  goric_c <- goric_res[["ratio.gw"]][1,2]


  return(list(bf_c    = bf_c,
              goric_c = goric_c))
}







