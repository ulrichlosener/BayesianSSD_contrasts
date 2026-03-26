get_neff_contrast <- function(model) {
  
  # --- Fixed design ---
  X <- model.matrix(model)
  n_obs <- nrow(X)
  
  # --- Random-effects components (correct way) ---
  Z      <- lme4::getME(model, "Z")
  Lambda <- lme4::getME(model, "Lambda")
  sigma2 <- sigma(model)^2
  
  # --- Marginal covariance ---
  V <- as.matrix(Z %*% Lambda %*% t(Lambda) %*% t(Z) +
                   sigma2 * diag(n_obs))
  
  W <- diag(diag(V))  # independence assumption
  
  # --- Information matrices ---
  info_dep   <- crossprod(X, solve(V, X))
  info_indep <- crossprod(X, solve(W, X))
  
  var_dep   <- solve(info_dep)
  var_indep <- solve(info_indep)
  
  w_vals <- diag(var_indep / var_dep)
  w <- mean(w_vals)
  
  n_eff <- w * n_obs
  
  return(n_eff)
}
