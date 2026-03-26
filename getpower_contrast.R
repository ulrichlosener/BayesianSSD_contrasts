
get_power_contrast <- function(m=1000, 
                               N=100, 
                               t.points=c(0,2,6), 
                               hypothesis="immediate_0 < immediate_1",
                               betas=c(0, 0.1, 0.5, 0, .5, 1.5),
                               var.u0=.1,
                               var.u1=.1,
                               var.e=.2,
                               cov=.01,
                               fraction=1,
                               seed=NULL,
                               attrition=FALSE,
                               threshold=10){
  
  
  future::plan(future::sequential)  # Reset plan to avoid unexpected leftover parallel behavior
  future::plan(future::multisession, workers = future::availableCores() - 1)  # Use all but one core
  
  Ns <- rep(N, m)  # object to use lapply on, with first argument for the function (N)
  
  # Run simulation m times
  res <- future.apply::future_lapply(
    Ns,
    function(ss){
      get_bf_contrast(ss,
                       hypothesis=hypothesis, t.points=t.points,
                       var.u0=var.u0, var.u1=var.u1, cov=cov, var.e=var.e,
                       betas=betas, fraction=fraction, attrition=attrition)
    },
    future.seed = TRUE
  )
  
  future::plan(future::sequential) # Reset plan to avoid unexpected parallel behavior later
  
  bfs <- lapply(res, "[[", 1)
  pmps <- lapply(res, "[[", 2)
  bf12 <- lapply(res, "[[", 3)
  goric12 <- lapply(res, "[[", 4)
  goric_c <- lapply(res, "[[", 5)
  goric_weight <- lapply(res, "[[", 6)
  simple <- lapply(res, "[[", 7)
  
  power_bf <- sum(unlist(bfs) > threshold, na.rm = T)/m
  power_pmps <- sum(unlist(pmps) > threshold, na.rm = T)/m
  power_bf12 <- sum(unlist(bf12) > threshold, na.rm = T)/m
  power_goric12 <- sum(unlist(goric12) > threshold, na.rm = T)/m
  power_goric_c <- sum(unlist(goric_c) > threshold, na.rm = T)/m
  power_goric_weight <- sum(unlist(goric_weight) > threshold, na.rm = T)/m
  prop_simplified <- sum(unlist(simple))/m
  
  return(list(power_bf = power_bf,
              power_pmps = power_pmps,
              power_bf12 = power_bf12,
              power_goric12 = power_goric12,
              power_goric_c = power_goric_c,
              power_goric_weight = power_goric_weight,
              prop_simplified = prop_simplified))

}
