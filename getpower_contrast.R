
get_power_contrast <- function(m=1000, 
                              N=100, 
                              t.points=c(0,1,5), 
                              hypothesis="immediate_diff < retention_diff < overall_diff",
                              betas=c(0, 0.1, 0.5, 0, .5, 1.5),
                              var.u0=.1,
                              var.u1=.1,
                              var.e=.2,
                              cov=.01,
                              fraction=1,
                              seed=NULL,
                              PMPthres=.9, 
                              BFthres=5,
                              GORICthres=5){
  
  
  future::plan(future::sequential)  # Reset plan to avoid unexpected leftover parallel behavior
  future::plan(future::multisession, workers = future::availableCores() - 1)  # Use all but one core
  
  Ns <- rep(N, m)  # object to use lapply on with first argument for the function (N)
  
  # Run simulation m times
  res <- future.apply::future_lapply(
    Ns,
    function(ss){
      get_bf_contrast(ss,
                       hypothesis=hypothesis, t.points=t.points,
                       var.u0=var.u0, var.u1=var.u1, cov=cov, var.e=var.e,
                       betas=betas, fraction=fraction, seed=seed)
    },
    future.seed = TRUE
  )
  
  future::plan(future::sequential) # Reset plan to avoid unexpected parallel behavior later

  
  bfs <- lapply(res, "[[", 1)
  gorics <- lapply(res, "[[", 2)
  
  
  
  power_bf <- mean(bfs > BFthres, na.rm = T)
  power_goric <- mean(gorics > GORICthres, na.rm = T)
  
  
  return(list(power_bf = power_bf,
              power_goric = power_goric))

}
