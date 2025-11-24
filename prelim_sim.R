################################################################################
# SIMULATION BAYES FACTOR
################################################################################

thres <- 5

n_iter_outer <- 10
n_iter_inner <- 1000
bf <- rep(NA, n_iter_inner)


################################################################################
# effect size


res_ef <- rep(NA, n_iter_outer)

diff_a <- rep(0, n_iter_outer)
diff_b <- seq(0, .1, length.out=n_iter_outer)
diff_c <- seq(0, 1, length.out=n_iter_outer)

for(j in 1:n_iter_outer){
  for(i in 1:n_iter_inner){
    bf[i] <- getbf_contrast(N=500,
                            betas=c(0, 0, 0, diff_a[j], diff_b[j], diff_c[j]))
  }
  res_ef[j] <- sum(bf>thres, na.rm=T)/n_iter_inner
  print(j/n_iter_outer)
}


plot(res_ef, type="l", xlab = "effect size", ylab="power")



################################################################################
# Sample size


res_N <- rep(NA, n_iter_outer)

Ns <- round(seq(30, 500, length.out=n_iter_outer))
Ns[which(Ns %% 2 !=0)] <- Ns[which(Ns %% 2 !=0)]+1

for(j in 1:n_iter_outer){
  for(i in 1:n_iter_inner){
    bf[i] <- getbf_contrast(N=Ns[j],
                            betas=c(0, 0, 0, 0, .1, 1))
  }
  res_N[j] <- sum(bf>thres, na.rm = T)/n_iter_inner
  print(j/n_iter_outer)
}


plot(y=res_N, x=Ns, type="l", xlab = "N", ylab="power")


################################################################################
# SIMULATION POWER
################################################################################











