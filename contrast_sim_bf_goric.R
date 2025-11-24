################################################################################
## SIMULATION contrasts - BF and GORIC
################################################################################
library(ggplot2)

hyp1 <- "immediate_0 < immediate_1" # true
hyp2 <- "t1_diff = 0" # true
hyp3 <- "immediate_diff < overall_diff" # true
hyp4 <- "immediate_0 = retention_0 < immediate_1" # true
hyp5 <- "avg_diff > 0"
hyp6 <- "immediate_diff < overall_diff < 1"

bet <- c(0, 0.2, 0.4, 0, .5, .8)

n_iter <- 100
N_max <- 200

Ns <- round(seq(10, N_max, length.out=n_iter))
Ns[Ns%%2==1] <- Ns[Ns%%2==1]+1 # make all N even numbers

res <- list()

for(i in 1:n_iter){
  res[[i]] <- get_power_contrast(N=Ns[i], hypothesis = hyp6, betas = bet, m=1000)
  print(i/n_iter)
}

# extract results
bf_power <- unlist(lapply(res, "[[", 1))
goric_power <- unlist(lapply(res, "[[", 2))

dat6 <- data.frame(
  power = c(bf_power, goric_power),
  N = rep(Ns, 2),
  type = rep(c("bf", "goric"), each = length(Ns))
)

# plot results
p1 <- ggplot(data=dat1, aes(x=N, y=power, color=type)) +
  geom_line() +
  theme(legend.position = c(.9, .11)) +
  labs(title=paste("hypothesis = ", hyp1))

p2 <- ggplot(data=dat2, aes(x=N, y=power, color=type)) +
  geom_line() +
  theme(legend.position = "none") +
  labs(title=paste("hypothesis = ", hyp2))

p3 <- ggplot(data=dat3, aes(x=N, y=power, color=type)) +
  geom_line() +
  theme(legend.position = "none") +
  labs(title=paste("hypothesis = ", hyp3))

p4 <- ggplot(data=dat4, aes(x=N, y=power, color=type)) +
  geom_line() +
  theme(legend.position = "none") +
  labs(title=paste("hypothesis = ", hyp4))

p5 <- ggplot(data=dat5, aes(x=N, y=power, color=type)) +
  geom_line() +
  theme(legend.position = "none") +
  labs(title=paste("hypothesis = ", hyp5))

p6 <- ggplot(data=dat6, aes(x=N, y=power, color=type)) +
  geom_line() +
  theme(legend.position = "none") +
  labs(title=paste("hypothesis = ", hyp6))

library(gridExtra)
p <- grid.arrange(p1, p2, p3, p4, p5, p6)


