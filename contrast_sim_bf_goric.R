################################################################################
## SIMULATION contrasts - BF and GORIC #########################################
################################################################################

library(ggplot2)

# BF/GORIC SIMULATION ##########################################################

# extract estimates
t.points <- c(0,2,6)
hypothesis <- "t3_diff>0"

betas <- c(0, 0.2, 0.4, 0, .5, .8)

var.u0 <- .1
var.u1 <- .1
var.e <- .5
cov <- .01
seed <- 123
fraction <- 1
attrition <- F

N_max <- 1000

n_iter <- 500
n_iter_outer <- 40


bf_c <- rep(NA, n_iter)
gor_c <- rep(NA, n_iter)
m_bf_c <- rep(NA, n_iter_outer)
m_gor_c <- rep(NA, n_iter_outer)
med_bf_c <- rep(NA, n_iter_outer)
med_gor_c <- rep(NA, n_iter_outer)
med_bf <- rep(NA, n_iter_outer)
med_gor <- rep(NA, n_iter_outer)

Ns <- round(seq(30, N_max, length.out=n_iter_outer))
Ns[Ns%%2==1] <- Ns[Ns%%2==1]+1 # make all N even numbers


res_bf <- list()

for(j in 1:n_iter_outer){
  for(i in 1:n_iter){
    suppressMessages({
      res_bf[[i]] <- get_bf_contrast(N=Ns[j], 
                                     t.points=t.points, 
                                     hypothesis=hypothesis, 
                                     betas=betas, 
                                     var.u0=var.u0, 
                                     var.u1=var.u1, 
                                     var.e=var.e, 
                                     cov=cov, 
                                     attrition=attrition, 
                                     fraction=1)
    })
    bf_c[i] <- res_bf[[i]]$bf_c
    gor_c[i] <- res_bf[[i]]$goric_c
  }
  m_bf_c[j] <- mean(log(bf_c), na.rm=T)
  m_gor_c[j] <- mean(log(gor_c), na.rm=T)
  med_bf_c[j] <- median(bf_c, na.rm=T)
  med_gor_c[j] <- median(gor_c, na.rm=T)
  print(j/n_iter_outer)
}

sim_data <- data.frame(
  size = c(m_bf_c, m_gor_c),
  N = rep(Ns, 2),
  type = rep(c("bf", "goric"), each = length(Ns))
)

pdf(file="./Plots/bf_vs_gor.pdf")

ggplot(data=sim_data, aes(x=N, y=size, color=type)) +
  geom_line() +
  theme(legend.position = c(.9, .11)) +
  ylab("mean(log)") +
  labs(title="H: -.5 < x < .5  vs.  Hc")

dev.off()

# sim_data2 <- data.frame(
#   size = c(med_bf_c, med_gor_c),
#   N = rep(Ns, 2),
#   type = rep(c("bf", "goric"), each = length(Ns))
# )
# 
# ggplot(data=sim_data2, aes(x=N, y=size, color=type)) +
#   geom_line() +
#   theme(legend.position = c(.1, .2)) +
#   ylim(c(0, 1000))


# POWER SIMULATION 1 ###########################################################



hyp1 <- "immediate_0 < immediate_1" # equivalent to "immediate_diff > 0"
hyp2 <- "retention_0 < retention_1" # equivalent to "retention_diff > 0"
hyp3 <- "overall_0 < overall_1"     # equivalent to "overall_diff > 0"
hyp4 <- "avg_diff > 0"
hyp5 <- "0.2 < t3_diff < 0.6" 
hyp6 <- "t2_diff > 0" 
hyp7 <- "0 < t2_diff < t3_diff"
# hyp8 <- "immediate_1 > retention_1"


# Mirjams suggestion for hypotheses
hypo1 <- "t1_diff>0"



att <- list(absent=FALSE, present=c(1, .7, .6))

bet <- c(0, 0.2, 0.4, 0, .5, .8)

n_iter <- 100
N_max <- 250

Ns <- round(seq(30, N_max, length.out=n_iter))
Ns[Ns%%2==1] <- Ns[Ns%%2==1]+1 # make all N even numbers

res <- list()
bf_pow <- list()
gor_pow <- list()
prop_simple <- list()

pb <- txtProgressBar(min = 0, max = n_iter, style = 3)

for(j in 1:2){
  k <- 0
  for(i in 1:n_iter){
    k <- k+1
    suppressMessages({
      res[[i]] <- get_power_contrast(N=Ns[i], 
                                     hypothesis = hyp5, 
                                     betas = bet, 
                                     m = 10000, 
                                     threshold = 5, 
                                     attrition = att[[j]],
                                     seed=123)
    })
    setTxtProgressBar(pb, k)
  }
  bf_pow[[j]] <- unlist(lapply(res, "[[", 1))
  gor_pow[[j]] <- unlist(lapply(res, "[[", 5))
  prop_simple[[j]] <- unlist(lapply(res, "[[", 7))
  print(j)
}

dat5 <- data.frame(
  power = unlist(c(bf_pow, gor_pow)),
  N = rep(Ns, 4),
  type = rep(c("bf", "goric"), each = length(Ns)*2),
  attrition = rep(c("absent", "present"), each=length(Ns))
)

# saveRDS(dat5, file="simdat_hyp5")

# plot results
p1 <- ggplot(data=dat1, aes(x=N, y=power, color=type, linetype = attrition)) +
        geom_line() +
        theme(legend.position = c(.9, .2)) +
        labs(title=paste("hypothesis: ", hyp1))

p2 <- ggplot(data=dat2, aes(x=N, y=power, color=type, linetype = attrition)) +
  geom_line() +
  theme(legend.position = "none") +
  labs(title=paste("hypothesis: ", hyp2))

p3 <- ggplot(data=dat3, aes(x=N, y=power, color=type, linetype = attrition)) +
  geom_line() +
  theme(legend.position = "none") +
  labs(title=paste("hypothesis: ", hyp3))

p4 <- ggplot(data=dat4, aes(x=N, y=power, color=type, linetype = attrition)) +
  geom_line() +
  theme(legend.position = "none") +
  labs(title=paste("hypothesis: ", hyp4))

p5 <- ggplot(data=dat5, aes(x=N, y=power, color=type, linetype = attrition)) +
  geom_line() +
  theme(legend.position = "none") +
  labs(title=paste("hypothesis: ", hyp5))

p6 <- ggplot(data=dat6, aes(x=N, y=power, color=type)) +
  geom_line() +
  theme(legend.position = "none") +
  labs(title=paste("hypothesis: ", hyp6))

p7 <- ggplot(data=dat7, aes(x=N, y=power, color=type)) +
  geom_line() +
  theme(legend.position = "none") +
  labs(title=paste("hypothesis: ", hyp7))

p8 <- ggplot(data=dat8, aes(x=N, y=power, color=type)) +
  geom_line() +
  theme(legend.position = "none") +
  labs(title=paste("hypothesis: ", hyp8))



pdf("p1_contrast_sim.pdf", width=8, height=6)
p1
dev.off()

library(gridExtra)
p <- grid.arrange(p1, p2, p3, p4, p5, p6, p7, p8)



# POWER SIMULATION 2 ###########################################################

h1 <- "overall_0 < overall_1"
h2 <- "overall_0 = overall_1"
# h3 <- "overall_0 > overall_1"
  
bet <- c(0, 0.2, 0.4, 0, .5, .8)

n_iter <- 100
N_max <- 200

Ns <- round(seq(10, N_max, length.out=n_iter))
Ns[Ns%%2==1] <- Ns[Ns%%2==1]+1 # make all N even numbers

res <- list()

pb <- txtProgressBar(min = 0, max = n_iter, style = 3)
k <- 0

for(i in 1:n_iter){
  k <- k+1
  suppressMessages({
    res[[i]] <- get_power_contrast(N=Ns[i], hypothesis = list(h1, h2), betas = bet, m=1000, threshold = .8, var.e = .001)
  })
  setTxtProgressBar(pb, k)
}

bf_power123 <- unlist(lapply(res, "[[", 2))
goric_power123 <- unlist(lapply(res, "[[", 6))


dat2.2 <- data.frame(
  power = c(bf_power123, goric_power123),
  N = rep(Ns, 2),
  type = rep(c("bf", "goric"), each = length(Ns))
)


p2.2 <- ggplot(data=dat2.2, aes(x=N, y=power, color=type)) +
  geom_line() +
  theme(legend.position = c(.9, .11)) +
  labs(title=paste("hypothesis: ", hyp1))
