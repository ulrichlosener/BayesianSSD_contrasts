################################################################################
## SIMULATION contrasts - BF and GORIC #########################################
################################################################################

library(ggplot2)
library(gridExtra)
library(ggpubr)
### SECOND EDITION: all hyps in the table

# already here: 
# initial_diff - hyp1
# maintenance_diff - hyp2
# overall_diff - hyp3
# average_diff - hyp4
# t3_diff - hyp6

# missing:
# initial_0
# initial_1
# maintenance_0
# maintenance_1
# overall_0
# overall_1
# t1_diff
# t2_diff

hyp8 <- "initial_0 > 0"
hyp9 <- "initial_1 > 0"
hyp10 <- "maintenance_0 > 0"
hyp11 <- "maintenance_1 > 0"
hyp12 <- "overall_0 > 0"
hyp13 <- "overall_1 > 0"
hyp14 <- "t1_diff > 0"
hyp15 <- "t2_diff > 0"

# POWER SIMULATION 1 ###########################################################
hyp1 <- "initial_diff > 0" # equivalent to "initial_0 < initial_1"
hyp2 <- "maintenance_diff > 0" # equivalent to "maintenance_0 < maintenance_1"
hyp3 <- "overall_0 < overall_1"     # equivalent to "overall_diff > 0"
hyp4 <- "avg_diff > 0"
hyp5 <- "0.2 < t3_diff < 0.6" 
hyp6 <- "t3_diff > 0" 
hyp7 <- "0 < t2_diff < t3_diff"

att <- list(absent=FALSE, present=c(1, .7, .6))

bet <- c(0, 0.2, 0.4, 0, .5, .8)

n_iter <- 25
N_max <- 250

Ns <- round(seq(30, N_max, length.out=n_iter))
Ns[Ns%%2==1] <- Ns[Ns%%2==1]+1 # make all N even numbers

# Ns <- unique(dat12$N)

res <- list()
bf_pow <- list()
gor_pow <- list()
prop_simple <- list()

pb <- txtProgressBar(min = 0, max = n_iter, style = 3)
set.seed(123)
for(j in 1:2){
  k <- 0
  for(i in 1:n_iter){
    k <- k+1
    suppressMessages({
      res[[i]] <- get_power_contrast(N=Ns[i], 
                                     hypothesis = hyp13, 
                                     betas = bet, 
                                     m = 10000, 
                                     threshold = 5, 
                                     attrition = att[[j]],
                                     seed=123,
                                     method="both")
    })
    setTxtProgressBar(pb, k)
  }
  bf_pow[[j]] <- unlist(lapply(res, "[[", 1))
  gor_pow[[j]] <- unlist(lapply(res, "[[", 5))
  prop_simple[[j]] <- unlist(lapply(res, "[[", 7))
  print(j)
}

dat13 <- data.frame(
  power = unlist(c(bf_pow, gor_pow)),
  N = rep(Ns, 4),
  method = rep(c("BF", "GORICA"), each = length(Ns)*2),
  attrition = rep(rep(c("absent", "present"), each=length(Ns)),2)
)

saveRDS(dat13, file="simdat_hyp13.R")

dat1 <- readRDS("simdat_hyp1")
dat1 <- dat1[seq(1, nrow(dat1), by=2),]
dat2 <- readRDS("simdat_hyp2")
dat2 <- dat2[seq(1, nrow(dat2), by=2),]
dat3 <- readRDS("simdat_hyp3")
dat3 <- dat3[seq(1, nrow(dat3), by=2),]
dat4 <- readRDS("simdat_hyp4")
dat4 <- dat4[seq(1, nrow(dat4), by=2),]
dat5 <- readRDS("simdat_hyp5")
dat5 <- dat5[seq(1, nrow(dat5), by=2),]
dat6 <- readRDS("simdat_hyp6")
dat6 <- dat6[seq(1, nrow(dat6), by=2),]
dat7 <- readRDS("simdat_hyp7")
dat7 <- dat7[seq(1, nrow(dat7), by=2),]

dat8 <- readRDS("simdat_hyp8")
dat8 <- dat8[seq(1, nrow(dat8), by=2),]
dat9 <- readRDS("simdat_hyp9")
dat9 <- dat9[seq(1, nrow(dat9), by=2),]
dat10 <- readRDS("simdat_hyp10")
dat10 <- dat10[seq(1, nrow(dat10), by=2),]
dat11 <- readRDS("simdat_hyp11")
dat11 <- dat11[seq(1, nrow(dat11), by=2),]
dat12 <- readRDS("simdat_hyp12.R")
dat12 <- dat12[seq(1, nrow(dat12), by=2),]
dat13 <- readRDS("simdat_hyp13")
dat13 <- dat13[seq(1, nrow(dat13), by=2),]
dat14 <- readRDS("simdat_hyp14")
dat14 <- dat14[seq(1, nrow(dat14), by=2),]
dat15 <- readRDS("simdat_hyp15")
dat15 <- dat15[seq(1, nrow(dat15), by=2),]

# rename columns for new data frames
names(dat8)[names(dat8) == "method"] <- "type"   
names(dat9)[names(dat9) == "method"] <- "type"   
names(dat10)[names(dat10) == "method"] <- "type"   
names(dat11)[names(dat11) == "method"] <- "type"   
names(dat12)[names(dat12) == "method"] <- "type"   
names(dat13)[names(dat13) == "method"] <- "type"   
names(dat14)[names(dat14) == "method"] <- "type"   
names(dat15)[names(dat15) == "method"] <- "type"   

dat8[dat8 == "bf"] <- "BF"
dat8[dat8 == "goric"] <- "GORICA"


# plot results
p1 <- ggplot(data=dat1, aes(x=N, y=power, color=type, linetype = attrition)) +
  geom_line() +
  labs(title=paste(hyp1)) +
  ylim(0,1) +
  theme_bw()

p2 <- ggplot(data=dat2, aes(x=N, y=power, color=type, linetype = attrition)) +
  geom_line(show.legend=F) +
  labs(title=paste(hyp2)) +
  ylim(0,1) +
  theme_bw()

p3 <- ggplot(data=dat3, aes(x=N, y=power, color=type, linetype = attrition)) +
  geom_line(show.legend=F) +
  labs(title=paste(hyp3)) +
  ylim(0,1) +
  theme_bw()

p4 <- ggplot(data=dat4, aes(x=N, y=power, color=type, linetype = attrition)) +
  geom_line(show.legend=F) +
  labs(title=paste(hyp4)) +
  ylim(0,1) +
  theme_bw()

p5 <- ggplot(data=dat5, aes(x=N, y=power, color=type, linetype = attrition)) +
  geom_line(show.legend=F) +
  labs(title=paste(hyp5)) +
  ylim(0,1) +
  theme_bw()

p6 <- ggplot(data=dat6, aes(x=N, y=power, color=type, linetype = attrition)) +
  geom_line(show.legend=F) +
  labs(title=paste(hyp6)) +
  ylim(0,1) +
  theme_bw()

p7 <- ggplot(data=dat7, aes(x=N, y=power, color=type, linetype = attrition)) +
  geom_line(show.legend=F) +
  labs(title=paste(hyp7)) +
  ylim(0,1) +
  theme_bw()


# pdf("contrast_grid.pdf")
  ggarrange(p1, p2, p3, p4, p6, p7, 
          common.legend = TRUE, 
          legend = "bottom", 
          nrow = 3,
          ncol = 2)
# dev.off()

library(latex2exp)
  
p8 <- ggplot(data=dat8, aes(x=N, y=power, color=method, linetype = attrition)) +
  geom_line(linewidth = 1) +
  labs(title=expression("initial_C > 0 vs. H"[c]),
       subtitle=TeX((r'($\Delta \beta = 0.2$)'))) +
  ylim(0,1) +
  theme_bw()

p9 <- ggplot(data=dat9, aes(x=N, y=power, color=method, linetype = attrition)) +
  geom_line(linewidth = 1) +
  labs(title=expression("initial_E > 0 vs. H"[c]),
       subtitle=TeX((r'($\Delta \beta = 0.5$)'))) +
  ylim(0,1) +
  theme_bw()

p10 <- ggplot(data=dat10, aes(x=N, y=power, color=method, linetype = attrition)) +
  geom_line(linewidth = 1) +
  labs(title=expression("maintenance_C > 0 vs. H"[c]),
       subtitle=TeX((r'($\Delta \beta = 0.2$)'))) +
  ylim(0,1) +
  theme_bw()

p11 <- ggplot(data=dat11, aes(x=N, y=power, color=method, linetype = attrition)) +
  geom_line(linewidth = 1) +
  labs(title=expression("maintenance_E > 0 vs. H"[c]),
       subtitle=TeX((r'($\Delta \beta = 0.3$)'))) +
  ylim(0,1) +
  theme_bw()

p12 <- ggplot(data=dat12, aes(x=N, y=power, color=method, linetype = attrition)) +
  geom_line(linewidth = 1) +
  labs(title=expression("overall_C > 0 vs. H"[c]),
       subtitle=TeX((r'($\Delta \beta = 0.4$)'))) +
  ylim(0,1) +
  theme_bw()

p13 <- ggplot(data=dat13, aes(x=N, y=power, color=method, linetype = attrition)) +
  geom_line(linewidth = 1) +
  labs(title=expression("overall_E > 0 vs. H"[c]),
       subtitle=TeX((r'($\Delta \beta = 0.8$)'))) +
  ylim(0,1) +
  theme_bw()

p14 <- ggplot(data=dat14, aes(x=N, y=power, color=type, linetype = attrition)) +
  geom_line(linewidth = 1) +
  labs(title=paste(hyp14)) +
  ylim(0,1) +
  theme_bw()

p15 <- ggplot(data=dat15, aes(x=N, y=power, color=method, linetype = attrition)) +
  geom_line(linewidth = 1) +
  labs(title=paste(hyp15)) +
  ylim(0,1) +
  theme_bw()

# pdf("contrast_grid_new.pdf")
ggarrange(p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11, p12, p13, p14, p15,
          common.legend = TRUE, 
          legend = "bottom", 
          nrow = 4,
          ncol = 4)
# dev.off()

pdf("contrast_grid_1-6.pdf", width = 9)
  ggarrange(p8, p9, p10, p11, p12, p13,
            common.legend = TRUE, 
            legend = "bottom", 
            nrow = 2,
            ncol = 3)
dev.off()

# pdf("contrast_grid.pdf", height=8, width=10)
#   grid.arrange(p1, p2, p3, p4, p6, p7)
# dev.off()

# POWER SIMULATION 2 ###########################################################


hyp_pair1 <- list("initial_0 < initial_1", "initial_0 = initial_1")#
hyp_pair2 <- list("maintenance_0 < maintenance_1", "maintenance_0 = maintenance_1")#
hyp_pair3 <- list("overall_0 < overall_1", "overall_0 = overall_1")#
hyp_pair4 <- list("avg_diff > 0", "avg_diff = 0")#
# hyp_pair5 <- list("t3_diff = 0.4", "0.2 < t3_diff < 0.6")
hyp_pair6 <- list("t3_diff > 0", "t3_diff = 0")#
hyp_pair7 <- list("0 < t2_diff < t3_diff", "0 = t2_diff = t3_diff")#

hyp_pair8 <- list("initial_0 > 0", "initial_0 = 0")
hyp_pair9 <- list("initial_1 > 0", "initial_1 = 0")
hyp_pair10 <- list("maintenance_0 > 0", "maintenance_0 = 0")
hyp_pair11 <- list("maintenance_1 > 0", "maintenance_1 = 0")
hyp_pair12 <- list("overall_0 > 0", "overall_0 = 0")
hyp_pair13 <- list("overall_1 > 0","overall_1 = 0")
hyp_pair14 <- "-0.2 < t1_diff < 0.2"
hyp_pair15 <- list("t2_diff > 0", "t2_diff = 0")

bet <- c(0, 0.2, 0.4, 0, .5, .8)

n_iter <- 25
N_max <- 250

Ns <- round(seq(30, N_max, length.out=n_iter))
Ns[Ns%%2==1] <- Ns[Ns%%2==1]+1 # make all N even numbers
att <- list(absent=FALSE, present=c(1, .7, .6))

res <- list()
bf_pow12 <- list()
gor_pow12 <- list()
prop_simple12 <- list()

pb <- txtProgressBar(min = 0, max = n_iter, style = 3)
k <- 0

for(j in 1:2){
  k <- 0
  for(i in 1:n_iter){
    k <- k+1
    suppressMessages({
      res[[i]] <- get_power_contrast(N = Ns[i], 
                                     hypothesis = hyp_pair12, 
                                     betas = bet, 
                                     m = 1000, 
                                     threshold = 5,
                                     attrition = att[[j]],
                                     var.u0 = .1,
                                     var.u1 = .1,
                                     var.e = .2,
                                     cov = .01,
                                     t.points = c(0,2,6), 
                                     method = "both",
                                     seed = 123)
    })
    setTxtProgressBar(pb, k)
  }
  bf_pow12[[j]] <- unlist(lapply(res, "[[", 3))
  gor_pow12[[j]] <- unlist(lapply(res, "[[", 4))
  prop_simple12[[j]] <- unlist(lapply(res, "[[", 7))
  print(j)
}  

dat2.12 <- data.frame(
  power = unlist(c(bf_pow12, gor_pow12)),
  N = rep(Ns, 4),
  method = rep(c("bf", "goric"), each = length(Ns)*2),
  attrition = rep(c("absent", "present"), each=length(Ns))
)

# saveRDS(dat2.12, file="simdat2.15_small.R")

dat2.1 <- readRDS("simdat2.1.R")
dat2.1 <- dat2.1[seq(1, nrow(dat2.1), by=2),]
dat2.2 <- readRDS("simdat2.2.R")
dat2.2 <- dat2.2[seq(1, nrow(dat2.2), by=2),]
dat2.3 <- readRDS("simdat2.3")
dat2.3 <- dat2.3[seq(1, nrow(dat2.3), by=2),]
dat2.4 <- readRDS("simdat2.4")
dat2.4 <- dat2.4[seq(1, nrow(dat2.4), by=2),]
dat2.6 <- readRDS("simdat2.6")
dat2.6 <- dat2.6[seq(1, nrow(dat2.6), by=2),]
dat2.7 <- readRDS("simdat2.7.R")
dat2.7 <- dat2.7[seq(1, nrow(dat2.7), by=2),]


dat2.8 <- readRDS("simdat2.8.R")
dat2.8 <- dat2.8[seq(1, nrow(dat2.8), by=2),]
dat2.9 <- readRDS("simdat2.9.R")
dat2.9 <- dat2.8[seq(1, nrow(dat2.9), by=2),]
dat2.10 <- readRDS("simdat2.10.R")
dat2.10 <- dat2.10[seq(1, nrow(dat2.10), by=2),]
dat2.11 <- readRDS("simdat2.11.R")
dat2.11 <- dat2.11[seq(1, nrow(dat2.11), by=2),]
dat2.12 <- readRDS("simdat2.12.R")
dat2.12 <- dat2.12[seq(1, nrow(dat2.12), by=2),]
dat2.13 <- readRDS("simdat2.13.R")
dat2.13 <- dat2.13[seq(1, nrow(dat2.13), by=2),]
dat2.15 <- readRDS("simdat2.15.R")
dat2.15 <- dat2.15[seq(1, nrow(dat2.15), by=2),]


p2.1 <- ggplot(data=dat2.1, aes(x=N, y=power, color=type, linetype = attrition)) +
  geom_line(linewidth = 1) +
  theme_bw() +
  labs(title="initial_diff > 0 vs. initial_diff = 0") +
  ylim(0,1)

p2.2 <- ggplot(data=dat2.2, aes(x=N, y=power, color=type, linetype = attrition)) +
  geom_line(linewidth = 1) +
  theme_bw() +
  labs(title="maintenance_diff > 0 vs. maintenance_diff = 0") +
  ylim(0,1)

p2.3 <- ggplot(data=dat2.3, aes(x=N, y=power, color=type, linetype = attrition)) +
  geom_line(linewidth = 1) +
  theme_bw() +
  labs(title="overall_diff > 0 vs. overall_diff = 0") +
  ylim(0,1)

p2.4 <- ggplot(data=dat2.4, aes(x=N, y=power, color=type, linetype = attrition)) +
  geom_line(linewidth = 1) +
  theme_bw() +
  labs(title="avg_diff > 0 vs. avg_diff = 0") +
  ylim(0,1)

p2.6 <- ggplot(data=dat2.6, aes(x=N, y=power, color=type, linetype = attrition)) +
  geom_line(linewidth = 1) +
  theme_bw() +
  labs(title="t3_diff > 0 vs. t3_diff = 0") +
  ylim(0,1)

p2.7 <- ggplot(data=dat2.7, aes(x=N, y=power, color=type, linetype = attrition)) +
  geom_line(linewidth = 1) +
  theme_bw() +
  labs(title="0 < t2_diff < t3_diff vs. 0 = t2_diff = t3_diff") +
  ylim(0,1)


# pdf("contrast_grid2.pdf")
ggarrange(p2.1, p2.2, p2.3, p2.4, p2.6, p2.7, 
          common.legend = TRUE, 
          legend = "bottom", 
          nrow = 3,
          ncol = 2)
# dev.off()

p2.8 <- ggplot(data=dat2.8, aes(x=N, y=power, color=method, linetype = attrition)) +
  geom_line(linewidth = 1) +
  theme_bw() +
  labs(title=expression("initial_C > 0 vs. H"[0]),
       subtitle=TeX((r'($\Delta \beta = 0.2$)'))) +
  ylim(0,1)

p2.9 <- ggplot(data=dat2.9, aes(x=N, y=power, color=method, linetype = attrition)) +
  geom_line(linewidth = 1) +
  theme_bw() +
  labs(title=expression("initial_E > 0 vs. H"[0]),
       subtitle=TeX((r'($\Delta \beta = 0.5$)'))) +
  ylim(0,1)

p2.10 <- ggplot(data=dat2.10, aes(x=N, y=power, color=method, linetype = attrition)) +
  geom_line(linewidth = 1) +
  theme_bw() +
  labs(title=expression("maintenance_C > 0 vs. H"[0]),
       subtitle=TeX((r'($\Delta \beta = 0.2$)'))) +
  ylim(0,.8)

p2.11 <- ggplot(data=dat2.11, aes(x=N, y=power, color=method, linetype = attrition)) +
  geom_line(linewidth = 1) +
  theme_bw() +
  labs(title=expression("maintenance_E > 0 vs. H"[0]),
       subtitle=TeX((r'($\Delta \beta = 0.3$)'))) +
  ylim(0,1)

p2.12 <- ggplot(data=dat2.12, aes(x=N, y=power, color=method, linetype = attrition)) +
  geom_line(linewidth = 1) +
  theme_bw() +
  labs(title=expression("overall_C > 0 vs. H"[0]),
       subtitle=TeX((r'($\Delta \beta = 0.4$)'))) +
  ylim(0,1)

p2.13 <- ggplot(data=dat2.13, aes(x=N, y=power, color=method, linetype = attrition)) +
  geom_line(linewidth = 1) +
  theme_bw() +
  labs(title=expression("overall_E > 0 vs. H"[0]),
       subtitle=TeX((r'($\Delta \beta = 0.8$)'))) +
  ylim(0,1)

p2.15 <- ggplot(data=dat2.15, aes(x=N, y=power, color=method, linetype = attrition)) +
  geom_line(linewidth = 1) +
  theme_bw() +
  labs(title=hyp_pair15) +
  ylim(0,1)


dat2.8[dat2.8 == "bf"] <- "BF"
dat2.8[dat2.8 == "goric"] <- "GORICA"
 

pdf("contrast_grid2_1-6.pdf", width = 9)
ggarrange(p2.8, p2.9, p2.10, p2.11, p2.12, p2.13,
          common.legend = TRUE, 
          legend = "bottom", 
          nrow = 2,
          ncol = 3)
dev.off()


# b fraction SIMULATION ########################################################

hyp_pair_b <- list("initial_0 < initial_1", "initial_0 = initial_1")
hyp_pair_b <- list("initial_0 > 0", "initial_0 = 0") # alternative
  
bet <- c(0, 0.2, 0.4, 0, .5, .8)

n_iter <- 25
N_max <- 250

Ns <- round(seq(30, N_max, length.out=n_iter))
Ns[Ns%%2==1] <- Ns[Ns%%2==1]+1 # make all N even numbers
att <- list(absent=FALSE, present=c(1, .7, .6))

res <- list()
bf_pow12 <- list()
gor_pow12 <- list()
prop_simple12 <- list()

pb <- txtProgressBar(min = 0, max = n_iter, style = 3)
k <- 0

for(j in 1:2){
  k <- 0
  for(i in 1:n_iter){
    k <- k+1
    suppressMessages({
      res[[i]] <- get_power_contrast(N = Ns[i], 
                                     hypothesis = hyp_pair_b, 
                                     betas = bet, 
                                     m = 10000, 
                                     threshold = 5,
                                     attrition = att[[j]],
                                     var.u0 = .1,
                                     var.u1 = .1,
                                     var.e = .2,
                                     cov = .01,
                                     t.points = c(0,2,6), 
                                     method = "bf",
                                     seed = 123,
                                     fraction = 2)
    })
    setTxtProgressBar(pb, k)
  }
  bf_pow12[[j]] <- unlist(lapply(res, "[[", 3))
  gor_pow12[[j]] <- unlist(lapply(res, "[[", 4))
  prop_simple12[[j]] <- unlist(lapply(res, "[[", 7))
  print(j)
}  

dat_b2 <- data.frame(
  power = unlist(c(bf_pow12, gor_pow12)),
  N = rep(Ns, 4),
  method = rep(c("bf", "goric"), each = length(Ns)*2),
  attrition = rep(c("absent", "present"), each=length(Ns))
)



dat_b1 <- dat_b[dat_b$group == "b1",]
dat_b1 <- dat_b1[seq(1, nrow(dat_b1), by=2),]

dat_gor <- dat_b[dat_b$group == "gorica", ]
dat_gor <- dat_gor[seq(1, nrow(dat_gor), by=2),]

dat_b2 <- dat_b2[1:50,]
dat_b3 <- dat_b3[1:50,]

dat_b <- rbind(dat_b1, dat_b2, dat_b3, dat_gor)

dat_b$group <- rep(c("b1", "b2", "b3", "gorica"), each=50)

library(dplyr)

dat_b <- bind_rows(
  dat_b1[1:100, ]     |> mutate(group = "b1"),
  dat_b2[1:100, ]     |> mutate(group = "b2"),
  dat_b3[1:100, ]     |> mutate(group = "b3"),
  dat_b1[101:200, ]   |> mutate(group = "gorica")
)

ggplot(dat_b, aes(x = N, y = power,
                      color = group,
                      linetype = attrition)) +
  geom_line() +
  scale_linetype_manual(
    values = c("absent" = "solid",
               "present" = "dashed")
  ) +
  theme_minimal()

# saveRDS(dat_b, file="dat_b.R")

dat_b <- readRDS("dat_b.R")



pdf("contrast_plot_b_v2.pdf", width=9)
ggplot(data=dat_b, aes(x=N, y=power, color=group, linetype = attrition)) +
  geom_line(linewidth = 1) +
  theme_bw() +
  labs(title=expression("initial_C > 0 vs. H"[0]),
       subtitle=TeX((r'($\Delta \beta = 0.2$)'))) +
  ylim(0,.8) +
  scale_color_manual(name = "method",
                     labels = c(expression("BF using b"[1]), 
                                expression("BF using b"[2]), 
                                expression("BF using b"[3]),
                                "GORICA"),
                     values = c("#F8766D", "#00BA38", "#619CFF", "#C77CFF"))
dev.off()






# BF/GORIC SIMULATION ##########################################################

# extract estimates
t.points <- c(0,2,6)
hypothesis <- "0 < t3_diff < 0.8"
hypothesis <- list("t3_diff > 0", "t3_diff=0")

betas <- c(0, 0.2, 0.4, 0, 0.5, 0.8)

var.u0 <- .1
var.u1 <- .1
var.e <- .5
cov <- .01
seed <- 123
fraction <- 1
attrition <- F

N_max <- 500

n_iter <- 1000
n_iter_outer <- 10


bf_c <- rep(NA, n_iter)
gor_c <- rep(NA, n_iter)
m_bf_c <- rep(NA, n_iter_outer)
m_gor_c <- rep(NA, n_iter_outer)
med_bf_c <- rep(NA, n_iter_outer)
med_gor_c <- rep(NA, n_iter_outer)
med_bf <- rep(NA, n_iter_outer)
med_gor <- rep(NA, n_iter_outer)
prop_bf <- rep(NA, n_iter_outer)
prop_gor <- rep(NA, n_iter_outer)

Ns <- round(seq(30, N_max, length.out=n_iter_outer))
Ns[Ns%%2==1] <- Ns[Ns%%2==1]+1 # make all N even numbers


res_bf <- list()

pb <- txtProgressBar(min = 0, max = n_iter_outer, style = 3)
k <- 0

for(j in 1:n_iter_outer){
  k <- k+1
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
    bf_c[i] <- res_bf[[i]]$bf12
    gor_c[i] <- res_bf[[i]]$goric12
  }
  m_bf_c[j] <- mean(log(bf_c), na.rm=T)
  m_gor_c[j] <- mean(log(gor_c), na.rm=T)
  med_bf_c[j] <- median(bf_c, na.rm=T)
  med_gor_c[j] <- median(gor_c, na.rm=T)
  prop_bf[j] <- sum(bf_c>10, na.rm=T)/length(bf_c)
  prop_gor[j] <- sum(gor_c>10, na.rm=T)/length(gor_c)
  setTxtProgressBar(pb, k)
}

sim_data <- data.frame(
  size = c(prop_bf, prop_gor),
  N = rep(Ns, 2),
  method = rep(c("bf", "goric"), each = length(Ns))
)

#pdf(file="./Plots/bf_vs_gor.pdf")

ggplot(data=sim_data, aes(x=N, y=size, color=method)) +
  geom_line() +
  theme(legend.position = c(.9, .11)) +
  ylab("mean(log)") +
  labs(title="H1: t3_diff > 0, H2: t3_diff=0")

dev.off()

# sim_data2 <- data.frame(
#   size = c(med_bf_c, med_gor_c),
#   N = rep(Ns, 2),
#   method = rep(c("bf", "goric"), each = length(Ns))
# )
# 
# ggplot(data=sim_data2, aes(x=N, y=size, color=method)) +
#   geom_line() +
#   theme(legend.position = c(.1, .2)) +
#   ylim(c(0, 1000))

