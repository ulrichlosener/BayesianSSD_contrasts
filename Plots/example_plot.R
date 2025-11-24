################################################################################
## EXAMPLE GRAPH PRE, POST, FOLLOW-UP
################################################################################

library(ggplot2)

t <- c(0, 2, 6)
y <- c(10, 9.5, 5, 4, 3, 1)
cond <- 0:1

dat <- data.frame(t=rep(c(0,2,6), each=2), condition=as.factor(0:1), y=c(10, 9.8, 5, 4, 3.8, .8))

ex_plot <- ggplot(data=dat, aes(x=t, y=y, color=condition, shape=condition)) +
        geom_line(linewidth = 1) +
        geom_point(size=3) +
        ylim(c(0, 10)) +
        ylab("outcome") +
        xlab("time") +
        theme_minimal() +
        theme(legend.position=c(.75, .75),
              legend.background = element_rect(size=0.5, linetype="solid"))

ggsave("example_plot.pdf",ex_plot, units="in", scale=1)
