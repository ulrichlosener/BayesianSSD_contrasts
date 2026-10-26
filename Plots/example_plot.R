################################################################################
## EXAMPLE GRAPH PRE, POST, FOLLOW-UP
################################################################################


library(ggplot2)

t <- c(0, 2, 6)
y <- c(0, 0, 0.2, .5, .4, .8)

cond <- c("control", "experimental")

ex_dat <- data.frame(
  t = rep(t, each = 2),
  condition = factor(rep(cond, 3)),
  y = y
)

# Variable indicating which line segment each observation belongs to
ex_dat$segment <- ifelse(ex_dat$t < 2, "solid", "dotted")

# plot
ex_plot <- ggplot(ex_dat, aes(x = t, y = y,
                   color = condition,
                   shape = condition,
                   group = condition)) +
  # Solid  
  geom_line(
    data = subset(ex_dat, t <= 2),
    linewidth = 1
  ) +
  # Dotted  
  geom_line(
    data = subset(ex_dat, t >= 2),
    linewidth = 1,
    linetype = "dotted"
  ) +
  geom_point(size = 3) +
  # Value labels
  geom_text(
    aes(label = y),
    vjust = -0.8,
    show.legend = FALSE
  ) +
  ylim(0, 0.85) +
  labs(x = "Time", y = "Quality of Life") +
  theme_minimal() +
  theme(
    legend.position = c(.2, .75),
    legend.background = element_rect(
      linewidth = 0.5,
      linetype = "solid"
    )
  )

ggsave("example_plot.pdf", ex_plot, units="in", scale=1)
