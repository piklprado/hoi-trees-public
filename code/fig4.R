# Reproduce Fig. 4

# raw cummulative effects
cumm.eff.raw <-
  dat.org.scale %>%
  # turn observations into long format
  gather(Predictor, Value, -non.X.colnames) %>%
  # also do not take OTHERS' cumulative effects into account
  filter(!str_detect(Predictor, "OTHERS")) %>% 
  # remove zero neighbour DBH before joining params
  # this will avoid the cummulative effects to bias towards zero
  filter(Value != 0) %>% 
  # assign i, j, k to neighbour ID
  mutate(Param.type = ifelse(Predictor %in% neigh, "alpha", "beta"),
         Neigh1 = ifelse(word(Predictor, 1, sep = "_")==Focal, "i", "j"),
         Neigh2 = ifelse(word(Predictor, 2, sep = "_")==Focal, "i", 
                         ifelse(word(Predictor, 2, sep = "_")=="2", Neigh1,
                                ifelse(Neigh1=="i", "j", "k"))),
         Neigh12 = paste0("i", Neigh1, ifelse(is.na(Neigh2), "", Neigh2))) %>% 
  # join parameters
  left_join(beta.coef.backscale) %>%
  # multiply observed neighbour count to parameter estimates
  mutate(eff_median = Value * Estimate,
         eff_lci    = Value * Q2.5,
         eff_uci    = Value * Q97.5)

# summarise raw effects to individual-focal-neighbourhood level
cumm.eff.i.f.n <-
  cumm.eff.raw %>%
  group_by(Plot, Year, Tag.tree, Focal, Predictor, Param.type, Neigh12) %>%
  summarise_at(vars(eff_median, eff_lci, eff_uci), sum) %>% 
  ungroup()

# summarise to focal species level
cumm.eff.ijk <-
  cumm.eff.i.f.n %>%
  group_by(Focal, Param.type, Neigh12) %>%
  summarise(
    median = median(eff_median),
    lci50  = quantile(eff_median, probs = 0.25),
    uci50  = quantile(eff_median, probs = 0.75),
    lci95  = quantile(eff_median, probs = 0.025),
    uci95  = quantile(eff_median, probs = 0.975)
  ) %>%
  mutate(Param_neigh = paste(Param.type, Neigh12, sep = "_")) %>% 
  # convert to exponential scale        
  mutate_at(vars(median, contains("ci")), exp)

cumm.eff.alpha.beta_intra <-
  cumm.eff.i.f.n %>%
  filter(Neigh12 %in% c("ii", "iii", "iij", "iji")) %>% 
  # halve beta_iij and beta_iji before summation
  mutate_at(vars(starts_with("eff_")), 
            ~ifelse(Neigh12 %in% c("iij", "iji"), ./2, .)) %>% 
  # summation
  group_by(Plot, Year, Focal, Tag.tree, Param.type) %>%
  summarise_at(vars(starts_with("eff_")), ~exp(sum(.))) %>% 
  group_by(Focal, Param.type) %>% 
  summarise(
    median = median(eff_median),
    lci50  = quantile(eff_median, probs = 0.25),
    uci50  = quantile(eff_median, probs = 0.75),
    lci95  = quantile(eff_median, probs = 0.025),
    uci95  = quantile(eff_median, probs = 0.975)
  )
cumm.eff.alpha.beta_inter <-
  cumm.eff.i.f.n %>%
  filter(Neigh12 %in% c("ij", "iij", "iji", "ijj", "ijk")) %>% 
  # halve beta_iij and beta_iji before summation
  mutate_at(vars(starts_with("eff_")), 
            ~ifelse(Neigh12 %in% c("iij", "iji"), ./2, .)) %>% 
  # summation
  group_by(Plot, Year, Focal, Tag.tree, Param.type) %>%
  summarise_at(vars(starts_with("eff_")), ~exp(sum(.))) %>% 
  group_by(Focal, Param.type) %>% 
  summarise(
    median = median(eff_median),
    lci50  = quantile(eff_median, probs = 0.25),
    uci50  = quantile(eff_median, probs = 0.75),
    lci95  = quantile(eff_median, probs = 0.025),
    uci95  = quantile(eff_median, probs = 0.975)
  )

cumm.eff.param <-
  cumm.eff.i.f.n %>%
  group_by(Plot, Year, Focal, Tag.tree, Param.type) %>%
  summarise_at(vars(starts_with("eff_")), ~exp(sum(.))) %>% 
  group_by(Focal, Param.type) %>% 
  summarise(
    median = median(eff_median),
    lci50  = quantile(eff_median, probs = 0.25),
    uci50  = quantile(eff_median, probs = 0.75),
    lci95  = quantile(eff_median, probs = 0.025),
    uci95  = quantile(eff_median, probs = 0.975)
  ) %>% 
  bind_rows(
    cumm.eff.i.f.n %>%
      group_by(Plot, Year, Focal, Tag.tree) %>%
      summarise_at(vars(starts_with("eff_")), ~exp(sum(.))) %>% 
      group_by(Focal) %>% 
      summarise(
        median = median(eff_median),
        lci50  = quantile(eff_median, probs = 0.25),
        uci50  = quantile(eff_median, probs = 0.75),
        lci95  = quantile(eff_median, probs = 0.025),
        uci95  = quantile(eff_median, probs = 0.975)
      ) %>% 
      mutate(Param.type = "Total")
  )

# focal species with HOI model as best supported
fig4.focal.labs <- c("ARCHCL", "PRUNPO", "GIRONE", "GARCPA")

# contours 
a.sim <- exp(seq(log(0.1), log(100), length.out = 100))
ab.sim <- data.frame(a = a.sim, 
                     b =  1 / a.sim, 
                     b0.5 = 0.5 / a.sim,
                     b2 = 2 / a.sim)
cumm.eff.alpha.beta_intra.wide <- 
  cumm.eff.alpha.beta_intra %>% 
  pivot_wider(names_from = Param.type,
              values_from = median:uci95) %>% 
  mutate(Lab = Focal %in% fig4.focal.labs)
cumm.eff.alpha.beta_inter.wide <- 
  cumm.eff.alpha.beta_inter %>% 
  pivot_wider(names_from = Param.type,
              values_from = median:uci95) %>% 
  mutate(Lab = Focal %in% fig4.focal.labs)
cumm.eff.param.wide <- 
  cumm.eff.param %>% 
  pivot_wider(names_from = Param.type,
              values_from = median:uci95) %>% 
  mutate(Lab = Focal %in% fig4.focal.labs)

# truncate 95th percentiles to plot limit
# to draw as arrow heads
cumm.eff.plot.lower.lim <- 0.5
cumm.eff.plot.upper.lim <- 2
cumm.eff.alpha.beta_intra.wide2 <- 
  cumm.eff.alpha.beta_intra.wide %>% 
  select(contains("median"), contains("ci95")) %>% 
  mutate_at(vars(contains("lci95")), 
            ~ifelse(. < cumm.eff.plot.lower.lim,
                    cumm.eff.plot.lower.lim,
                    -1000)) %>% 
  mutate_at(vars(contains("uci95")), 
            ~ifelse(. > cumm.eff.plot.upper.lim,
                    cumm.eff.plot.upper.lim,
                    1000))
cumm.eff.alpha.beta_inter.wide2 <- 
  cumm.eff.alpha.beta_inter.wide %>% 
  select(contains("median"), contains("ci95")) %>% 
  mutate_at(vars(contains("lci95")), 
            ~ifelse(. < cumm.eff.plot.lower.lim,
                    cumm.eff.plot.lower.lim,
                    -1000)) %>% 
  mutate_at(vars(contains("uci95")), 
            ~ifelse(. > cumm.eff.plot.upper.lim,
                    cumm.eff.plot.upper.lim,
                    1000))
cumm.eff.param.wide2 <- 
  cumm.eff.param.wide %>% 
  select(contains("median"), contains("ci95")) %>% 
  mutate_at(vars(contains("lci95")), 
            ~ifelse(. < cumm.eff.plot.lower.lim,
                    cumm.eff.plot.lower.lim,
                    -1000)) %>% 
  mutate_at(vars(contains("uci95")), 
            ~ifelse(. > cumm.eff.plot.upper.lim,
                    cumm.eff.plot.upper.lim,
                    1000))

cumm.beta_intra.alpha.plot <- 
  ggplot(cumm.eff.alpha.beta_intra.wide, aes(median_alpha, median_beta)) +
  geom_line(data = ab.sim, aes(a, b), linetype = 2) +
  geom_line(data = ab.sim, aes(a, b0.5), linetype = 3) +
  geom_line(data = ab.sim, aes(a, b2), linetype = 3) +
  geom_errorbarh(aes(xmin = lci50_alpha, xmax = uci50_alpha), colour = "darkgrey", height = 0, size = 2) +
  geom_errorbarh(aes(xmin = lci95_alpha, xmax = uci95_alpha), colour = "darkgrey", height = 0) +
  geom_errorbar(aes(ymin = lci50_beta, ymax = uci50_beta), colour = "darkgrey", width = 0, size = 2) +
  geom_errorbar(aes(ymin = lci95_beta, ymax = uci95_beta), colour = "darkgrey", width = 0) +
  geom_point(aes(fill = Lab), size = 2, pch = 21, colour = "black") +
  # add arrows if extending beyond plot limits
  geom_segment(data = cumm.eff.alpha.beta_intra.wide2,
               aes(x = lci95_alpha+0.01, y = median_beta,
                   xend = lci95_alpha, yend = median_beta), 
               colour = "darkgrey",
               arrow = arrow(length = unit(0.05, "inches"))) +
  geom_segment(data = cumm.eff.alpha.beta_intra.wide2,
               aes(x = median_alpha, y = lci95_beta+0.01,
                   xend = median_alpha, yend = lci95_beta), 
               colour = "darkgrey",
               arrow = arrow(length = unit(0.05, "inches"))) +
  geom_segment(data = cumm.eff.alpha.beta_intra.wide2,
               aes(x = uci95_alpha-0.01, y = median_beta,
                   xend = uci95_alpha, yend = median_beta), 
               colour = "darkgrey",
               arrow = arrow(length = unit(0.05, "inches"))) +
  geom_segment(data = cumm.eff.alpha.beta_intra.wide2,
               aes(x = median_alpha, y = uci95_beta-0.01,
                   xend = median_alpha, yend = uci95_beta), 
               colour = "darkgrey",
               arrow = arrow(length = unit(0.05, "inches"))) +
  annotate("text", x = 1, y = 0.515, label = expression(0.5 %*% ""), 
           angle = -45, hjust = 1, vjust = 0, size = 3) +
  annotate("text", x = 1.95, y = 0.52, label = expression(1 %*% ""), 
           angle = -45, hjust = 1, vjust = 0, size = 3) +
  annotate("text", x = 1.95, y = 1.05, label = expression(2 %*% ""), 
           angle = -45, hjust = 1, vjust = 0, size = 3) +
  scale_y_log10(breaks = seq(0.5, 2, 0.1),
                labels = c(0.5, rep("", 4), 1, rep("", 9), 2)) +
  scale_x_log10(breaks = seq(0.5, 2, 0.1),
                labels = c(0.5, rep("", 4), 1, rep("", 9), 2)) +
  scale_fill_manual(values = c("white", "blue")) +
  coord_equal(ylim = c(cumm.eff.plot.lower.lim, cumm.eff.plot.upper.lim), 
              xlim = c(cumm.eff.plot.lower.lim, cumm.eff.plot.upper.lim), 
              expand = FALSE) +
  theme_classic() +
  theme(plot.margin = margin(14, 10, 2, 2),
        axis.title = element_blank(),
        legend.position = "none")

cumm.beta_inter.alpha.plot <- 
  ggplot(cumm.eff.alpha.beta_inter.wide, aes(median_alpha, median_beta)) +
  geom_line(data = ab.sim, aes(a, b), linetype = 2) +
  geom_line(data = ab.sim, aes(a, b0.5), linetype = 3) +
  geom_line(data = ab.sim, aes(a, b2), linetype = 3) +
  geom_errorbarh(aes(xmin = lci50_alpha, xmax = uci50_alpha), colour = "darkgrey", height = 0, size = 2) +
  geom_errorbarh(aes(xmin = lci95_alpha, xmax = uci95_alpha), colour = "darkgrey", height = 0) +
  geom_errorbar(aes(ymin = lci50_beta, ymax = uci50_beta), colour = "darkgrey", width = 0, size = 2) +
  geom_errorbar(aes(ymin = lci95_beta, ymax = uci95_beta), colour = "darkgrey", width = 0) +
  geom_point(aes(fill = Lab), size = 2, pch = 21, colour = "black") +
  # add arrows if extending beyond plot limits
  geom_segment(data = cumm.eff.alpha.beta_inter.wide2,
               aes(x = lci95_alpha+0.01, y = median_beta,
                   xend = lci95_alpha, yend = median_beta), 
               colour = "darkgrey",
               arrow = arrow(length = unit(0.05, "inches"))) +
  geom_segment(data = cumm.eff.alpha.beta_inter.wide2,
               aes(x = median_alpha, y = lci95_beta+0.01,
                   xend = median_alpha, yend = lci95_beta), 
               colour = "darkgrey",
               arrow = arrow(length = unit(0.05, "inches"))) +
  geom_segment(data = cumm.eff.alpha.beta_inter.wide2,
               aes(x = uci95_alpha-0.01, y = median_beta,
                   xend = uci95_alpha, yend = median_beta), 
               colour = "darkgrey",
               arrow = arrow(length = unit(0.05, "inches"))) +
  geom_segment(data = cumm.eff.alpha.beta_inter.wide2,
               aes(x = median_alpha, y = uci95_beta-0.01,
                   xend = median_alpha, yend = uci95_beta), 
               colour = "darkgrey",
               arrow = arrow(length = unit(0.05, "inches"))) +
  annotate("text", x = 1, y = 0.515, label = expression(0.5 %*% ""), 
           angle = -45, hjust = 1, vjust = 0, size = 3) +
  annotate("text", x = 1.95, y = 0.52, label = expression(1 %*% ""), 
           angle = -45, hjust = 1, vjust = 0, size = 3) +
  annotate("text", x = 1.95, y = 1.05, label = expression(2 %*% ""), 
           angle = -45, hjust = 1, vjust = 0, size = 3) +
  scale_y_log10(breaks = seq(0.5, 2, 0.1),
                labels = c(0.5, rep("", 4), 1, rep("", 9), 2)) +
  scale_x_log10(breaks = seq(0.5, 2, 0.1),
                labels = c(0.5, rep("", 4), 1, rep("", 9), 2)) +
  scale_fill_manual(values = c("white", "blue")) +
  coord_equal(ylim = c(cumm.eff.plot.lower.lim, cumm.eff.plot.upper.lim), 
              xlim = c(cumm.eff.plot.lower.lim, cumm.eff.plot.upper.lim), 
              expand = FALSE) +
  theme_classic() +
  theme(plot.margin = margin(14, 10, 2, 2),
        axis.title = element_blank(),
        legend.position = "none")

cumm.beta.alpha.plot <- 
  ggplot(cumm.eff.param.wide, aes(median_alpha, median_beta)) +
  geom_line(data = ab.sim, aes(a, b), linetype = 2) +
  geom_line(data = ab.sim, aes(a, b0.5), linetype = 3) +
  geom_line(data = ab.sim, aes(a, b2), linetype = 3) +
  geom_errorbarh(aes(xmin = lci50_alpha, xmax = uci50_alpha), colour = "darkgrey", height = 0, size = 2) +
  geom_errorbarh(aes(xmin = lci95_alpha, xmax = uci95_alpha), colour = "darkgrey", height = 0) +
  geom_errorbar(aes(ymin = lci50_beta, ymax = uci50_beta), colour = "darkgrey", width = 0, size = 2) +
  geom_errorbar(aes(ymin = lci95_beta, ymax = uci95_beta), colour = "darkgrey", width = 0) +
  geom_point(aes(fill = Lab), size = 2, pch = 21, colour = "black") +
  # add arrows if extending beyond plot limits
  geom_segment(data = cumm.eff.param.wide2,
               aes(x = lci95_alpha+0.01, y = median_beta,
                   xend = lci95_alpha, yend = median_beta), 
               colour = "darkgrey",
               arrow = arrow(length = unit(0.05, "inches"))) +
  geom_segment(data = cumm.eff.param.wide2,
               aes(x = median_alpha, y = lci95_beta+0.01,
                   xend = median_alpha, yend = lci95_beta), 
               colour = "darkgrey",
               arrow = arrow(length = unit(0.05, "inches"))) +
  geom_segment(data = cumm.eff.param.wide2,
               aes(x = uci95_alpha-0.01, y = median_beta,
                   xend = uci95_alpha, yend = median_beta), 
               colour = "darkgrey",
               arrow = arrow(length = unit(0.05, "inches"))) +
  geom_segment(data = cumm.eff.param.wide2,
               aes(x = median_alpha, y = uci95_beta-0.01,
                   xend = median_alpha, yend = uci95_beta), 
               colour = "darkgrey",
               arrow = arrow(length = unit(0.05, "inches"))) +
  annotate("text", x = 1, y = 0.515, label = expression(0.5 %*% ""), 
           angle = -45, hjust = 1, vjust = 0, size = 3) +
  annotate("text", x = 1.95, y = 0.52, label = expression(1 %*% ""), 
           angle = -45, hjust = 1, vjust = 0, size = 3) +
  annotate("text", x = 1.95, y = 1.05, label = expression(2 %*% ""), 
           angle = -45, hjust = 1, vjust = 0, size = 3) +
  # geom_text_repel(aes(median_alpha, median_beta, label = Focal)) +
  scale_y_log10(breaks = seq(0.5, 2, 0.1),
                labels = c(0.5, rep("", 4), 1, rep("", 9), 2)) +
  scale_x_log10(breaks = seq(0.5, 2, 0.1),
                labels = c(0.5, rep("", 4), 1, rep("", 9), 2)) +
  scale_fill_manual(values = c("white", "blue")) +
  coord_equal(ylim = c(cumm.eff.plot.lower.lim, cumm.eff.plot.upper.lim), 
              xlim = c(cumm.eff.plot.lower.lim, cumm.eff.plot.upper.lim), 
              expand = FALSE) +
  theme_classic() +
  theme(plot.margin = margin(14, 10, 2, 2),
        axis.title = element_blank(),
        legend.position = "none")

fig4 <- 
  ggarrange(cumm.beta_intra.alpha.plot, 
            cumm.beta_inter.alpha.plot, 
            cumm.beta.alpha.plot,
            nrow = 1,
            labels = 
              c("(a) Intraspecific",
                "(b) Interspecific",
                "(c) All"),
            font.label = list(size = 12),
            hjust = 0, vjust = 1)
fig4 <- annotate_figure(
  fig4,
  left = text_grob("Cumulative instantaneous\neffects of HOIs", size = 11, rot = 90),
  bottom = text_grob(
    "Cumulative instantaneous effects of direct interactions",
    size = 11,
    vjust = 0
  )
)

fig4