# Reproduce Fig. 5

source("code/D_over_time.R")

# Diameter over time simulation -------------------------------------------

# using Lai et al. 2020 predictions
N_init_lai2020_q5 <- readRDS("data/N_init_lai2020_q5.rds")
N_init_lai2020_median <- readRDS("data/N_init_lai2020_median.rds")
N_init_lai2020_q95 <- readRDS("data/N_init_lai2020_q95.rds")

focals_vec_func <- function(N_init) {
  as.character(rep(spp_list$Abbreviation[match(names(N_init), spp_list$Species)], N_init))
}
focals_lai2020_q5 <- focals_vec_func(N_init_lai2020_q5$pred.focal.adj)
focals_lai2020_median <- focals_vec_func(N_init_lai2020_median$pred.focal.adj)
focals_lai2020_q95 <- focals_vec_func(N_init_lai2020_q95$pred.focal.adj)

D_sim_lai2020_median <-
  simulateD_focal_garden_sep_mod(
    model.n = null.mlm,
    model.a = alpha.mlm,
    model.b = beta.mlm,
    t.step = 365*2, t.max = 2,
    focals = focals_lai2020_median
  )
D_sim_lai2020_q5 <-
  simulateD_focal_garden_sep_mod(
    model.n = null.mlm,
    model.a = alpha.mlm,
    model.b = beta.mlm,
    t.step = 365*2, t.max = 2,
    focals = focals_lai2020_q5
  )
D_sim_lai2020_q95 <-
  simulateD_focal_garden_sep_mod(
    model.n = null.mlm,
    model.a = alpha.mlm,
    model.b = beta.mlm,
    t.step = 365*2, t.max = 2,
    focals = focals_lai2020_q95
  )


# Make plot ---------------------------------------------------------------

t_lim <- 2
D_lim <- 3.7

# average growth all else at their means
avg.growth <- 
  ranef(beta.mlm)$Focal[,,"a_Intercept"] %>%
  exp() %>% 
  as.data.frame %>% 
  rownames_to_column("Focal") %>% 
  select(Focal,
         Intercept = Estimate,
         Intercept_Q2.5 = Q2.5,
         Intercept_Q97.5 = Q97.5) %>% 
  right_join(cumm.eff.param)

# combine simulation across recruitment scenarios
D_sim_lai2020_raw <- 
  bind_rows(D_sim_lai2020_q5 %>% mutate(R = "Low"),
            D_sim_lai2020_median %>% mutate(R = "Median recruitment"),
            D_sim_lai2020_q95 %>% mutate(R = "High")) %>% 
  mutate(Focal = as.factor(Focal),
         R = fct_relevel(R, c("Low", "Median recruitment", "High"))) %>% 
  filter(Pred.type != "Null") %>% 
  mutate(Pred.type = 
           fct_recode(Pred.type, 
                      "Direct-interaction-only" = "Direct",
                      "HOI-inclusive" = "HOI")) %>% 
  group_by(Focal, Pred.type, t)

D_sim_lai2020 <- 
  D_sim_lai2020_raw %>% 
  group_by(Focal, R, Pred.type, t) %>%
  summarise(D = median(D))

# name labels at the end of simulation
D_sim_lai2020_label <-
  D_sim_lai2020_raw %>%
  ungroup() %>% 
  filter(t == t_lim) %>% 
  group_by(Focal, Pred.type, R, t) %>%
  summarise(D_median = median(D),
            D_max = max(D))
D_sim_lai2020_label_order <-
  avg.growth %>% distinct(Focal, Intercept) %>%  arrange(Intercept) %>% pull(Focal)

# density of size distributions
D_sim_lai2020_dens <- 
  D_sim_lai2020_raw %>%
  ungroup() %>% 
  filter(t == t_lim) %>% 
  mutate(t = 0)  # force density plots to be placed at t = 0

D_sim_lai2020 <- 
  D_sim_lai2020 %>% 
  ungroup %>% 
  mutate(Focal = fct_relevel(Focal, D_sim_lai2020_label_order))

D_sim_lai2020_direct_plot <- 
  ggplot(data = D_sim_lai2020 %>% 
           filter(Pred.type=="Direct-interaction-only"), 
         aes(t, D)) +
  facet_grid(R ~ .) +
  geom_half_violin(data = D_sim_lai2020_dens %>% 
                     filter(Pred.type=="Direct-interaction-only"),
                   aes(t, D),
                   draw_quantiles = c(0.25, 0.5, 0.75),
                   scale = "count", side = "r", 
                   fill = "lightgrey", colour = "white",
                   inherit.aes = FALSE) +
  geom_line(aes(colour = Focal)) +
  geom_text_repel(data = D_sim_lai2020_label %>% 
                    filter(Pred.type=="Direct-interaction-only"), 
                  aes(1.01 * t, D_median, label = Focal),
                  size = 2.5, direction = "y", nudge_x = 0.5, hjust = 0,
                  box.padding = 0.1,
                  colour = "darkgrey",
                  segment.color = "grey",
                  segment.size = 0.2,
                  inherit.aes = FALSE) +
  labs(y = "Diameter at breast height (cm)", x = "Time (yr)") +
  scale_colour_viridis_d(option = "cividis") +
  coord_cartesian(expand = FALSE) +
  scale_x_continuous(breaks = seq(0, t_lim, 1), limits = c(0, t_lim+1)) +
  scale_y_continuous(breaks = seq(0, D_lim, 1), limits = c(0.9, D_lim)) +
  theme_bw() +
  theme(legend.position = "none",
        panel.grid.minor = element_blank(),
        strip.text = element_blank(), 
        panel.spacing = unit(0, "cm"),
        plot.margin = margin(12, 0.5, 0.5, 10),
        panel.grid = element_line(size = 0.1),
        axis.title = element_text(size = 10))
D_sim_lai2020_hoi_plot <- 
  ggplot(data = D_sim_lai2020 %>% 
           filter(Pred.type=="HOI-inclusive"), 
         aes(t, D)) +
  facet_grid(R ~ .) +
  geom_half_violin(data = D_sim_lai2020_dens %>% 
                     filter(Pred.type=="HOI-inclusive"),
                   aes(t, D), 
                   draw_quantiles = c(0.25, 0.5, 0.75),
                   scale = "count", side = "r", 
                   fill = "lightgrey", colour = "white",
                   inherit.aes = FALSE) +
  geom_line(aes(colour = Focal)) +
  geom_text_repel(data = D_sim_lai2020_label %>% 
                    filter(Pred.type=="HOI-inclusive"), 
                  aes(1.01 * t, D_median, label = Focal),
                  size = 2.5, direction = "y", nudge_x = 0.5, hjust = 0,
                  box.padding = 0.1,
                  colour = "darkgrey",
                  segment.color = "grey",
                  segment.size = 0.2,
                  inherit.aes = FALSE) +
  labs(y = "Diameter at breast height (cm)", x = "Time (yr)") +
  scale_colour_viridis_d(option = "cividis") +
  coord_cartesian(expand = FALSE) +
  scale_x_continuous(breaks = seq(0, t_lim, 1), limits = c(0, t_lim+1)) +
  scale_y_continuous(breaks = seq(0, D_lim, 1), limits = c(0.9, D_lim)) +
  theme_bw() +
  theme(legend.position = "none",
        panel.grid.minor = element_blank(),
        strip.text = element_blank(), 
        panel.spacing = unit(0, "cm"),
        plot.margin = margin(12, 0.5, 0.5, 10),
        panel.grid = element_line(size = 0.1),
        axis.title = element_text(size = 10))

# calculate ratio between direct-only vs. HOI-inclusive
D_ratio_raw <- 
  D_sim_lai2020_raw %>% 
  pivot_wider(names_from = Pred.type,
              values_from = D) %>% 
  group_by(t, R) %>% 
  mutate(Ratio = `HOI-inclusive` / `Direct-interaction-only`,
         Focal = fct_relevel(Focal, D_sim_lai2020_label_order)) %>%
  select(-`Direct-interaction-only`, -`HOI-inclusive`)
D_ratio <- 
  D_ratio_raw %>%
  group_by(Focal, t, R) %>%
  summarise(Ratio = median(Ratio))
D_ratio_label <-
  D_ratio %>%
  filter(t == t_lim) %>% 
  group_by(Focal, R, t) %>%
  summarise(Ratio = median(Ratio))
D_ratio_plot <- 
  ggplot(data = D_ratio, aes(t, Ratio)) +
  facet_grid(R ~ .) +
  geom_line(aes(colour = Focal)) +
  geom_segment(aes(x = 0, y = 1, xend = 2, yend = 1), linetype = 2) +
  geom_text_repel(data = D_ratio_label,
                  aes(1.01 * t, Ratio, label = Focal),
                  size = 2.5, direction = "y", nudge_x = 0.5, hjust = 0,
                  box.padding = 0.1,
                  colour = "darkgrey",
                  segment.color = "grey",
                  segment.size = 0.2,
                  inherit.aes = FALSE) +
  labs(y = "Diameter ratio", x = "Time (yr)") +
  scale_colour_viridis_d(option = "cividis") +
  coord_cartesian(expand = FALSE) +
  scale_x_continuous(breaks = seq(0, t_lim, 1), limits = c(0, t_lim+1)) +
  scale_y_continuous(breaks = seq(0.5, 1.5, 0.25), limits = c(0.5, 1.6)) +
  theme_bw() +
  theme(legend.position = "none",
        panel.grid.minor = element_blank(),
        strip.background = element_blank(), 
        strip.text = element_text(size = 11),
        panel.spacing = unit(0, "cm"),
        plot.margin = margin(12, 0.5, 0.5, 10),
        panel.grid = element_line(size = 0.1),
        axis.title = element_text(size = 10))

D_sim_plot <- 
  ggarrange(D_sim_lai2020_direct_plot, 
            D_sim_lai2020_hoi_plot,
            D_ratio_plot,
            ncol = 3,
            labels = 
              c(
                "(a) Direct-interaction-only",
                "(b) HOI-inclusive",
                "(c) Ratio between (b) and (a)"
              ),
            font.label = list(size = 10),
            hjust = -0.1, vjust = 1.1,
            widths = c(1, 1, 1.1)
  )

D_sim_plot