# Reproduce Fig. 2

source("code/get_tidy_median_qi.R")

X.scale <- 
  data.frame(Scale = in.dat$scale) %>% 
  rownames_to_column("Predictor")
dat.org.scale <- generate.model.matrix(mandai, scale.X = FALSE)$in.dat
beta.coef <- get.tidy.median.qi(beta.mlm, "HOI")
beta.coef.backscale <- 
  beta.coef %>% 
  select(Focal, Param.type, Predictor=Param, Neigh12, Estimate, Q2.5, Q97.5) %>% 
  left_join(X.scale) %>% 
  mutate_at(vars(c("Estimate", "Q2.5", "Q97.5")), ~./Scale)

coef_posterior <- 
  beta.mlm %>% 
  gather_draws(r_Focal__a[Focal, Predictor]) %>% 
  # change signs of coefs
  mutate(.value = - .value) %>% 
  ungroup %>% 
  filter(Predictor != "Intercept") %>% 
  mutate(Param.type = ifelse(Predictor %in% neigh, "alpha", "beta"),
         Neigh1 = ifelse(word(Predictor, 1, sep = "_")==Focal, "i", "j"),
         Neigh2 = ifelse(word(Predictor, 2, sep = "_")==Focal, "i", 
                         ifelse(word(Predictor, 2, sep = "_")=="2", Neigh1,
                                ifelse(Neigh1=="i", "j", "k"))),
         Neigh12 = paste0("i", Neigh1, ifelse(is.na(Neigh2), "", Neigh2)))

alpha_ii_coefs <- 
  coef_posterior %>% 
  filter(Param.type == "alpha",
         Neigh12 == "ii") %>% 
  mutate(Param_neigh = paste(Param.type, Neigh12, sep = "_"))  %>% 
  select(Focal, Predictor, Param_neigh, .chain, .iteration, .draw, Est_alpha = .value)

alpha_ij_coefs <- 
  coef_posterior %>% 
  filter(Param.type == "alpha",
         Neigh12 == "ij") %>% 
  mutate(Param_neigh = paste(Param.type, Neigh12, sep = "_"))  %>% 
  select(Focal, Predictor, Param_neigh, .chain, .iteration, .draw, Est_alpha = .value)

beta_all_coefs <- 
  coef_posterior %>% 
  filter(Param.type == "beta") %>% 
  # halve beta_iij and beta_iji before summation
  mutate_at(vars(.value),
            ~ifelse(Neigh12 %in% c("iij", "iji"), ./2, .)) %>%
  mutate(Param_neigh = paste(Param.type, Neigh12, sep = "_"))  %>% 
  select(Focal, Predictor, Param_neigh, .chain, .iteration, .draw, Est_beta = .value) %>% 
  # sloppy trick to match betas to alphas by defining direct neighbour
  separate(Predictor, c("Predictor", "Predictor_2"), sep = "_")

beta_all_coefs_flip <- 
  beta_all_coefs %>% 
  rename(Predictor_2 = Predictor,
         Predictor = Predictor_2)
beta_all_coefs_comb <- bind_rows(beta_all_coefs, beta_all_coefs_flip)

intra_coefs_comp <- 
  alpha_ii_coefs %>% 
  left_join(beta_all_coefs_comb %>% 
              filter(Param_neigh %in% c("beta_iii", "beta_iij", "beta_iji")),
            by = c("Focal" = "Focal", 
                   "Predictor" = "Predictor",
                   ".chain" = ".chain",
                   ".iteration" = ".iteration",
                   ".draw" = ".draw")) %>% 
  group_by(Focal, Predictor, Predictor_2) %>% 
  summarise(Est_alpha_median = median(Est_alpha),
            Est_alpha_lci95  = quantile(Est_alpha, probs = 0.025),
            Est_alpha_uci95  = quantile(Est_alpha, probs = 0.975),
            Est_alpha_lci50  = quantile(Est_alpha, probs = 0.25),
            Est_alpha_uci50  = quantile(Est_alpha, probs = 0.75),
            Est_beta_median  = median(Est_beta),
            Est_beta_lci95   = quantile(Est_beta, probs = 0.025),
            Est_beta_uci95   = quantile(Est_beta, probs = 0.975),
            Est_beta_lci50   = quantile(Est_beta, probs = 0.25),
            Est_beta_uci50   = quantile(Est_beta, probs = 0.75)) %>% 
  # key to distinguish effects that involve OTHERS
  mutate(has_OTHERS = as.numeric(Predictor == "OTHERS" | Predictor_2 == "OTHERS")) %>% 
  # remove OTHERS for main text
  filter(Predictor   != "OTHERS",
         Predictor_2 != "OTHERS")

inter_coefs_comp <- 
  alpha_ij_coefs %>% 
  left_join(beta_all_coefs_comb %>% 
              filter(Param_neigh %in% c("beta_ijj", "beta_iij", "beta_iji", "beta_ijk")),
            by = c("Focal" = "Focal", 
                   "Predictor" = "Predictor",
                   ".chain" = ".chain",
                   ".iteration" = ".iteration",
                   ".draw" = ".draw")) %>% 
  group_by(Focal, Predictor, Predictor_2) %>% 
  summarise(Est_alpha_median = median(Est_alpha),
            Est_alpha_lci95  = quantile(Est_alpha, probs = 0.025),
            Est_alpha_uci95  = quantile(Est_alpha, probs = 0.975),
            Est_alpha_lci50  = quantile(Est_alpha, probs = 0.25),
            Est_alpha_uci50  = quantile(Est_alpha, probs = 0.75),
            Est_beta_median  = median(Est_beta),
            Est_beta_lci95   = quantile(Est_beta, probs = 0.025),
            Est_beta_uci95   = quantile(Est_beta, probs = 0.975),
            Est_beta_lci50   = quantile(Est_beta, probs = 0.25),
            Est_beta_uci50   = quantile(Est_beta, probs = 0.75)) %>% 
  # key to distinguish effects that involve OTHERS
  mutate(has_OTHERS = as.numeric(Predictor == "OTHERS" | Predictor_2 == "OTHERS")) %>% 
  # remove OTHERS for main text
  filter(Predictor   != "OTHERS",
         Predictor_2 != "OTHERS")

# proportions(table(intra_coefs_comp$Est_alpha_median >= 0))
# proportions(table(intra_coefs_comp$Est_beta_median >= 0))
# proportions(table(inter_coefs_comp$Est_alpha_median >= 0))
# proportions(table(inter_coefs_comp$Est_beta_median >= 0))

# Plots
xlims <- 
  range(
    c(range(intra_coefs_comp$Est_alpha_median),
      range(inter_coefs_comp$Est_alpha_median))
  )
ylims <- 
  range(
    c(range(intra_coefs_comp$Est_beta_median),
      range(inter_coefs_comp$Est_beta_median))
  )
xylims <- range(xlims, ylims)
xylims <- xylims + 0.05 * sign(xylims)  # expand lims

# truncate 95th percentiles to plot limit
# to draw as arrow heads
intra_coefs_comp2 <- 
  intra_coefs_comp %>% 
  mutate_at(vars(contains("lci95")), 
            ~ifelse(. < xylims[1],
                    xylims[1],
                    -1000)) %>% 
  mutate_at(vars(contains("uci95")), 
            ~ifelse(. > xylims[2],
                    xylims[2],
                    1000))
inter_coefs_comp2 <- 
  inter_coefs_comp %>% 
  mutate_at(vars(contains("lci95")), 
            ~ifelse(. < xylims[1],
                    xylims[1],
                    -1000)) %>% 
  mutate_at(vars(contains("uci95")), 
            ~ifelse(. > xylims[2],
                    xylims[2],
                    1000))

p1 <- 
  ggplot() +
  geom_errorbar(data = intra_coefs_comp,
                aes(Est_alpha_median, 
                    ymin = Est_beta_lci95, ymax = Est_beta_uci95),
                colour = "darkgrey", width = 0) +
  geom_errorbar(data = intra_coefs_comp,
                aes(Est_alpha_median, 
                    ymin = Est_beta_lci50, ymax = Est_beta_uci50),
                colour = "darkgrey", width = 0, size = 2) +
  geom_errorbarh(data = intra_coefs_comp,
                 aes(y = Est_beta_median, 
                     xmin = Est_alpha_lci95, xmax = Est_alpha_uci95),
                 colour = "darkgrey", height = 0) +
  geom_errorbarh(data = intra_coefs_comp,
                 aes(y = Est_beta_median, 
                     xmin = Est_alpha_lci50, xmax = Est_alpha_uci50),
                 colour = "darkgrey", height = 0, size = 2) + 
  # add arrows if extending beyond plot limits
  geom_segment(data = intra_coefs_comp2,
               aes(x = Est_alpha_lci95+0.01, y = Est_beta_median,
                   xend = Est_alpha_lci95, yend = Est_beta_median), 
               colour = "darkgrey",
               arrow = arrow(length = unit(0.05, "inches"))) +
  geom_segment(data = intra_coefs_comp2,
               aes(x = Est_alpha_median, y = Est_beta_lci95+0.01,
                   xend = Est_alpha_median, yend = Est_beta_lci95), 
               colour = "darkgrey",
               arrow = arrow(length = unit(0.05, "inches"))) +
  geom_segment(data = intra_coefs_comp2,
               aes(x = Est_alpha_uci95-0.01, y = Est_beta_median,
                   xend = Est_alpha_uci95, yend = Est_beta_median), 
               colour = "darkgrey",
               arrow = arrow(length = unit(0.05, "inches"))) +
  geom_segment(data = intra_coefs_comp2,
               aes(x = Est_alpha_median, y = Est_beta_uci95-0.01,
                   xend = Est_alpha_median, yend = Est_beta_uci95), 
               colour = "darkgrey",
               arrow = arrow(length = unit(0.05, "inches"))) +
  geom_vline(xintercept = 0, linetype = 2) +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_point(data = intra_coefs_comp,
             aes(Est_alpha_median, Est_beta_median, shape = as.character(has_OTHERS)), 
             size = 2, colour = "white", fill = "black") +
  scale_shape_manual(values = c(21, 4)) +
  labs(x = expression(alpha[ii]), 
       # y = expression(alpha[ii] + beta[ii.]~bar(A[.]))) +
       y = expression(beta[ii.])) +
  coord_equal(expand = FALSE, xlim = xylims, ylim = xylims) +
  theme_classic() +
  theme(plot.margin = margin(14, 10, 2, 2),
        axis.title = element_text(size = 12),
        legend.position = "none")

p2 <- 
  ggplot(inter_coefs_comp) +
  geom_errorbar(aes(Est_alpha_median, 
                    ymin = Est_beta_lci95, ymax = Est_beta_uci95),
                colour = "darkgrey", width = 0) +
  geom_errorbar(aes(Est_alpha_median,
                    ymin = Est_beta_lci50, ymax = Est_beta_uci50),
                colour = "darkgrey", width = 0, size = 2) +
  geom_errorbarh(aes(y = Est_beta_median, 
                     xmin = Est_alpha_lci95, xmax = Est_alpha_uci95),
                 colour = "darkgrey", height = 0) +
  geom_errorbarh(aes(y = Est_beta_median, 
                     xmin = Est_alpha_lci50, xmax = Est_alpha_uci50),
                 colour = "darkgrey", height = 0, size = 2) +
  # add arrows if extending beyond plot limits
  geom_segment(data = inter_coefs_comp2,
               aes(x = Est_alpha_lci95+0.01, y = Est_beta_median,
                   xend = Est_alpha_lci95, yend = Est_beta_median), 
               colour = "darkgrey",
               arrow = arrow(length = unit(0.05, "inches"))) +
  geom_segment(data = inter_coefs_comp2,
               aes(x = Est_alpha_median, y = Est_beta_lci95+0.01,
                   xend = Est_alpha_median, yend = Est_beta_lci95), 
               colour = "darkgrey",
               arrow = arrow(length = unit(0.05, "inches"))) +
  geom_segment(data = inter_coefs_comp2,
               aes(x = Est_alpha_uci95-0.01, y = Est_beta_median,
                   xend = Est_alpha_uci95, yend = Est_beta_median), 
               colour = "darkgrey",
               arrow = arrow(length = unit(0.05, "inches"))) +
  geom_segment(data = inter_coefs_comp2,
               aes(x = Est_alpha_median, y = Est_beta_uci95-0.01,
                   xend = Est_alpha_median, yend = Est_beta_uci95), 
               colour = "darkgrey",
               arrow = arrow(length = unit(0.05, "inches"))) +
  geom_vline(xintercept = 0, linetype = 2) +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_point(aes(Est_alpha_median, Est_beta_median, shape = as.character(has_OTHERS)), 
             size = 2, colour = "white", fill = "black") +
  scale_shape_manual(values = c(21, 4)) +
  labs(x = expression(alpha[ij]), 
       y = expression(beta[ij.])) +
  coord_equal(expand = FALSE, xlim = xylims, ylim = xylims) +
  theme_classic() +
  theme(plot.margin = margin(14, 10, 2, 2),
        axis.title = element_text(size = 12),
        legend.position = "none")

alpha_beta_coef_comp <- 
  ggarrange(
    p1, p2,
    nrow = 1,
    labels =
      c("(a) Intraspecific",
        "(b) Interspecific"),
    font.label = list(size = 12),
    hjust = -0.1,
    vjust = 1
  )

alpha_beta_coef_comp
