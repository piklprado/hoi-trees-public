# Reproduce Fig. 3

source("code/predict_dbh_example.R")

# population-level growth curve
newdat.pop <- 
  data.frame(DBH_mi = seq(0.01, max(mandai$DBH_mi), length.out = 100) / sd(mandai$DBH_mi)) %>% 
  mutate(Estimate = apply(
    inv_modulus_trans(
      predict(beta.mlm, newdata = ., re_formula = NA, summary = FALSE), 
      0.55
    ), 2, median)) %>% 
  # backtransfomations
  mutate(DBH_mi = DBH_mi * sd(mandai$DBH_mi)) 

# focal-level prediction
cumm.eff.growth.curve <- growth.curve.cumm.pred(beta.mlm, mandai, DBH.range = "species")

# calculate peak growth reduction 
peak.growth <- 
  cumm.eff.growth.curve %>% 
  group_by(Focal, Pred.type) %>% 
  filter(Estimate == max(Estimate)) %>% 
  select(Focal, Pred.type, Estimate) %>% 
  pivot_wider(names_from = Pred.type,
              values_from = Estimate) %>% 
  mutate(Reduction = `D + HOI_intra + HOI_inter` - Max) %>% 
  arrange(Reduction)
smallD.growth <- 
  cumm.eff.growth.curve %>% 
  group_by(Focal, Pred.type) %>% 
  filter(DBH_mi == min(DBH_mi),
         Pred.type == "Max") %>% 
  select(Focal, Estimate) %>% 
  arrange(Estimate)

cumm.eff.growth.curve.plot <- 
  ggplot(cumm.eff.growth.curve %>% 
           filter(Pred.type != "D + HOI_intra")) +
  facet_wrap(~ Focal, nrow = 2, scale = "free") +
  geom_line(data = newdat.pop, aes(DBH_mi, Estimate),
            colour = "grey", size = 0.5) +
  geom_line(aes(DBH_mi, Estimate, colour = Pred.type), size = 0.3) +
  labs(x = expression(paste("Diameter at breast height (cm)")), 
       y = expression(paste("Absolute instantaneous growth rate (cm ", yr^"-1", ")"))) +
  scale_colour_manual(values = c("black", "red", "blue"),
                      name = "Biotic background",
                      aesthetics = c("colour", "fill")) +
  scale_y_continuous(expand = expansion(0,0), limits = c(0, NA)) +
  scale_x_continuous(expand = expansion(0,0), limits = c(0, NA)) +
  theme_classic() +
  theme(legend.position = "none",
        axis.text = element_text(size = 8),
        strip.background = element_rect(fill = NA, colour = NA))

cumm.eff.growth.curve.plot