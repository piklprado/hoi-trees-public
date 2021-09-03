# Reproduce Fig. 1

r2_plot <- 
  ggplot(r2.comp, aes(Species, Estimate)) +
  geom_errorbar(aes(ymin = Q2.5, ymax = Q97.5, colour = Model), 
                position = position_dodge(width = 0.5), width = 0) +
  geom_point(aes(colour = Model), 
             position = position_dodge(width = 0.5), width = 0, pch = 21, fill = "white", size = 2) +
  scale_colour_manual(values = c("grey", "red", "blue")) +
  ylim(0, 1) +
  labs(x = "", y = expression(paste("Bayes ", italic(R)^2))) +
  coord_flip() +
  theme_classic() +
  theme(axis.title.y = element_blank(),
        axis.text = element_text(size = 6))

waic_plot <-
  ggplot(dWAIC.boot) +
  geom_hline(yintercept = c(-2, 2), linetype = 2, colour = "grey") +
  geom_hline(yintercept = 0, colour = "grey") +
  geom_errorbar(aes(ymin = lci95, ymax = uci95, x = Focal, colour = Model), 
                width = 0, position = position_dodge(width = 0.5)) +
  geom_errorbar(aes(ymin = lci50, ymax = uci50, x = Focal, colour = Model), 
                width = 0, size = 2, position = position_dodge(width = 0.5)) +
  geom_point(aes(Focal, median, colour = Model), 
             pch = 21, fill = "white", size = 2, position = position_dodge(width = 0.5)) +
  coord_flip() +
  geom_text(data = dWAIC.actual, 
            aes(Focal, 1.4*min(dWAIC.boot$lci95), label = round(dWAIC.beta,1)),
            colour = "blue", hjust = "right", size = 2) +
  geom_text(data = dWAIC.actual, 
            aes(Focal, 1.01*min(dWAIC.boot$lci95), label = round(dWAIC.alpha,1)),
            colour = "red", hjust = "right", size = 2) +
  labs(x = "", y = expression(Delta~italic(WAIC))) +
  scale_colour_manual(values = c("red", "blue")) +
  ylim(min(dWAIC.boot$lci95)*1.7, NA) +
  theme_classic() +
  theme(axis.title.y = element_blank(),
        axis.text = element_text(size = 6),
        legend.position = "none")

loo_plot <-
  ggplot(dloo.boot) +
  geom_hline(yintercept = c(-2, 2), linetype = 2, colour = "grey") +
  geom_hline(yintercept = 0, colour = "grey") +
  geom_errorbar(aes(ymin = lci95, ymax = uci95, x = Focal, colour = Model), 
                width = 0, position = position_dodge(width = 0.5)) +
  geom_errorbar(aes(ymin = lci50, ymax = uci50, x = Focal, colour = Model), 
                width = 0, size = 2, position = position_dodge(width = 0.5)) +
  geom_point(aes(Focal, median, colour = Model), 
             pch = 21, fill = "white", size = 2, position = position_dodge(width = 0.5)) +
  coord_flip() +
  geom_text(data = dloo.actual, 
            aes(Focal, 1.4*min(dloo.boot$lci95), label = round(dloo.beta,1)),
            colour = "blue", hjust = "right", size = 2) +
  geom_text(data = dloo.actual, 
            aes(Focal, 1.01*min(dloo.boot$lci95), label = round(dloo.alpha,1)),
            colour = "red", hjust = "right", size = 2) +
  labs(x = "", y = expression(Delta~italic(LOOIC))) +
  scale_colour_manual(values = c("red", "blue")) +
  ylim(min(dloo.boot$lci95)*1.7, NA) +
  theme_classic() +
  theme(axis.title.y = element_blank(),
        axis.text = element_text(size = 6),
        legend.position = "none")

model_comp_plot <- 
  ggarrange(r2_plot, waic_plot, loo_plot,
            ncol = 3,
            labels = paste0("(", letters[1:3], ")"),
            font.label = list(size = 10),
            hjust = 0, vjust = 0.5, 
            common.legend = TRUE)

model_comp_plot
