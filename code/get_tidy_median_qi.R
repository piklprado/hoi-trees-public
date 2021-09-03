get.tidy.median.qi <- function(model, model.type) {
  out <- 
    reshape::melt(coefficients(model, robust = TRUE, probs = c(0.025, 0.975))$Focal, 
                  varnames = c("Focal", "Stat", "Param")) %>% 
    filter(str_detect(Param, "a_")) %>% 
    pivot_wider(names_from = Stat,
                values_from = value) %>% 
    mutate(Model = model.type,
           Param = str_remove(Param, "a_"),
           Param.type = ifelse(str_detect(Param, "Intercept"), Param, 
                               ifelse(Param %in% neigh, "alpha", "beta")),
           Neigh1 = ifelse(Param %in% c("Intercept", "sigma"), NA, 
                           ifelse(word(Param, 1, sep = "_")==Focal, "i", "j")),
           Neigh2 = ifelse(word(Param, 2, sep = "_")==Focal, "i", 
                           ifelse(word(Param, 2, sep = "_")=="2", Neigh1,
                                  ifelse(Neigh1=="i", "j", "k"))),
           Neigh12 = paste0("i", Neigh1, ifelse(is.na(Neigh2), "", Neigh2))) %>% 
    mutate(Param = factor(Param, levels = rev(c("Intercept", setdiff(unique(Param), "Intercept")))))
  return(out)
}
