# script to summarise Bayes R2

null.R2 <- bayes_R2(null.mlm)
null.R2.focal <- 
  lapply(focals, function(x) {
    newdat <- 
      in.dat$in.dat %>% 
      filter(Focal == x)
    bayes_R2(null.mlm, newdata = newdat)
  })
alpha.R2 <- bayes_R2(alpha.mlm)
alpha.R2.focal <- 
  lapply(focals, function(x) {
    newdat <- 
      in.dat$in.dat %>% 
      filter(Focal == x)
    bayes_R2(alpha.mlm, newdata = newdat)
  })
beta.R2 <- bayes_R2(beta.mlm)
beta.R2.focal <- 
  lapply(focals, function(x) {
    newdat <- 
      in.dat$in.dat %>% 
      filter(Focal == x)
    bayes_R2(beta.mlm, newdata = newdat)
  })

r2.comp <- 
  rbind(null.R2, do.call(rbind, null.R2.focal)) %>% 
  as.data.frame %>% 
  bind_cols(Species = c("All", focals), 
            Model = rep("Null", length(focals) + 1), 
            .) %>% 
  bind_rows(
    rbind(alpha.R2, do.call(rbind, alpha.R2.focal)) %>% 
      as.data.frame %>% 
      bind_cols(Species = c("All", focals), 
                Model = rep("Direct", length(focals) + 1), 
                .)
  ) %>% 
  bind_rows(
    rbind(beta.R2, do.call(rbind, beta.R2.focal)) %>% 
      as.data.frame %>% 
      bind_cols(Species = c("All", focals), 
                Model = rep("HOI", length(focals) + 1), 
                .) 
  ) %>% 
  mutate(Species = 
           factor(Species, 
                  levels = c("All", 
                             rev(focals[order(do.call(rbind, beta.R2.focal)[,1])]))),
         Model = factor(Model, 
                        levels = c("Null", "Direct", "HOI"),
                        labels = c("Null",
                                   "Direct-interaction-only",
                                   "HOI-inclusive")))