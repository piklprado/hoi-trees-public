# script to summarise loo

null.loo <-  loo(null.mlm, pointwise = TRUE)
alpha.loo <- loo(alpha.mlm, pointwise = TRUE)
beta.loo <-  loo(beta.mlm, pointwise = TRUE)

# Pointwise loo
null.loo.pointwise  <- null.loo$pointwise[, "looic"]
alpha.loo.pointwise <- alpha.loo$pointwise[, "looic"]
beta.loo.pointwise  <- beta.loo$pointwise[, "looic"]
dloo.pointwise <- 
  beta.mlm$data %>% 
  select(Focal) %>% 
  mutate(null = null.loo.pointwise,
         alpha = alpha.loo.pointwise,
         beta  = beta.loo.pointwise) %>% 
  mutate(dloo.alpha = alpha - null,
         dloo.beta = beta - null)
dloo.actual.all <-
  dloo.pointwise %>% 
  mutate(Focal = "All") %>%
  group_by(Focal) %>% 
  summarise(dloo.alpha = sum(dloo.alpha),
            dloo.beta = sum(dloo.beta))
dloo.actual.focal <-
  dloo.pointwise %>% 
  group_by(Focal) %>% 
  summarise(dloo.alpha = sum(dloo.alpha),
            dloo.beta = sum(dloo.beta))
dloo.actual <- bind_rows(dloo.actual.all, dloo.actual.focal)

# Bootstrap pointwise loo
focal.N <- 
  dloo.pointwise %>% 
  group_by(Focal) %>% 
  summarise(N = n())
loo.compare.n <- min(focal.N$N)
loo.boot.n <- 1000

# sum dloo by sampling 75 data points without replacement for n times
dloo.boot.func <- function(x, boot.n, sample.size) {
  dloo.boot = numeric(boot.n)
  for (i in seq_len(boot.n)) {
    dloo.boot[i] = sum(sample(x, loo.compare.n, replace = TRUE))
  }
  out = c(median = median(dloo.boot),
          lci50 = as.numeric(quantile(dloo.boot, probs = 0.25)),
          uci50 = as.numeric(quantile(dloo.boot, probs = 0.75)),
          lci95 = as.numeric(quantile(dloo.boot, probs = 0.025)),
          uci95 = as.numeric(quantile(dloo.boot, probs = 0.975)))
  return(out)
}

dloo.boot.all <-
  dloo.pointwise %>%
  mutate(Focal = "All") %>%
  group_by(Focal) %>% 
  summarise(median.alpha = dloo.boot.func(dloo.alpha, 1000, loo.compare.n)["median"],
            lci50.alpha = dloo.boot.func(dloo.alpha, 1000, loo.compare.n)["lci50"],
            uci50.alpha = dloo.boot.func(dloo.alpha, 1000, loo.compare.n)["uci50"],
            lci95.alpha = dloo.boot.func(dloo.alpha, 1000, loo.compare.n)["lci95"],
            uci95.alpha = dloo.boot.func(dloo.alpha, 1000, loo.compare.n)["uci95"],
            median.beta = dloo.boot.func(dloo.beta, 1000, loo.compare.n)["median"],
            lci50.beta = dloo.boot.func(dloo.beta, 1000, loo.compare.n)["lci50"],
            uci50.beta = dloo.boot.func(dloo.beta, 1000, loo.compare.n)["uci50"],
            lci95.beta = dloo.boot.func(dloo.beta, 1000, loo.compare.n)["lci95"],
            uci95.beta = dloo.boot.func(dloo.beta, 1000, loo.compare.n)["uci95"])
dloo.boot.focal <- 
  dloo.pointwise %>% 
  group_by(Focal) %>% 
  summarise(median.alpha = dloo.boot.func(dloo.alpha, 1000, loo.compare.n)["median"],
            lci50.alpha = dloo.boot.func(dloo.alpha, 1000, loo.compare.n)["lci50"],
            uci50.alpha = dloo.boot.func(dloo.alpha, 1000, loo.compare.n)["uci50"],
            lci95.alpha = dloo.boot.func(dloo.alpha, 1000, loo.compare.n)["lci95"],
            uci95.alpha = dloo.boot.func(dloo.alpha, 1000, loo.compare.n)["uci95"],
            median.beta = dloo.boot.func(dloo.beta, 1000, loo.compare.n)["median"],
            lci50.beta = dloo.boot.func(dloo.beta, 1000, loo.compare.n)["lci50"],
            uci50.beta = dloo.boot.func(dloo.beta, 1000, loo.compare.n)["uci50"],
            lci95.beta = dloo.boot.func(dloo.beta, 1000, loo.compare.n)["lci95"],
            uci95.beta = dloo.boot.func(dloo.beta, 1000, loo.compare.n)["uci95"])
dloo.boot <- 
  bind_rows(dloo.boot.all, dloo.boot.focal) %>% 
  pivot_longer(cols = -Focal,
               names_to = "Summary",
               values_to = "dloo") %>% 
  separate(Summary, c("Summary", "Model"), "\\.") %>% 
  pivot_wider(names_from = Summary,
              values_from = dloo)
dloo.boot$Focal <-
  factor(dloo.boot$Focal, 
         levels = dloo.actual$Focal[order(dloo.actual$dloo.beta)])
