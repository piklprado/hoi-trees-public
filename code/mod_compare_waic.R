# script to summarise WAIC

nullLL <- log_lik(null.mlm)
alphaLL <- log_lik(alpha.mlm)
betaLL <- log_lik(beta.mlm)
null.waic <- waic(nullLL)
alpha.waic <- waic(alphaLL)
beta.waic <- waic(betaLL)

# Pointwise WAIC
null.waic.pointwise <- apply(nullLL, 2, function(x) waic(as.matrix(x))$estimates["waic","Estimate"])
alpha.waic.pointwise <- apply(alphaLL, 2, function(x) waic(as.matrix(x))$estimates["waic","Estimate"])
beta.waic.pointwise <- apply(betaLL, 2, function(x) waic(as.matrix(x))$estimates["waic","Estimate"])
dWAIC.pointwise <- 
  beta.mlm$data %>% 
  select(Focal) %>% 
  mutate(null = null.waic.pointwise,
         alpha = alpha.waic.pointwise,
         beta  = beta.waic.pointwise) %>% 
  mutate(dWAIC.alpha = alpha - null,
         dWAIC.beta = beta - null)
dWAIC.actual.all <-
  dWAIC.pointwise %>% 
  mutate(Focal = "All") %>%
  group_by(Focal) %>% 
  summarise(dWAIC.alpha = sum(dWAIC.alpha),
            dWAIC.beta = sum(dWAIC.beta))
dWAIC.actual.focal <-
  dWAIC.pointwise %>% 
  group_by(Focal) %>% 
  summarise(dWAIC.alpha = sum(dWAIC.alpha),
            dWAIC.beta = sum(dWAIC.beta))
dWAIC.actual <- bind_rows(dWAIC.actual.all, dWAIC.actual.focal)

# Bootstrap pointwise WAIC
focal.N <- 
  dWAIC.pointwise %>% 
  group_by(Focal) %>% 
  summarise(N = n())
waic.compare.n <- min(focal.N$N)
waic.boot.n <- 1000

# sum dWAIC by sampling 75 data points without replacement for n times
dWAIC.boot.func <- function(x, boot.n, sample.size) {
  dWAIC.boot = numeric(boot.n)
  for (i in seq_len(boot.n)) {
    dWAIC.boot[i] = sum(sample(x, waic.compare.n, replace = TRUE))
  }
  out = c(median = median(dWAIC.boot),
          lci50 = as.numeric(quantile(dWAIC.boot, probs = 0.25)),
          uci50 = as.numeric(quantile(dWAIC.boot, probs = 0.75)),
          lci95 = as.numeric(quantile(dWAIC.boot, probs = 0.025)),
          uci95 = as.numeric(quantile(dWAIC.boot, probs = 0.975)))
  return(out)
}

dWAIC.boot.all <-
  dWAIC.pointwise %>%
  mutate(Focal = "All") %>%
  group_by(Focal) %>% 
  summarise(median.alpha = dWAIC.boot.func(dWAIC.alpha, 1000, waic.compare.n)["median"],
            lci50.alpha = dWAIC.boot.func(dWAIC.alpha, 1000, waic.compare.n)["lci50"],
            uci50.alpha = dWAIC.boot.func(dWAIC.alpha, 1000, waic.compare.n)["uci50"],
            lci95.alpha = dWAIC.boot.func(dWAIC.alpha, 1000, waic.compare.n)["lci95"],
            uci95.alpha = dWAIC.boot.func(dWAIC.alpha, 1000, waic.compare.n)["uci95"],
            median.beta = dWAIC.boot.func(dWAIC.beta, 1000, waic.compare.n)["median"],
            lci50.beta = dWAIC.boot.func(dWAIC.beta, 1000, waic.compare.n)["lci50"],
            uci50.beta = dWAIC.boot.func(dWAIC.beta, 1000, waic.compare.n)["uci50"],
            lci95.beta = dWAIC.boot.func(dWAIC.beta, 1000, waic.compare.n)["lci95"],
            uci95.beta = dWAIC.boot.func(dWAIC.beta, 1000, waic.compare.n)["uci95"])
dWAIC.boot.focal <- 
  dWAIC.pointwise %>% 
  group_by(Focal) %>% 
  summarise(median.alpha = dWAIC.boot.func(dWAIC.alpha, 1000, waic.compare.n)["median"],
            lci50.alpha = dWAIC.boot.func(dWAIC.alpha, 1000, waic.compare.n)["lci50"],
            uci50.alpha = dWAIC.boot.func(dWAIC.alpha, 1000, waic.compare.n)["uci50"],
            lci95.alpha = dWAIC.boot.func(dWAIC.alpha, 1000, waic.compare.n)["lci95"],
            uci95.alpha = dWAIC.boot.func(dWAIC.alpha, 1000, waic.compare.n)["uci95"],
            median.beta = dWAIC.boot.func(dWAIC.beta, 1000, waic.compare.n)["median"],
            lci50.beta = dWAIC.boot.func(dWAIC.beta, 1000, waic.compare.n)["lci50"],
            uci50.beta = dWAIC.boot.func(dWAIC.beta, 1000, waic.compare.n)["uci50"],
            lci95.beta = dWAIC.boot.func(dWAIC.beta, 1000, waic.compare.n)["lci95"],
            uci95.beta = dWAIC.boot.func(dWAIC.beta, 1000, waic.compare.n)["uci95"])
dWAIC.boot <- 
  bind_rows(dWAIC.boot.all, dWAIC.boot.focal) %>% 
  pivot_longer(cols = -Focal,
               names_to = "Summary",
               values_to = "dWAIC") %>% 
  separate(Summary, c("Summary", "Model"), "\\.") %>% 
  pivot_wider(names_from = Summary,
              values_from = dWAIC)
dWAIC.boot$Focal <-
  factor(dWAIC.boot$Focal, 
         levels = dWAIC.actual$Focal[order(dWAIC.actual$dWAIC.beta)])