modulus_trans <- function(g, lambda) {
  ifelse(g >= 0, g^lambda, -((- g)^lambda))
}

inv_modulus_trans <- function(g_t, lambda) {
  ifelse(g_t >= 0, g_t^(1/lambda), -((- g_t)^(1/lambda)))
}

# lambda <- seq(0.1, 0.8, 0.01)
# lambda <- seq(-2, 2, 0.01)

# car::yjPower(out$G, lambda = 0.55)

# skewness_test <- sapply(lambda, function(x) {
#   g_trans <- modulus_trans(g = out$G, lambda = x)
#   # g_trans <- car::yjPower(out$G, lambda = x)
#   moments::skewness(g_trans)
# })
#   
# plot(lambda, skewness_test)
# abline(h = 0)
# 
# lambda_best <- lambda[which(abs(skewness_test) == min(abs(skewness_test), na.rm = TRUE))]
# 
# G_modulus <- modulus_trans(g = out$G, lambda = lambda_best)
# G_yjPower <- car::yjPower(out$G, lambda = lambda_best)
# 
# plot(out$G, G_modulus)
# plot(out$G, G_yjPower)
# plot(G_modulus, G_yjPower)
# hist(G_modulus, breaks = 200)
# hist(G_yjPower, breaks = 200)
# hist(out$G, breaks = 200)
# 
# dat <- data.frame(G_modulus = G_modulus, 
#                   G_yjPower = G_yjPower, 
#                   G = out$G,
#                   sp = out$Focal)
# dat2 <- dat[dat$G >= 0, ]
# 
# hist(dat2$G_modulus, breaks = 100)
# abline(v = c(0.1, 0.2, 0.3, 0.4))
# 
# hist(dat2$G_yjPower, breaks = 100)
# 
# 
# library(brms)
# 
# m0 <- brm(G ~ 1, data = dat, family = gaussian(), core = 4)
# m1 <- brm(G ~ 1, data = dat, family = student(), core = 4)
# m2 <- brm(G_modulus ~ 1, data = dat, family = gaussian(), core = 4)
# m3 <- brm(G_yjPower ~ 1, data = dat, family = gaussian(), core = 4)
# m3 <- brm(G_yjPower ~ 1, data = dat, family = gaussian(), core = 4)
# 
# mcmc_dens(m0, "b_Intercept")
# mcmc_dens(m1, "b_Intercept")
# 
# brms::pp_check(m0)
# brms::pp_check(m1) + lims(x = c(-2, 2))
# brms::pp_check(m2)
# brms::pp_check(m3)
