# Script to fit models

library(tidyverse)
library(brms)
library(future)
plan(multiprocess)

# auxillary functions
source("code/generate_model_matrix.R")
source("code/modulus_transformation.R")


# Read data ---------------------------------------------------------------

mandai <- read_csv("data/mandai_clean.csv")

# some constant variables for the loop below
focals <- as.character(sort(unique(mandai$Focal)))
non.X.colnames <- c("Plot", "Year", "Tag.tree", "Focal", "DBH_mi", "G")

# manually generate model matrix
in.dat <- generate.model.matrix(mandai, scale.X = TRUE) 

# some transformation prior to modelling
in.dat$in.dat <- 
  in.dat$in.dat %>% 
  # scale DBH and apply modulus transformation to response
  mutate(DBH_mi = DBH_mi / sd(DBH_mi),
         G = modulus_trans(G, lambda = 0.55))


# Fit models --------------------------------------------------------------
# Priors
nlprior <-
  prior(normal(0, 10), class = "b", nlpar = "a") +
  prior(normal(0, 10), class = "b", nlpar = "b", lb = 0) +
  prior(normal(0, 10), class = "b", nlpar = "c", lb = 0) +
  prior(student_t(3, 0, 1), class = "sd", nlpar = "a") +
  prior(student_t(3, 0, 1), class = "sd", nlpar = "b") +
  prior(student_t(3, 0, 1), class = "sd", nlpar = "c")

# Null model without biotic interaction
null.mlm <-
  brm(
    bf(G ~ (DBH_mi ^ b) * exp(a - c * DBH_mi),
       a ~ 1 + (1 || Focal) + (1 || Plot) + (1 || Year),
       b ~ 1 + (1 || Focal),
       c ~ 1 + (1 || Focal),
       sigma ~ 1 + (1 || Focal),
       nl = TRUE),
    data = in.dat$in.dat,
    family = gaussian(),
    prior = nlprior,
    inits = 0,
    warmup = 3000,
    iter = 4000,
    cores = 4,
    control = list(adapt_delta = 0.99,
                   max_treedepth = 15)
  )

# Direct-interaction-only model
alpha.part <- paste(c(sort(focals), "OTHERS"), collapse = " + ")
a.alpha.formula <- 
  as.formula(paste0(
    "a ~ 1 + ", 
    paste0("(1 +", alpha.part, "|| Focal) + (1 || Plot) + (1 || Year)")))
alpha.mlm <-
  brm(
    bf(G ~ (DBH_mi ^ b) * exp(a - c * DBH_mi),
       a.alpha.formula,
       b ~ 1 + (1 || Focal),
       c ~ 1 + (1 || Focal),
       sigma ~ 1 + (1 || Focal),
       nl = TRUE),
    data = in.dat$in.dat,
    family = gaussian(),
    prior = nlprior,
    inits = 0,
    warmup = 3000,
    iter = 4000,
    cores = 4,
    control = list(adapt_delta = 0.99,
                   max_treedepth = 15)
  )

# HOI-inclusive model
beta.part <- 
  paste(setdiff(colnames(in.dat$in.dat), non.X.colnames), collapse = " + ")
a.beta.formula <-
  as.formula(paste0(
    "a ~ 1 + ",
    paste0("(1 +", beta.part, "|| Focal) + (1 || Plot) + (1 || Year)")))
beta.mlm <-
  brm(
    bf(G ~ (DBH_mi ^ b) * exp(a - c * DBH_mi), 
       a.beta.formula, 
       b ~ 1 + (1 || Focal),
       c ~ 1 + (1 || Focal),
       sigma ~ 1 + (1 || Focal),
       nl = TRUE),
    data = in.dat$in.dat,
    family = gaussian(),
    prior = nlprior,
    inits = 0,
    warmup = 3000,
    iter = 4000,
    cores = 4,
    control = list(adapt_delta = 0.99,
                   max_treedepth = 15)
  )
