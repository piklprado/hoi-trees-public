# Scripts to reproduce results in 
# Lai, Chong, Yee, Mayfield and Stouffer (2021) Ecology
# https://doi.org/10.1101/2020.09.16.300616 
# the Bayesian inference, figures, and diameter simulations will take a while 
# please have copious amount of hot beverage in hand

library(tidyverse)
library(brms)
library(future)
plan(multiprocess)
library(tidybayes)
library(ggpubr)


# Read data and fit model -------------------------------------------------
source("code/fit_models.R")


# Model comparison --------------------------------------------------------
source("code/compare_models.R")


# Figure 1 ----------------------------------------------------------------
source("code/fig1.R")


# Figure 2 ----------------------------------------------------------------
source("code/fig2.R")


# Figure 3 ----------------------------------------------------------------
source("code/fig3.R")


# Figure 4 ----------------------------------------------------------------
source("code/fig4.R")


# Figure 5 ----------------------------------------------------------------
source("code/fig5.R")
