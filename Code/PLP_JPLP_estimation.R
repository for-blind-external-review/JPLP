library(dplyr)
library(rstan)
library(broom)

# Read in data
PLP = readRDS('Data/PLP.rds')
JPLP = readRDS('Data/JPLP.rds')


# Hamiltonian Monte Carlo to sample from the posterior distribution
fit_PLP = stan("nhpp_plp_hierarchical.stan",
               chains = 4, iter = 5000, warmup = 1000, data = PLP)
fit_JPLP = stan("jplp_hierarchical.stan",
                chains = 4, iter = 5000, warmup = 1000, data = JPLP)


# Get posterior parameter estimate statistics
tidy(fit_NHPP, conf.int = T, rhat = T, ess = T)
tidy(fit_JPLP, conf.int = T, rhat = T, ess = T)

# Trace plots
stan_trace(fit_NHPP, pars = c('mu0', 'sigma0', 'beta'), ncol = 1)
stan_trace(fit_JPLP, pars = c('mu0', 'sigma0', 'beta', 'kappa'), ncol = 1)