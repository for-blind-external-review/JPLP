pacman::p_load(dplyr, rstan, broom)
source('Functions/NHPP_functions.R')
set.seed(123)
df_NHPP = sim_hier_nhpp(D = 10, beta = 1.2)
PLP = df_NHPP$hier_dat

source('Functions/JPLP_functions.R')
set.seed(123)
df_JPLP = sim_hier_JPLP(D = 10, beta = 1.2)
JPLP = df_JPLP$stan_dt

saveRDS(PLP, 'Data/PLP.rds')
saveRDS(JPLP, 'Data/JPLP.rds')







