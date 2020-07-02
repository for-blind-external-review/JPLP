# simulating PLP - time truncated case
sim_plp_tau = function(tau = 30,
                       beta = 1.5,
                       theta = 10){
  # initialization
  s = 0; t = 0
  while (max(t) <= tau) {
    u <- runif(1)
    s <- s - log(u)
    t_new <- theta*s^(1/beta)
    t <- c(t, t_new)
  }
  t = t[c(-1, -length(t))]

  return(t)
}

# simulating PLP - failure truncated case
sim_plp_n = function(mean_n, beta, theta){
  N = rpois(1, mean_n)
  u = runif(N, 0, 1)
  n_logu = -log(u)
  s = cumsum(n_logu)
  Delta_t = theta*s^(1/beta)
  return(Delta_t)
}

# simulate multiple NHPPs - failure truncated case
sim_mul_plp_n = function(n_shift = 20,
                         shift_len_mean = 20, shift_len_sd = 5,
                         theta = 10, beta = 2, mean_n = 5){
  #end_time = rnorm(n_shift, shift_len_mean, shift_len_sd)#to be deleted
  t_list = list()
  len_list = list()
  end_time1 = list()

  for (i in 1:n_shift) {
    t_list[[i]] = sim_plp_n(mean_n, beta, theta)
    len_list[[i]] = length(t_list[[i]])
    end_time1[[i]] = ifelse(length(t_list[[i]]) == 0, 0,
                            max(t_list[[i]]))
  }

  event_dat = data.frame(
    shift_id = rep(1:n_shift, unlist(len_list)),
    event_time = Reduce(c, t_list)
  )

  start_end_dat = data.frame(
    shift_id = 1:n_shift,
    start_time = rep(0, n_shift),
    end_time = Reduce(c, end_time1)#to be changed
  )
  start_end_dat$end_time[start_end_dat$end_time == 0] =
    mean(start_end_dat$end_time) # deal with no events

  return(list(event_dat = event_dat,
              start_end_dat = start_end_dat,
              shift_length = unlist(len_list)))
}

# simulate multiple NHPPs - time truncated case
sim_mul_plp_tau = function(n_shift = 20,
                           shift_len_mean = 20, shift_len_sd = 5,
                           theta = 10, beta = 2, mean_n = 5){
  tau_vector = rnorm(n_shift, shift_len_mean, shift_len_sd)#difference1

  t_list = list()
  len_list = list()
  # end_time1 = list() # not needed for time truncated case


  for (i in 1:n_shift) {
    t_list[[i]] = sim_plp_tau(tau_vector[i], beta, theta)
    len_list[[i]] = length(t_list[[i]])
  }

  event_dat = data.frame(
    shift_id = rep(1:n_shift, unlist(len_list)),
    event_time = Reduce(c, t_list)
  )

  start_end_dat = data.frame(
    shift_id = 1:n_shift,
    start_time = rep(0, n_shift),
    end_time = tau_vector #difference2
  )

  return(list(event_dat = event_dat,
              start_end_dat = start_end_dat,
              shift_length = unlist(len_list)))
}

plot_events = function(event_dat, start_end_dat, cross_size = 2){
  p = event_dat %>%
    ggplot(aes(x = event_time, y = shift_id)) +
    geom_point(alpha = 0.8, shape = 4, color = 'red', size = cross_size) +
    scale_y_continuous("shift ID",
                       labels = as.character(start_end_dat$shift_id),
                       breaks = start_end_dat$shift_id)+
    xlab('Time to event (minutes)') +
    geom_segment(data = start_end_dat,
                 aes(x = start_time, xend = end_time,
                     y = shift_id, yend = shift_id),
                 lineend = 'butt',
                 arrow = arrow(length = unit(0.2, "cm"))) +
    theme_classic()
  return(p)
}


sim_hier_plp_tau = function(N, beta = 1.5, theta){
  t_list = list()
  len_list = list()
  tau_vector = rnorm(N, 10, 1.3)

  for (i in 1:N) {
    t_list[[i]] = sim_plp_tau(tau_vector[i], beta = beta, theta = theta[i])
    len_list[[i]] = length(t_list[[i]])
  }

  event_dat = data.frame(
    shift_id = rep(1:N, unlist(len_list)),
    event_time = Reduce(c, t_list)
  )

  start_end_dat = data.frame(
    shift_id = 1:N,
    start_time = rep(0, N),
    end_time = tau_vector #difference2
  )

  return(list(event_dat = event_dat,
              start_end_dat = start_end_dat,
              shift_length = unlist(len_list)))
}


sim_hier_nhpp = function(
  beta = 1.5, # Shape parameter for PLP
  D = 10, # the number of drivers
  K = 3, # the number of predictor variables
  group_size_lambda = 10,
  mu0 = 0.2, # Hyperparameters: mean
  sigma0 = 0.5, # Hyperparameters: s.e.
  R_K = c(1, 0.3, 0.2)# 2. Fixed-effects parameters
)
{
  # 1. Random-effect intercepts
  r_0D = rnorm(D, mean = mu0, sd = sigma0)

  # 3. The number of shifts in the $d$-th driver: $N_{d}$
  N_K = rpois(D, group_size_lambda)
  N = sum(N_K) # the total number of obs
  id = rep(1:D, N_K)

  # 4. Generate data: x_1, x_2, .. x_K
  sim1 = function(group_sizes = N_K)
  {
    ntot = sum(group_sizes)

    int1 = rep(1, ntot)
    x1 = rnorm(ntot, 1, 1)
    x2 = rgamma(ntot, 1, 1)
    x3 = rpois(ntot, 2)

    return(data.frame(int1, x1, x2, x3))
  }
  X = sim1(N_K)

  # 5. Scale parameters of a NHPP
  # 5a. parameter matrix: P
  P = cbind(r0 = rep(r_0D, N_K),
            t(replicate(N, R_K)))
  M_logtheta = P*X

  # returned parameter for each observed shift
  theta_vec = exp(rowSums(M_logtheta))

  df = sim_hier_plp_tau(N = N, beta = beta, theta = theta_vec)

  hier_dat = list(
    N = nrow(df$event_dat),
    K = K,
    S = nrow(df$start_end_dat),
    D = max(id),
    id = id, # driver index at shift level
    tau = df$start_end_dat$end_time,
    event_time = df$event_dat$event_time,
    group_size = df$shift_length, #the number of events in each shift
    X_predictors = X[,2:4]
  )

  true_params = list(
    mu0 = mu0, sigma0 = sigma0,
    r0 = r_0D, r1_rk = R_K,
    beta = beta,
    theta = theta_vec
  )

  return(list(hier_dat = hier_dat, true_params = true_params))
}

pull_use = function(var = "theta", est_obj = f){
  z = est_obj %>%
    broom::tidy() %>%
    filter(grepl(var, term))
  return(z)
}


plot_est = function(data, var = "beta", hline_var = 1.5){
  p = data %>%
    filter(term == var) %>%
    ggplot(aes(id, est_mean)) +
    geom_point() +
    geom_line(linetype = "dashed", color = "red")+
    geom_errorbar(aes(ymax = est_mean + 1.96*est_sd,
                      ymin = est_mean - 1.96*est_sd),
                  width = 1)+
    geom_segment(aes(x = 10, xend = 100,
                     y = hline_var, yend = hline_var),
                 color = "green")+
    scale_x_continuous(breaks = c(0, 10, 25, 50, 75, 100),
                       labels = c("0", "10", "25", "50", "75", "100")) +
    labs(x = "The number of drivers (random effects)",
         y = var) +
    theme_bw()
  return(p)
}