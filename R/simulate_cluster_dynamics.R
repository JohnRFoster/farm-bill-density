## functions to simulate pig population abundance

simulate_cluster_dynamics <- function(
    property_data,
    property_num,
    county_num,
    phi_mu,
    psi_phi,
    X,
    beta_p,
    start_density,
    method_lookup,
    effort_csv
){

  require(dplyr)
  require(tidyr)

  # mean litters per year from VerCauteren et al. 2019 pg 64
  data_litters_per_year <- c(1, 2, 0.86, 1, 2.28, 2.9, 0.49, 0.85, 1.57)

  # mean litter size year from VerCauteren et al. 2019 pg 63
  data_litter_size <- c(5.6, 6.1, 5.6, 6.1, 4.2, 5.0, 5.0, 6.5, 5.5, 6.8,
                        5.6, 5.9, 4.9, 5.1, 4.5, 4.7, 5.3, 5.7, 7.4, 8.4,
                        4.7, 4.9, 3.0, 3.0, 4.8, 4.8, 4.2, 5.4, 4.7, 5.2, 5.4)

  data_litter_size <- round(data_litter_size)

  # ecological process
  process_model <- function(N, zeta, a_phi, b_phi){
    phi <- rbeta(1, a_phi, b_phi)
    lambda <- N * zeta / 2 + N * phi
    rpois(1, lambda)
  }
  a_phi <- phi_mu * psi_phi
  b_phi <- (1 - phi_mu) * psi_phi

  # mean pigs produced per 28 days assuming 1 litter per year
  zeta <- mean(data_litter_size)*28*1/365

  survey_area <- property_data$area
  log_survey_area <- log(survey_area)
  log_rho <- log(method_lookup$rho)
  log_gamma <- log(method_lookup$gamma)
  p_unique <- method_lookup$p_unique

  source("R/functions_removal.R")
  effort_data <- readr::read_csv(effort_csv) |>
    suppressMessages()

  removal_effort <- property_data$effort
  sample_occasions <- removal_effort$sample_occasions

  # storage
  time_vec <- min(sample_occasions):max(sample_occasions)
  n_time <- length(time_vec)
  # pop_growth <- matrix(0, n_property, n_pp-1)
  take <- tibble()

  # spin-up period of 6 primary periods
  N_spin <- ceiling(survey_area * start_density) # initial abundance
  for(i in 1:6){
    N_spin <- process_model(N_spin, zeta, a_phi, b_phi)
  }

  extinct <- N_spin == 0
  too_dense <- N_spin / survey_area >= 10

  if(extinct | too_dense) return(NULL)

  N <- numeric(n_time)
  N[1] <- N_spin
  for(t in seq_len(n_time)){

    # determine if we sample during this primary period
    pp <- time_vec[t]
    rem <- pp %in% sample_occasions

    if(rem){
      sample_effort <- removal_effort |> filter(sample_occasions == pp)

      # determine order of removal events
      removal_order <- determine_removal_order(sample_effort)

      # conduct removals
      take_t <- conduct_removals(N[t], removal_order, effort_data, log_survey_area, X, beta_p, pp,
                                 log_rho, log_gamma, p_unique, method_lookup)
      take <- bind_rows(take, take_t)

      # how many pigs are left?
      z <- N[t] - sum(take_t$take)
      assertthat::assert_that(z >= 0,
                              msg = paste("More pigs removed than are alive! Time = ", t))
    } else {
      z <- N[t]
    }

    # simulate population dynamics
    if(t < n_time){
      N[t + 1] <- process_model(z, zeta, a_phi, b_phi)
    }
  }

  zero_take <- sum(take$take, na.rm = TRUE) == 0
  if(zero_take) return(NULL)

  N_tb <- tibble(
    PPNum = time_vec,
    N = N
  )

  # return all information
  # NA values for take information indicate no removal event occurred in that primary period
  all_info <- left_join(N_tb, take) |>
    group_by(PPNum) |>
    mutate(sum_take = sum(take)) |>
    ungroup() |>
    mutate(property = property_num,
           county = county_num,
           property_area = survey_area,
           obs_flag = if_else(is.na(take), 0, 1)) |>
    suppressMessages()

  pp_keep <- all_info |>
    select(PPNum, sum_take) |>
    distinct() |>
    mutate(sum_take = if_else(is.na(sum_take), 0, sum_take),
           cTake = cumsum(sum_take)) |>
    filter(cTake > 0) |>
    pull(PPNum)

  if(length(pp_keep) == 1){
    return(NULL) # need at least 2 primary periods for model to work
  } else {
    condition_first_capture <- all_info |>
      filter(PPNum %in% pp_keep)

    return(condition_first_capture)
  }

}
