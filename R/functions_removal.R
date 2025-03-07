determine_removal_order <- function(effort_df){
  # need to make effort_df long to remove methods not used in multimethod properties
  methods <- effort_df |>
    select(starts_with("method")) |>
    pivot_longer(cols = everything(),
                 names_to = "method_n",
                 values_to = "method",
                 values_drop_na = TRUE)
  reps <- effort_df |>
    select(starts_with("n_reps")) |>
    pivot_longer(cols = everything(),
                 names_to = "method_n",
                 values_to = "n_reps",
                 values_drop_na = TRUE)

  # randomly determine order of removal events
  tibble(
    method = rep(methods$method, times = reps$n_reps),
    order = sample.int(sum(reps$n_reps))
  ) |>
    arrange(order) |>
    pull(method)

}

generate_take_data <- function(m, effort_data){
  x <- effort_data |>
    filter(method == m)

  tc <- sample(x$trap_count, 1)
  y <- x |> filter(trap_count == tc)
  e <- sample(y$effort, 1)
  effort_per <- e / tc

  tibble(effort_per = effort_per, trap_count = tc)
}

conduct_removals <- function(N, removal_order, effort_data, log_survey_area, X, beta_p, pp,
                             log_rho, log_gamma, p_unique, method_lookup, alphs_project, alpha_cluster){

  n_reps <- length(removal_order)
  log_theta <- C <- rep(NA, n_reps)
  E <- tibble()

  for(r in seq_along(removal_order)){
    method_r <- removal_order[r]

    m <- method_lookup |>
      filter(method_name == method_r) |> pull(idx)

    e <- generate_take_data(m, effort_data)

    effort_per <- pull(e, effort_per)
    log_effort_per <- log(effort_per)
    trap_count <- pull(e, trap_count)
    n_trap_m1 <- trap_count - 1

    if(m <= 3){ # firearms, fixed wing and helicopter
      log_potential_area <- log_rho[m] + log_effort_per
    } else {
      log_potential_area <- log(pi) +
        (2 * (log_rho[m] + log_effort_per -
                log(exp(log_gamma[m-3]) + effort_per))) +
        log(1 + (p_unique[m-3] * n_trap_m1))
    }

    assertthat::assert_that(!is.na(log_potential_area),
                            msg = "NA log potential area!")

    # probability of capture, given that an individual is in the surveyed area
    log_theta[r] <- log(
      boot::inv.logit(
        beta_p[m, 1] +
          beta_p[m, 2] * X[1] +
          beta_p[m, 3] * X[2] +
          beta_p[m, 4] * X[3]) +
        alpha_project +
        alpha_cluster) +
      min(0, log_potential_area - log_survey_area)

    if(r == 1){
      p <- exp(log_theta[1])
      n_avail <- N
    } else {
      p <- exp(log_theta[1] + sum(log(1 - exp(log_theta[1:(r-1)]))))
      n_avail <- N - sum(C[1:(r-1)])
    }

    # remove
    C[r] <- min(rpois(1, n_avail * p), n_avail)

    et <- tibble(
      PPNum = pp,
      take = C[r],
      method = method_r,
      effort_per = effort_per,
      trap_count = trap_count,
      order = r,
      theta = exp(log_theta[r]),
      p = p,
      potential_area = exp(log_potential_area)
    )
    E <- bind_rows(E, et)

  }
  E
}


