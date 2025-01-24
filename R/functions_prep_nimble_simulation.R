
prep_nimble <- function(take){

  require(dplyr)
  require(tidyr)

  source("R/functions_misc.R")

  all_cluster_pp <- make_all_pp(take, "cluster") |> mutate(m_id = 1:n())
  all_prop_pp <- make_all_pp(take, "property") |> mutate(n_id = 1:n())

  cluster_property_lookup <- take |>
    select(cluster, property) |>
    distinct()

  all_prop_pp2 <- left_join(all_prop_pp, cluster_property_lookup)
  all_cluster_pp2 <- all_cluster_pp |> rename(timestep_m = timestep)

  all_timepoints <- left_join(all_prop_pp2, all_cluster_pp2) |>
    arrange(property, timestep)

  # START WITH ALL TIMEPOINTS SO THAT EVERYTIME THERE IS A CLUSTER ESTIMATE
  # ALL PROPERTIES WITHIN THAT CLUSTER GET AN ESTIMATE

  # NEED TO DEFINE
  # n_single_property_clusters
  # n_time_single_property_clusters
  # n_multi_property_clusters
  # n_time_multi_property_clusters

  n_clusters <- max(all_cluster_pp$cluster)
  n_properties <- max(all_prop_pp$property)

  n_time_prop <- all_prop_pp |>
    group_by(property) |>
    filter(timestep == max(timestep)) |>
    pull(timestep)

  assertthat::are_equal(length(n_time_prop), n_properties)

  n_time_clust <- all_cluster_pp |>
    group_by(cluster) |>
    filter(timestep == max(timestep)) |>
    pull(timestep)

  assertthat::are_equal(length(n_time_clust), n_clusters)

  mH <- all_cluster_pp |>
    select(-PPNum) |>
    pivot_wider(names_from = timestep,
                values_from = m_id) |>
    select(-cluster)

  assertthat::are_equal(nrow(mH), n_clusters)

  nH <- all_prop_pp |>
    select(-PPNum) |>
    pivot_wider(names_from = timestep,
                values_from = n_id) |>
    select(-property)

  assertthat::are_equal(nrow(nH), n_properties)



  clust_prop_join <- left_join(all_prop_pp2, all_cluster_pp2)
  nmH <- clust_prop_join |>
    select(property, timestep, m_id) |>
    pivot_wider(names_from = timestep,
                values_from = m_id) |>
    select(-property)

  assertthat::are_equal(max(nmH, na.rm = TRUE), max(mH, na.rm = TRUE))
  assertthat::are_equal(nrow(nmH), n_properties)


  take_timestep <- take |>
    select(property, PPNum) |>
    distinct() |>
    group_by(property) |>
    mutate(timestep = 1:n()) |>
    ungroup()

  take <- left_join(take, take_timestep) |>
    mutate(start = 0,
           end = 0)

  pb <- txtProgressBar(max = nrow(take), style = 1)
  for (i in 1:nrow(take)) {
    if (take$order[i] > 1) {
      idx <- which(#take$county == take$county[i] &
                     take$property == take$property[i] &
                     take$timestep == take$timestep[i] &
                     take$order < take$order[i])
      take$start[i] <- idx[1]
      take$end[i] <- idx[length(idx)]
      assertthat::are_equal(idx, take$start[i]:take$end[i])
    }
    setTxtProgressBar(pb, i)
  }
  close(pb)

  y_sum <- take |>
    group_by(cluster, property, PPNum) |>
    mutate(ysum = cumsum(take) - take) |>
    ungroup() |>
    pull(ysum)

  assertthat::are_equal(length(y_sum), nrow(take))

  sum_take <- take |>
    group_by(cluster, PPNum) |>
    summarise(sum_take = sum(take)) |>
    ungroup()

  removed_timestep <- left_join(all_cluster_pp, sum_take) |>
    mutate(sum_take = if_else(is.na(sum_take), 0, sum_take)) |>
    select(-PPNum, -m_id) |>
    pivot_wider(names_from = timestep,
                values_from = sum_take) |>
    select(-cluster)

  assertthat::are_equal(dim(removed_timestep), dim(mH))

  tH <- take |>
    select(property, PPNum)

  nH_p <- left_join(tH, all_prop_pp) |> pull(n_id)

  cluster_areas <- take |>
    select(property, cluster_area_km2) |>
    distinct() |>
    pull(cluster_area_km2)

  assertthat::are_equal(length(cluster_areas), n_properties)

  property_areas <- take |>
    select(property, property_area_km2) |>
    distinct() |>
    pull(property_area_km2)

  assertthat::are_equal(length(property_areas), n_properties)

  # mean litter size year from VerCauteren et al. 2019 pg 63
  data_litter_size <- c(5.6, 6.1, 5.6, 6.1, 4.2, 5.0, 5.0, 6.5, 5.5, 6.8,
                        5.6, 5.9, 4.9, 5.1, 4.5, 4.7, 5.3, 5.7, 7.4, 8.4,
                        4.7, 4.9, 3.0, 3.0, 4.8, 4.8, 4.2, 5.4, 4.7, 5.2, 5.4)

  data_litter_size <- round(data_litter_size)

  X <- take |>
    select(c_road_den, c_rugged, c_canopy) |>
    as.matrix()

  constants <- list(
    n_cluster = n_clusters,
    n_property = n_properties,
    n_survey = nrow(take),
    n_ls = length(data_litter_size),
    n_method = length(unique(take$method)),
    n_time_prop = n_time_prop,
    n_time_clust = n_time_clust,
    n_first_survey = length(which(take$order == 1)),
    n_not_first_survey = length(which(take$order != 1)),
    mH = as.matrix(mH),
    nH = as.matrix(nH),
    nmH = as.matrix(nmH),
    start = take$start,
    end = take$end,
    y_sum = y_sum,
    rem = as.matrix(removed_timestep),
    first_survey = which(take$order == 1),
    not_first_survey = which(take$order != 1),
    nH_p = nH_p,
    log_pi = log(pi),
    m_p = ncol(X),
    method = as.numeric(as.factor(take$method)),
    property = take$property,
    pp_len = 28,
    phi_mu_a = 3.23,
    phi_mu_b = 0.2,
    log_rho_mu = rep(0, 5),
    log_rho_tau = rep(0.1, 5),
    p_mu_mu = rep(0, 2),
    p_mu_tau = rep(1, 2),
    log_gamma_mu = rep(0, 2),
    log_gamma_tau = rep(0.1, 2),
    beta1_mu = rep(0, 5),
    beta1_tau = rep(1, 5),
    beta_p_mu = rep(0, 15),
    beta_p_tau = rep(1, 15),
    psi_shape = 1,
    psi_rate = 0.1,
    log_nu_mu = 2,
    log_nu_tau = 1

  )

  data <- list(
    y = take$take,
    J = data_litter_size,
    X_p = X,
    effort_per = take$effort_per,
    log_effort_per = log(take$effort_per),
    n_trap_m1 = take$trap_count - 1,
    log_survey_area_km2 = log(take$property_area_km2),
    log_cluster_area = log(cluster_areas),
    log_property_area = log(property_areas)
  )

  return(
    list(
      constants = constants,
      data = data
    )
  )
}

nimble_inits <- function(constants_nimble, data_nimble, start_density, buffer = 1000){

  params <- read_csv("../pigs-property/data/posterior_95CI_range_all.csv") |>
    suppressMessages()

  draw_value <- function(x){
    params |>
      filter(grepl(x, node)) |>
      pull(med) |>
      round(2) |>
      jitter()
  }

  phi_mu <- draw_value("phi_mu")
  psi_phi <- draw_value("psi_phi")
  nu <- draw_value("nu")
  zeta <- 28 * nu / 365
  beta_p <- matrix(draw_value("beta_p"), 5, 3, byrow = TRUE)
  beta1 <- matrix(draw_value("beta1"), 5, 1)
  beta <- cbind(beta1, beta_p)
  omega <- draw_value("omega")
  gamma <- draw_value("gamma")
  log_gamma <- log(gamma)
  rho <- draw_value("rho")
  log_rho <- log(rho)

  with(append(constants_nimble, data_nimble), {

    a <- phi_mu * psi_phi
    b <- (1 - phi_mu) * psi_phi
    M <- phi <- lambda <- rep(NA, max(mH, na.rm = TRUE))
    M_init <- rep(NA, n_cluster)
    i=1
    for(i in 1:n_cluster){

      M_init[i] <- round(exp(log_cluster_area[i])*start_density + sum(rem[i, ], na.rm = TRUE))
      M[mH[i, 1]] <- M_init[i]
      j=2
      for(j in 2:n_time_clust[i]){

        phi[mH[i, j-1]] <- rbeta(1, a, b)
        z <- M[mH[i, j-1]] - rem[i, j-1]
        z <- max(1, z)

        lambda[mH[i, j-1]] <- z * zeta / 2 + z * phi[mH[i, j-1]]

        M[mH[i, j]] <- rpois(1, lambda[mH[i, j-1]])
        if(is.na(M[mH[i, j]])) print(j); print(M[mH[i, j]])

      }
    }

    list(
      log_lambda_1 = log(n_init + buffer),
      beta_p = beta_p,
      beta1 = beta1,
      p_mu = p_mu,
      p_unique = boot::inv.logit(p_mu),
      phi_mu = phi_mu,
      psi_phi = psi_phi,
      a_phi = a,
      b_phi = b,
      N = N + buffer,
      # lambda = lambda + buffer,
      log_nu = log(mean_ls),
      nu = mean_ls,
      log_gamma = log_gamma,
      log_rho = log_rho,
      phi = phi,
      zeta = zeta,
      log_zeta = log(zeta)
    )

  })
}






