
prep_nimble <- function(take){

  require(dplyr)
  require(tidyr)

  source("R/functions_misc.R")

  all_time_ids <- make_all_pp(take)

  n_clusters <- max(all_time_ids$cluster)
  n_properties <- max(all_time_ids$property)

  n_time_prop <- all_time_ids |>
    group_by(property) |>
    mutate(timestep = 1:n()) |>
    filter(timestep == max(timestep)) |>
    pull(timestep)

  assertthat::are_equal(length(n_time_prop), n_properties)

  n_time_clust <- all_time_ids |>
    select(cluster, PPNum) |>
    distinct() |>
    group_by(cluster) |>
    mutate(timestep = 1:n()) |>
    filter(timestep == max(timestep)) |>
    pull(timestep)

  assertthat::are_equal(length(n_time_clust), n_clusters)

  mH <- all_time_ids |>
    select(cluster, PPNum, m_id) |>
    distinct() |>
    group_by(cluster) |>
    mutate(timestep = 1:n()) |>
    ungroup() |>
    select(-PPNum) |>
    pivot_wider(names_from = timestep,
                values_from = m_id) |>
    select(-cluster)

  assertthat::are_equal(nrow(mH), n_clusters)

  nH <- all_time_ids |>
    group_by(property) |>
    mutate(timestep = 1:n()) |>
    ungroup() |>
    select(property, timestep, n_id) |>
    pivot_wider(names_from = timestep,
                values_from = n_id)

  assertthat::are_equal(nrow(nH), n_properties)

  nmH <- all_time_ids |>
    group_by(property) |>
    mutate(timestep = 1:n()) |>
    ungroup() |>
    select(property, timestep, m_id) |>
    pivot_wider(names_from = timestep,
                values_from = m_id)

  assertthat::are_equal(max(nmH, na.rm = TRUE), max(mH, na.rm = TRUE))
  assertthat::are_equal(nrow(nmH), n_properties)

  n_props_per_cluster <- take |>
    select(cluster, property) |>
    distinct() |>
    group_by(cluster) |>
    mutate(n = 1:n()) |>
    filter(n == max(n)) |>
    ungroup()

  solo_properties <- n_props_per_cluster |>
    filter(n == 1) |>
    pull(property)

  n_solo_properties <- length(solo_properties)

  actual_clusters <- n_props_per_cluster |>
    filter(n > 1) |>
    pull(cluster)

  grouped_properties <- take |>
    filter(cluster %in% actual_clusters) |>
    pull(property) |>
    unique()

  n_grouped_properties <- length(grouped_properties)

  assertthat::are_equal(n_solo_properties + n_grouped_properties, n_properties)

  n_time_single_property_clusters <- all_time_ids |>
    filter(property %in% solo_properties) |>
    group_by(property) |>
    mutate(timestep = 1:n()) |>
    filter(timestep == max(timestep)) |>
    pull(timestep)

  nmH_single <- nmH |>
    filter(property %in% solo_properties) |>
    select(-property)

  nH_single <- nH |>
    filter(property %in% solo_properties) |>
    select(-property)

  n_time_multi_property_clusters <- all_time_ids |>
    filter(property %in% grouped_properties) |>
    group_by(property) |>
    mutate(timestep = 1:n()) |>
    filter(timestep == max(timestep)) |>
    pull(timestep)

  nmH_multi <- nmH |>
    filter(property %in% grouped_properties) |>
    select(-property)

  nH_multi <- nH |>
    filter(property %in% grouped_properties) |>
    select(-property)

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

  removed_timestep <- all_time_ids |>
    select(cluster, PPNum, m_id) |>
    distinct() |>
    group_by(cluster) |>
    mutate(timestep = 1:n()) |>
    ungroup() |>
    left_join(sum_take) |>
    mutate(sum_take = if_else(is.na(sum_take), 0, sum_take)) |>
    select(cluster, timestep, sum_take) |>
    pivot_wider(names_from = timestep,
                values_from = sum_take) |>
    select(-cluster)

  assertthat::are_equal(dim(removed_timestep), dim(mH))

  tH <- take |>
    select(property, PPNum)

  nH_p <- left_join(tH, all_time_ids) |> pull(n_id)

  assertthat::are_equal(length(nH_p),  nrow(take))

  n_observed_units <- take |>
    select(property, PPNum) |>
    distinct() |>
    nrow()

  assertthat::are_equal(length(unique(nH_p)), n_observed_units)

  cluster_areas <- take |>
    filter(!property %in% solo_properties) |>
    select(property, cluster_area_km2) |>
    distinct() |>
    pull(cluster_area_km2)

  assertthat::are_equal(length(cluster_areas), n_grouped_properties)

  property_areas <- take |>
    filter(!property %in% solo_properties) |>
    select(property, property_area_km2) |>
    distinct() |>
    pull(property_area_km2)

  assertthat::are_equal(length(property_areas), n_grouped_properties)

  # mean litter size year from VerCauteren et al. 2019 pg 63
  data_litter_size <- c(5.6, 6.1, 5.6, 6.1, 4.2, 5.0, 5.0, 6.5, 5.5, 6.8,
                        5.6, 5.9, 4.9, 5.1, 4.5, 4.7, 5.3, 5.7, 7.4, 8.4,
                        4.7, 4.9, 3.0, 3.0, 4.8, 4.8, 4.2, 5.4, 4.7, 5.2, 5.4)

  data_litter_size <- round(data_litter_size)

  X <- take |>
    select(c_road_den, c_rugged, c_canopy) |>
    as.matrix()

  n_method <- length(unique(take$method))

  constants <- list(
    n_cluster = n_clusters,
    n_survey = nrow(take),
    n_ls = length(data_litter_size),
    n_method = n_method,
    n_time_clust = n_time_clust,
    n_first_survey = length(which(take$order == 1)),
    n_not_first_survey = length(which(take$order != 1)),
    n_single_property_clusters = n_solo_properties,
    n_multi_property_clusters = n_grouped_properties,
    n_time_single_property_clusters = n_time_single_property_clusters,
    n_time_multi_property_clusters = n_time_multi_property_clusters,
    mH = as.matrix(mH),
    nH_single = as.matrix(nH_single),
    nH_multi = as.matrix(nH_multi),
    nmH_single = as.matrix(nmH_single),
    nmH_multi = as.matrix(nmH_multi),
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
    log_nu_tau = 1,
    n_betaP = n_method * ncol(X),
    beta_p_row = rep(1:n_method, each = ncol(X)),
    beta_p_col = rep(1:ncol(X), n_method)

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

  constants_nimble$phi_mu <- draw_value("phi_mu")
  constants_nimble$psi_phi <- draw_value("psi_phi")
  nu <- draw_value("nu")
  constants_nimble$zeta <- 28 * nu / 365
  constants_nimble$beta_p <- matrix(draw_value("beta_p"), 5, 3, byrow = TRUE)
  constants_nimble$beta1 <- draw_value("beta1")
  constants_nimble$omega <- draw_value("omega")
  constants_nimble$gamma <- draw_value("gamma")
  constants_nimble$rho <- draw_value("rho")

  # for testing
  # attach(constants_nimble)
  # attach(data_nimble)

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
        # if(is.na(M[mH[i, j]])) print(j); print(M[mH[i, j]])

      }
    }

    maxn1 <- max(nH_single, na.rm = TRUE)
    maxn2 <- max(nH_multi, na.rm = TRUE)
    N <- rep(NA, max(c(maxn1, maxn2)))

    for(i in 1:n_single_property_clusters){
      for(t in 1:n_time_single_property_clusters[i]){
        N[nH_single[i, t]] <- M[nmH_single[i, t]]
      }
    }

    # multiple property clusters
    # distribute based on average density
    for(i in 1:n_multi_property_clusters){
      for(t in 1:n_time_multi_property_clusters[i]){
        d <- log(M[nmH_multi[i, t]]) - log_cluster_area[i] + log_property_area[i]
        N[nH_multi[i, t]] <- rpois(1, exp(d))
      }
    }


    list(
      log_lambda_1 = log(M_init + buffer),
      beta_p = beta_p,
      beta1 = beta1,
      p_mu = boot::logit(omega),
      p_unique = omega,
      phi_mu = phi_mu,
      psi_phi = psi_phi,
      a_phi = a,
      b_phi = b,
      N = N + buffer,
      M = M + buffer,
      # lambda = lambda + buffer,
      log_nu = log(nu),
      nu = nu,
      log_gamma = log(gamma),
      log_rho = log(rho),
      phi = phi,
      zeta = zeta,
      log_zeta = log(zeta)
    )

  })
}






