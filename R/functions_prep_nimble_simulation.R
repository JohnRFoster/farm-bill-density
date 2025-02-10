
prep_nimble <- function(take, posterior_path){

  require(dplyr)
  require(tidyr)

  source("R/functions_misc.R")

  all_obs_ids <- take |>
    select(cluster, property, PPNum) |>
    mutate(mcmc_order = 1:n())

  all_time_ids <- make_all_pp(take)

  n_clusters <- max(all_obs_ids$cluster)
  n_properties <- max(all_obs_ids$property)

  cluster_timesteps <- all_time_ids |>
    select(cluster, PPNum) |>
    distinct() |>
    group_by(cluster) |>
    mutate(timestep = 1:n()) |>
    ungroup()

  assertthat::are_equal(nrow(cluster_timesteps), max(all_time_ids$m_id))

  n_time_clust <- cluster_timesteps |>
    group_by(cluster) |>
    filter(timestep == max(timestep)) |>
    pull(timestep)

  assertthat::assert_that(min(n_time_clust) >= 2)
  assertthat::are_equal(length(n_time_clust), n_clusters)

  mH <- all_time_ids |>
    left_join(cluster_timesteps) |>
    select(cluster, timestep, m_id) |>
    distinct() |>
    pivot_wider(names_from = timestep,
                values_from = m_id)

  assertthat::are_equal(mH$cluster, unique(all_obs_ids$cluster))
  assertthat::are_equal(nrow(mH), n_clusters)

  mH <- mH |> select(-cluster)

  property_timesteps <- all_time_ids |>
    select(property, PPNum) |>
    distinct() |>
    group_by(property) |>
    mutate(timestep = 1:n()) |>
    ungroup()

  nH <- all_time_ids |>
    left_join(property_timesteps) |>
    select(property, timestep, n_id) |>
    pivot_wider(names_from = timestep,
                values_from = n_id)

  assertthat::are_equal(nH$property, unique(all_obs_ids$property))
  assertthat::are_equal(nrow(nH), n_properties)

  nmH <- all_time_ids |>
    left_join(property_timesteps) |>
    select(property, timestep, m_id) |>
    pivot_wider(names_from = timestep,
                values_from = m_id)

  assertthat::are_equal(nmH$property, unique(all_obs_ids$property))
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
    left_join(property_timesteps) |>
    group_by(property) |>
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
                values_from = sum_take)

  assertthat::are_equal(removed_timestep$cluster, unique(all_obs_ids$cluster))

  removed_timestep <- removed_timestep |> select(-cluster)
  assertthat::are_equal(dim(removed_timestep), dim(mH))

  tH <- take |>
    select(property, PPNum)

  nH_p <- left_join(tH, all_time_ids)

  assertthat::are_equal(unique(nH_p$property), unique(all_obs_ids$property))
  assertthat::are_equal(unique(nH_p$cluster), unique(all_obs_ids$cluster))

  nH_p <- nH_p$n_id
  assertthat::are_equal(length(nH_p),  nrow(take))

  n_observed_units <- take |>
    select(property, PPNum) |>
    distinct() |>
    nrow()

  assertthat::are_equal(length(unique(nH_p)), n_observed_units)

  cluster_areas <- take |>
    filter(property %in% grouped_properties) |>
    select(property, cluster_area_km2) |>
    distinct() |>
    pull(cluster_area_km2)

  assertthat::are_equal(length(cluster_areas), n_grouped_properties)

  property_areas <- take |>
    filter(property %in% grouped_properties) |>
    select(property, property_area_km2) |>
    distinct() |>
    pull(property_area_km2)

  assertthat::are_equal(length(property_areas), n_grouped_properties)

  # mean litter size year from VerCauteren et al. 2019 pg 63
  data_litter_size <- c(5.6, 6.1, 5.6, 6.1, 4.2, 5.0, 5.0, 6.5, 5.5, 6.8,
                        5.6, 5.9, 4.9, 5.1, 4.5, 4.7, 5.3, 5.7, 7.4, 8.4,
                        4.7, 4.9, 3.0, 3.0, 4.8, 4.8, 4.2, 5.4, 4.7, 5.2, 5.4)

  data_litter_size <- round(data_litter_size)

  rds <- file.path(posterior_path, "posteriorSamples.rds")
  post <- read_rds(rds)

  get_vec <- function(df, node){
    df |>
      select(contains(node)) |>
      pivot_longer(cols = everything(),
                   names_to = "node") |>
      group_by(node) |>
      summarise(mu = mean(value),
                tau = 1/var(value))
  }

  log_rho <- get_vec(post, "log_rho")
  p_mu <- get_vec(post, "p_mu")
  log_gamma <- get_vec(post, "log_gamma")
  beta1 <- get_vec(post, "beta1")
  beta_p <- get_vec(post, "beta_p")
  log_nu <- get_vec(post, "log_nu")

  phi_mu <- get_vec(post, "phi_mu")
  mu <- phi_mu$mu
  v <- 1 / phi_mu$tau
  w <- ((mu * (1 - mu)) / v) - 1
  alpha <- mu * w
  beta <- (1 - mu) * w

  psi_phi <- get_vec(post, "psi_phi")
  mu <- psi_phi$mu
  v <- 1 / psi_phi$tau
  shape <- (mu^2) / v
  rate <- mu / v

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
    phi_mu_a = alpha,
    phi_mu_b = beta,
    log_rho_mu = log_rho$mu,
    log_rho_tau = log_rho$tau,
    p_mu_mu = p_mu$mu,
    p_mu_tau = p_mu$tau,
    log_gamma_mu = log_gamma$mu,
    log_gamma_tau = log_gamma$tau,
    beta1_mu = beta1$mu,
    beta1_tau = beta1$tau,
    beta_p_mu = beta_p$mu,
    beta_p_tau = beta_p$tau,
    psi_shape = shape,
    psi_rate = rate,
    log_nu_mu = log_nu$mu,
    log_nu_tau = log_nu$tau,
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

nimble_inits <- function(constants_nimble, data_nimble, start_density, buffer = 100){

  params <- read_rds("data/posteriorSamples.rds")
  params <- params |> slice(sample.int(nrow(params), 1))

  draw_value <- function(x){
    params |>
      select(contains(x)) |>
      pivot_longer(cols = everything()) |>
      pull(value) |>
      jitter()

  }

  constants_nimble$phi_mu <- draw_value("phi_mu")
  constants_nimble$psi_phi <- draw_value("psi_phi")
  log_nu <- draw_value("log_nu")
  constants_nimble$zeta <- 28 * exp(log_nu) / 365
  constants_nimble$beta_p <- matrix(draw_value("beta_p"), 5, 3)
  constants_nimble$beta1 <- draw_value("beta1")
  constants_nimble$omega <- draw_value("p_mu")
  constants_nimble$log_gamma <- draw_value("gamma")
  constants_nimble$log_rho <- draw_value("rho")

  # for testing
  # attach(constants_nimble)
  # attach(data_nimble)

  with(append(constants_nimble, data_nimble), {

    a <- phi_mu * psi_phi
    b <- (1 - phi_mu) * psi_phi
    M <- phi <- lambda <- rep(NA, max(mH, na.rm = TRUE))
    M_init <- rep(NA, n_cluster)
    # i=1
    for(i in 1:n_cluster){

      M_init[i] <- round(exp(log_cluster_area[i])*start_density + sum(rem[i, ], na.rm = TRUE))
      M[mH[i, 1]] <- M_init[i]
      # j=2
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
      p_mu = omega,
      p_unique =  boot::inv.logit(omega),
      phi_mu = phi_mu,
      psi_phi = psi_phi,
      a_phi = a,
      b_phi = b,
      N = N + buffer,
      M = M + buffer,
      # lambda = lambda + buffer,
      log_nu = log_nu,
      log_gamma = log_gamma,
      log_rho = log_rho,
      phi = phi,
      zeta = zeta,
      log_zeta = log(zeta)
    )

  })
}






