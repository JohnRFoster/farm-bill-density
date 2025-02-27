
simulate_cluster_dynamics <- function(start_density, cluster_props, properties){

  library(dplyr)
  library(tidyr)
  library(readr)

  source("R/functions_removal.R")
  source("R/functions_misc.R")

  effort_data <- read_csv("data/effort_data.csv") |>
    mutate(method_name = if_else(method_name == "TRAPS, CAGE", "TRAPS", method_name)) |>
    suppressMessages()

  method_lookup <- effort_data |>
    select(method, method_name) |>
    mutate(idx = method) |>
    distinct()

  start_density <- start_density
  hierarchy <- "Poisson"

  params <- read_csv("data/posterior_95CI_range_all.csv") |>
    suppressMessages()

  draw_value <- function(x){
    params |>
      filter(grepl(x, node)) |>
      pull(med) |>
      round(2)
  }

  phi_mu <- draw_value("phi_mu")
  psi_phi <- draw_value("psi_phi")
  nu <- draw_value("nu")
  beta_p <- matrix(draw_value("beta_p"), 5, 3, byrow = TRUE)
  beta1 <- matrix(draw_value("beta1"), 5, 1)
  beta <- cbind(beta1, beta_p)
  omega <- draw_value("omega")
  gamma <- draw_value("gamma")
  log_gamma <- log(gamma)
  rho <- draw_value("rho")
  log_rho <- log(rho)

  a_phi <- phi_mu * psi_phi
  b_phi <- (1 - phi_mu) * psi_phi
  zeta <- 28 * nu / 365

  # there can be up to m0 primary periods after all pigs are gone before we drop data
  m0 <- 6

  n_clusters <- max(cluster_props$cluster)

  p_unique <- omega

  cluster_area <- cluster_props |>
    group_by(cluster) |>
    summarise(cluster_area = sum(property_area_km2)) |>
    ungroup() |>
    pull(cluster_area)

  all_take <- all_pigs <- tibble()

  pb <- txtProgressBar(max = n_clusters, style = 1)
  for(i in seq_len(n_clusters)){

    properties_in_cluster <- cluster_props |>
      filter(cluster == i)

    if(properties_in_cluster$drop_flag[1] == 1) next

    area <- cluster_area[i]

    propertyIDs <- properties_in_cluster |>
      pull(propertyID)

    survey_area <- properties_in_cluster |>
      pull(property_area_km2)

    log_survey_area <- log(survey_area)

    n_props <- length(survey_area)
    X <- matrix(rnorm(3 * n_props), n_props, 3)

    single_property_cluster <- if_else(n_props == 1, TRUE, FALSE)

    Mspin <- rpois(1, area * start_density)

    for(t in 1:6){
      Mspin <- process_model(Mspin, zeta, a_phi, b_phi)
    }

    # skip clusters that go extinct
    if(Mspin == 0) next

    # distribute based on density
    d <- Mspin / area * survey_area

    N <- matrix(NA, n_props, n_pp)
    M <- M_actual <- rep(NA, n_pp)
    M[1] <- Mspin

    if(single_property_cluster){
      N[,1] <- d
      M_actual[1] <- d
    } else {
      N[,1] <- rpois(n_props, d)
      M_actual[1] <- sum(N[,1])
    }

    # if we choose to use a normal HB model
    # alpha <- log(d)
    # lambda <- rnorm(n_props, alpha, sigma)
    # N[,1] <- round(exp(lambda))

    property_data <- properties[propertyIDs]

    e_count <- 0
    take <- tibble()
    for(t in 1:(n_pp-1)){

      # reset catch
      Y <- numeric(n_props)

      for(j in seq_len(n_props)){
        tmp_data <- property_data[[j]]

        removal_effort <- tmp_data$effort
        sample_occasions <- removal_effort$sample_occasions

        remove_pigs <- t %in% sample_occasions

        if(remove_pigs){
          sample_effort <- removal_effort |> filter(sample_occasions == t)

          # determine order of removal events
          removal_order <- determine_removal_order(sample_effort)

          # conduct removals
          nn <- N[j, t]

          take_t <- conduct_removals(nn, removal_order, effort_data, log_survey_area[j],
                                     X[j,], beta, t,
                                     log_rho, log_gamma, p_unique, method_lookup) |>
            mutate(propertyID = propertyIDs[j],
                   property_area_km2 = survey_area[j],
                   cluster_area_km2 = area,
                   cluster = i,
                   N = nn,
                   M = M[t],
                   M_actual = M_actual[t],
                   c_road_den = X[j, 1],
                   c_rugged = X[j, 2],
                   c_canopy = X[j ,3])

          take <- bind_rows(take, take_t)
          Y[j] <- sum(take_t$take)

          # how many pigs are left?
          zi <- nn - Y[j]

          identifier <- paste0("\nCluster ", i, "\nProperty ", j, "\nTime ", t)
          if(is.na(Y[j])) print(identifier)

          assertthat::assert_that(zi >= 0,
                                  msg = paste("More pigs removed than are alive!", identifier))

        }

      }

      # need to use the realized cluster size in case all pigs are removed
      Z <- M_actual[t] - sum(Y)

      identifier <- paste0("\nCluster ", i, "\nTime ", t)
      assertthat::assert_that(Z >= 0,
                              msg = paste("More pigs removed from cluster than are alive!", identifier))

      if(t < n_pp){
        M[t+1] <- process_model(Z, zeta, a_phi, b_phi)

        # distribute based on density
        d <- M[t+1] / area * survey_area

        if(single_property_cluster){
          N[,t+1] <- d
          M_actual[t+1] <- d
        } else {
          N[,t+1] <- rpois(n_props, d)
          M_actual[t+1] <- sum(N[,t+1])
        }

        # alpha <- log(d)
        # lambda <- rnorm(n_props, alpha, sigma)
        # N[,t] <- round(exp(lambda))

        e_count <- if_else(M_actual[t+1] == 0, e_count + 1, e_count)
        if(e_count >= m0) break

      }


    }

    M_df <- tibble(
      M = M,
      M_actual = M_actual,
      PPNum = 1:n_pp,
      cluster = i
    )

    colnames(N) <- 1:n_pp
    N_df <- N |>
      as_tibble() |>
      mutate(propertyID = propertyIDs) |>
      pivot_longer(cols = -propertyID,
                   names_to = "PPNum",
                   values_to = "N") |>
      mutate(PPNum = as.numeric(PPNum))

    pigs_i <- left_join(M_df, N_df) |> suppressMessages()

    all_pigs <- bind_rows(all_pigs, pigs_i)
    all_take <- bind_rows(all_take, take)

    setTxtProgressBar(pb, i)

  }
  close(pb)

  good_clusters <- all_take |>
    select(cluster, PPNum) |>
    distinct() |>
    group_by(cluster) |>
    mutate(nt = 1:n()) |>
    filter(max(nt) > 1) |>
    pull(cluster) |>
    unique()

  all_take <- all_take |>
    filter(cluster %in% good_clusters) |>
    mutate(cluster = as.numeric(as.factor(cluster)),
           property = as.numeric(as.factor(propertyID))) |>
    arrange(property, PPNum, order)

  all_pigs <- all_pigs |>
    filter(!is.na(M),
           cluster %in% good_clusters,
           propertyID %in% unique(all_take$propertyID)) |>
    mutate(cluster = as.numeric(as.factor(cluster)),
           property = as.numeric(as.factor(propertyID))) |>
    arrange(property, PPNum)


  return(list(
    all_pigs = all_pigs,
    all_take = all_take
  ))
}


