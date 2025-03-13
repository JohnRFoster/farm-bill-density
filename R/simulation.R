
simulate_cluster_dynamics <- function(start_density, cluster_props, prop_ls){

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
  p_unique <- omega
  gamma <- draw_value("gamma")
  log_gamma <- log(gamma)
  rho <- draw_value("rho")
  log_rho <- log(rho)

  a_phi <- phi_mu * psi_phi
  b_phi <- (1 - phi_mu) * psi_phi
  zeta <- 28 * nu / 365

  # there can be up to m0 primary periods after all pigs are gone before we drop data
  m0 <- 6

  projects <- unique(cluster_props$project)
  clusters <- unique(cluster_props$cluster)

  n_projects <- length(projects)
  n_clusters <- length(clusters)

  # project and cluster random effects
  tau_project <- runif(1, 1e-4, 100)
  tau_cluster <- runif(1, 1e-4, 100)

  project_lookup <- tibble(
    project = projects,
    tau_project = tau_project,
    alpha_project = rnorm(n_projects, 0, 1/sqrt(tau_project))
  )

  cluster_lookup <- tibble(
    cluster = clusters,
    tau_cluster = tau_cluster,
    alpha_cluster = rnorm(n_clusters, 0, 1/sqrt(tau_cluster))
  )

  group_lookup <- cluster_props |>
    select(property, project, cluster) |>
    distinct() |>
    left_join(project_lookup) |>
    left_join(cluster_lookup)

  all_take <- all_pigs <- tibble()

  pb <- txtProgressBar(max = n_clusters, style = 1)
  for(i in seq_len(n_clusters)){

    prop_tmp_in_cluster <- cluster_props |>
      filter(cluster == clusters[i])

    if(prop_tmp_in_cluster$drop_flag[1] == 1) next

    area <- prop_tmp_in_cluster$cluster_area_km2[1]
    prop_tmp <- prop_tmp_in_cluster$property
    survey_area <- prop_tmp_in_cluster$property_area_km2
    project <- prop_tmp_in_cluster$project

    cp_group <- group_lookup |> filter(property %in% prop_tmp)

    project_re <- cp_group$alpha_project
    cluster_re <- cp_group$alpha_cluster

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

    N <- matrix(NA, n_props, n_pp)
    M <- M_actual <- rep(NA, n_pp)
    M[1] <- Mspin

    if(single_property_cluster){
      N[,1] <- Mspin
      M_actual[1] <- Mspin
    } else {

      # distribute based on density
      d <- Mspin / area * survey_area

      N[,1] <- rpois(n_props, d)
      M_actual[1] <- sum(N[,1])
    }

    # if we choose to use a normal HB model
    # alpha <- log(d)
    # lambda <- rnorm(n_props, alpha, sigma)
    # N[,1] <- round(exp(lambda))

    property_data <- prop_ls[prop_tmp]

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

          take_t <- conduct_removals(
            N = nn,
            removal_order = removal_order,
            effort_data = effort_data,
            log_survey_area = log_survey_area[j],
            X = X[j,],
            beta_p = beta,
            pp = t,
            log_rho = log_rho,
            log_gamma = log_gamma,
            p_unique = p_unique,
            method_lookup = method_lookup,
            alpha_project = project_re[j],
            alpha_cluster = cluster_re[j]) |>
            mutate(property = prop_tmp[j],
                   property_area_km2 = survey_area[j],
                   cluster_area_km2 = area,
                   cluster = i,
                   project = project[j],
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

        if(single_property_cluster){
          N[,t+1] <- Mspin
          M_actual[t+1] <- Mspin
        } else {

          # distribute based on density
          d <- Mspin / area * survey_area

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
      mutate(property = prop_tmp) |>
      pivot_longer(cols = -property,
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

  all_take <- left_join(all_take, group_lookup) |>
    filter(cluster %in% good_clusters) |>
    mutate(cluster = as.numeric(as.factor(cluster)),
           project = as.numeric(as.factor(project)),
           property = as.numeric(as.factor(property))) |>
    arrange(property, PPNum, order)

  all_pigs <- all_pigs |>
    filter(!is.na(M),
           cluster %in% good_clusters,
           property %in% unique(all_take$property)) |>
    mutate(cluster = as.numeric(as.factor(cluster)),
           project = as.numeric(as.factor(project)),
           property = as.numeric(as.factor(property))) |>
    arrange(property, PPNum)


  return(list(
    all_pigs = all_pigs,
    all_take = all_take
  ))
}


