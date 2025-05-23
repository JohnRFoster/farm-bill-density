
simulate_cluster_dynamics <- function(start_density, prop_ls, n_pp, include_project, include_cluster){

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

  # habitat_data <- read_csv("data/habitat_proportions.csv")

  hierarchy <- "Poisson"

  params <- read_rds("data/posteriorSamples.rds")
  params <- params |> slice(sample.int(nrow(params), 1))

  draw_value <- function(x){
    params |>
      select(contains(x)) |>
      pivot_longer(cols = everything()) |>
      pull(value)

  }

  phi_mu <- draw_value("phi_mu")
  psi_phi <- draw_value("psi_phi")
  log_nu <- draw_value("log_nu")
  nu <- exp(log_nu)
  # beta_p <- matrix(draw_value("beta_p"), 5, 3, byrow = TRUE)
  beta1 <- matrix(draw_value("beta1"), 5, 1)
  beta_p <- matrix(
    runif(20, -2, 2),
    nrow = 5,
    ncol = 4
  )
  beta <- cbind(beta1, beta_p)
  p_unique <- draw_value("p_mu")
  omega <- boot::inv.logit(p_unique)
  log_gamma <- draw_value("log_gamma")
  gamma <- exp(log_gamma)
  log_rho <- draw_value("log_rho")
  rho <- exp(log_rho)

  a_phi <- phi_mu * psi_phi
  b_phi <- (1 - phi_mu) * psi_phi
  zeta <- 28 * nu / 365

  # there can be up to m0 primary periods after all pigs are gone before we drop data
  m0 <- 6

  get_attribute <- function(ls, name){
    1:length(ls) |>
      map(
        \(x) pluck(ls[[x]], name)
      ) |>
      list_c()
  }

  group_lookup <- tibble(
    property = 1:length(prop_ls),
    projects = get_attribute(prop_ls, "project"),
    clusters = get_attribute(prop_ls, "cluster"),
    property_area = get_attribute(prop_ls, "property_area"),
    cluster_area = get_attribute(prop_ls, "cluster_area")
  ) |>
    filter(cluster_area >= 1.8)

  n_projects <- length(unique(group_lookup$projects))
  n_clusters <- length(unique(group_lookup$clusters))

  # project and cluster random effects if used
  if(include_project){
    tau_project <- runif(1, 0.001, 4)
    alpha_project <- rnorm(n_projects, 0, prec_2_sd(tau_project))
  } else {
    tau_project <- NA
    alpha_project <- numeric(n_projects)
  }

  project_lookup <- group_lookup |>
    select(projects) |>
    distinct() |>
    mutate(tau_project = tau_project,
           alpha_project = alpha_project)

  if(include_cluster){
    tau_cluster <- runif(1, 0.001, 4)
    alpha_cluster <- rnorm(n_clusters, 0, prec_2_sd(tau_cluster))
  } else {
    tau_cluster <- NA
    alpha_cluster <- numeric(n_clusters)
  }

  cluster_lookup <- group_lookup |>
    select(clusters) |>
    distinct() |>
    mutate(tau_cluster = tau_cluster,
           alpha_cluster = alpha_cluster)

  group_lookup <- group_lookup |>
    left_join(project_lookup) |>
    left_join(cluster_lookup)

  all_take <- all_pigs <- known_abundance <- tibble()
  cluster_vec <- unique(group_lookup$clusters)

  pb <- txtProgressBar(max = n_clusters, style = 1)
  for(i in seq_len(n_clusters)){

    prop_tmp_in_cluster <- group_lookup |>
      filter(clusters %in% cluster_vec[i])

    area <- prop_tmp_in_cluster$cluster_area[1]
    prop_tmp <- prop_tmp_in_cluster$property
    survey_area <- prop_tmp_in_cluster$property_area
    log_survey_area <- log(survey_area)
    project <- prop_tmp_in_cluster$projects
    project_re <- prop_tmp_in_cluster$alpha_project
    cluster_re <- prop_tmp_in_cluster$alpha_cluster

    n_props <- length(prop_tmp)

    X <- matrix(
      data = rnorm(n_props * 4, 0, 0.5),
      nrow = n_props,
      ncol = 4
    )

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

    tka <- tibble(
      project = project,
      cluster = cluster_vec[i],
      cluster_area_km2 = area,
      property = prop_tmp,
      property_area_km2 = survey_area,
      M = M[1],
      N = N[,1],
      PPNum = 1
    )

    known_abundance <- bind_rows(known_abundance, tka)

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
            p_unique = omega,
            method_lookup = method_lookup,
            alpha_project = project_re[j],
            alpha_cluster = cluster_re[j]) |>
            mutate(property = prop_tmp[j],
                   property_area_km2 = survey_area[j],
                   cluster_area_km2 = area,
                   cluster = cluster_vec[i],
                   project = project[j],
                   N = nn,
                   M = M[t],
                   M_actual = M_actual[t],
                   agriculture = X[j, 1],
                   forested = X[j, 2],
                   grassland = X[j ,3],
                   wetland = X[j, 4])

          take <- bind_rows(take, take_t)
          Y[j] <- sum(take_t$take)

          # how many pigs are left?
          zi <- nn - Y[j]

          identifier <- paste0("\nCluster ", i, "\nProperty ", j, "\nTime ", t)
          if(is.na(Y[j])) print(identifier)

          assertthat::assert_that(zi >= 0,
                                  msg = paste("More pigs removed than are alive!", identifier))

        } # remove pigs loop

      } # property loop

      # need to use the realized cluster size in case all pigs are removed
      Z <- max(c(M[t], M_actual[t])) - sum(Y)

      identifier <- paste0("\nCluster ", i, "\nTime ", t)
      assertthat::assert_that(Z >= 0,
                              msg = paste("More pigs removed from cluster than are alive!", identifier))

      if(t < n_pp){
        M[t+1] <- process_model(Z, zeta, a_phi, b_phi)

        if(single_property_cluster){
          N[,t+1] <- M[t+1]
          M_actual[t+1] <- M[t+1]
        } else {

          # distribute based on density
          d <- M[t+1] / area * survey_area

          N[,t+1] <- rpois(n_props, d)
          M_actual[t+1] <- sum(N[,t+1])
        }

        tka <- tibble(
          project = project,
          cluster = cluster_vec[i],
          cluster_area_km2 = area,
          property = prop_tmp,
          property_area_km2 = survey_area,
          M = M[t+1],
          N = N[,t+1],
          PPNum = t+1
        )

        known_abundance <- bind_rows(known_abundance, tka)

         # alpha <- log(d)
        # lambda <- rnorm(n_props, alpha, sigma)
        # N[,t] <- round(exp(lambda))

        e_count <- if_else(M[t+1] == 0, e_count + 1, e_count)
        if(e_count >= m0) break

      }


    } # time loop

    all_take <- bind_rows(all_take, take)

    setTxtProgressBar(pb, i)

  } # cluster loop
  close(pb)

  take_filter <- function(df, x){
    df |>
      group_by(.data[[x]]) |>
      summarise(total_take = sum(take)) |>
      ungroup() |>
      filter(total_take > 0) |>
      pull(.data[[x]])
  }

  project_filter <- all_take |>
    take_filter("project")

  cluster_filter <- all_take |>
    take_filter("cluster")

  good_projects1 <- all_take |>
    filter(project %in% project_filter) |>
    select(project, cluster) |>
    distinct() |>
    count(project) |>
    filter(n > 1) |>
    pull(project)

  good_projects2 <- all_take |>
    filter(project %in% good_projects1) |>
    select(project, PPNum) |>
    distinct() |>
    count(project) |>
    filter(n >= 2) |>
    pull(project)

  good_clusters <- all_take |>
    filter(project %in% good_projects2,
           cluster %in% cluster_filter) |>
    select(cluster, PPNum) |>
    distinct() |>
    count(cluster) |>
    filter(n >= 2) |>
    pull(cluster)

  all_take <- left_join(all_take, group_lookup) |>
    select(-projects, -clusters) |>
    filter(cluster %in% good_clusters,
           project %in% good_projects2) |>
    arrange(property, PPNum, order)

  all_pigs <- known_abundance |>
    filter(!is.na(M),
           cluster %in% good_clusters,
           project %in% good_projects2,
           property %in% unique(all_take$property)) |>
    arrange(property, PPNum)

  take_check <- all_take |>
    select(project, cluster, property) |>
    distinct()

  pigs_check <- all_pigs |>
    select(project, cluster, property) |>
    distinct()

  assertthat::assert_that(nrow(take_check) == nrow(pigs_check))

  assertthat::assert_that(all(take_check$project == pigs_check$project),
                          msg = "projects do not align")

  assertthat::assert_that(all(take_check$cluster == pigs_check$cluster),
                          msg = "clusters do not align")

  assertthat::assert_that(all(take_check$property == pigs_check$property),
                          msg = "properties do not align")

  group_fac <- take_check |>
    mutate(cluster_fac = as.numeric(as.factor(cluster)),
           project_fac = as.numeric(as.factor(project)),
           property_fac = as.numeric(as.factor(property)))

  take_return <- left_join(all_take, group_fac) |>
    select(-cluster, -project, -property) |>
    rename(cluster = cluster_fac,
           project = project_fac,
           property = property_fac)

  pigs_return <- left_join(all_pigs, group_fac) |>
    select(-cluster, -project, -property) |>
    rename(cluster = cluster_fac,
           project = project_fac,
           property = property_fac)


  return(list(
    all_pigs = pigs_return,
    all_take = take_return,
    beta = beta
  ))
}


