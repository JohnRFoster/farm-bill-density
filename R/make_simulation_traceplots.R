# --------------------------------------------------------------------
#
# create traceplots for simulations
#
# --------------------------------------------------------------------

library(dplyr)
library(tidyr)
library(readr)
library(ggplot2)
library(coda)

source("R/functions_misc.R")

# function to create traceplots from thinned posterior
trace_plot <- function(post, nodes_2_plot, thin = 5000){
  df <- post |>
    select(chain, all_of(nodes_2_plot)) |>
    group_by(chain) |>
    mutate(iteration = 1:n()) |>
    ungroup()

  posterior_mat <- df |>
    select(-chain, -iteration) |>
    as.matrix()

  print(apply(posterior_mat, 2, quantile, c(0.025, 0.5, 0.975)))

  total_iterations <- max(df$iteration)
  thin_interval <- floor(seq(1, total_iterations, length.out = thin))

  gg <- df |>
    filter(iteration %in% thin_interval) |>
    pivot_longer(cols = -c(iteration, chain),
                 names_to = "node") |>
    mutate(chain = as.character(chain)) |>
    ggplot() +
    aes(x = iteration, y = value, color = chain) +
    geom_line() +
    facet_wrap(~ node, scales = "free_y") +
    labs(x = "Iteration",
         y = "Value") +
    theme_bw()
  return(gg)
}

config_name <- "hpc_test"
config <- config::get(config = config_name)

include_project <- config$include_project
include_cluster <- config$include_cluster

message("Project random effect: ", include_project)
message("Cluster random effect: ", include_cluster)

if(include_project){
  model_dir <- if_else(include_cluster, "project_cluster", "project")
} else {
  model_dir <- if_else(include_cluster, "cluster", "base")
}


# simulations are stored here
out_dir <- file.path(config$project_dir, config$out_dir, config$dev_dir, model_dir)

sim_files <- list.files(out_dir)


for(i in seq_along(sim_files)){

  task_id <- sim_files[i]

  task_dir <- file.path(out_dir, task_id)
  message("\n\n====", task_dir, "====")

  out_file <- file.path(task_dir, "parameterSamples.rds")

  # traceplots_exist <- any(grepl("mcmcTimeseries_", list.files(task_dir)))
  # if(traceplots_exist) next

  if(file.exists(out_file)){
    params_mcmc_list <- read_rds(out_file)[[1]]
    n_chains <- length(params_mcmc_list)

    posterior <- tibble()
    message("Create posterior tibble...")
    pb <- txtProgressBar(min = 1, max = n_chains, style = 1)
    for(c in seq_len(n_chains)){
      chain_c <- as.matrix(params_mcmc_list[[c]]) |>
        as_tibble() |>
        mutate(chain = c)
      posterior <- bind_rows(posterior, chain_c)

      setTxtProgressBar(pb, c)
    }
    close(pb)
    message("  done")


  } else {
    next
  }

  nodes <- setdiff(colnames(posterior), "chain")
  n_plots_per_page <- 4   # want to put 4 nodes on a single plot
  idx <- rep(seq(1, ceiling(length(nodes) / n_plots_per_page), by = 1),
             each = n_plots_per_page)[1:length(nodes)]

  plots <- tibble(
    nodes = nodes,
    idx = idx
  )

  message("Creating traceplots...")
  # options(bitmapType = 'Xlib')
  pb <- txtProgressBar(min = 1, max = max(plots$idx), style = 1)
  for(j in seq_along(unique(plots$idx))){
    n2p <- plots |>
      filter(idx == j) |>
      pull(nodes)

    gg <- trace_plot(posterior, n2p)

    filename <- file.path(task_dir, paste0("mcmcTimeseries_", sprintf("%03d", j), ".pdf"))
    ggsave(filename, gg)

    setTxtProgressBar(pb, j)
  }
  close(pb)

  message("Parameters Done")

  params_samples <- posterior |>
    as_tibble() |>
    mutate(
      `gamma[1]` = exp(`log_gamma[1]`),
      `gamma[2]` = exp(`log_gamma[2]`),
      `omega[1]` = boot::inv.logit(`p_mu[1]`),
      `omega[2]` = boot::inv.logit(`p_mu[2]`),
      `rho[1]` = exp(`log_rho[1]`),
      `rho[2]` = exp(`log_rho[2]`),
      `rho[3]` = exp(`log_rho[3]`),
      `rho[4]` = exp(`log_rho[4]`),
      `rho[5]` = exp(`log_rho[5]`),
      nu = exp(log_nu),
    )

  data <- read_rds(file.path(task_dir, "simulationData.rds"))

  if(model_dir == "base"){
    params_actual <- data$params_actual |>
      filter(
        !grepl("alpha", node),
        !grepl("tau", node)
      )
  } else if (model_dir == "project"){
    params_actual <- data$params_actual |>
      filter(
        !grepl("alpha_cluster", node),
        !grepl("tau_cluster", node)
      )
  } else if (model_dir == "cluster"){
    params_actual <- data$params_actual |>
      filter(
        !grepl("alpha_project", node),
        !grepl("tau_project", node)
      )
  } else if (model_dir == "project_cluster"){
    params_actual <- data$params_actual
  }

  make_prior <- function(p, sd = NULL){
    n <- 10000

    if(grepl("beta", p)) x <- rnorm(n, 0, 1)
    if(grepl("rho", p)) x <- exp(rnorm(n, 0, 1))
    if(grepl("gamma", p)) x <- exp(rnorm(n, 0, prec_2_sd(3)))
    if(grepl("nu", p)) x <- exp(rnorm(n, 2, 1))
    if(grepl("omega", p)) x <- boot::inv.logit(rnorm(n, 0, 1))
    if(grepl("tau", p)) x <- rgamma(n, 1, 1)
    if(grepl("phi_mu", p)) x <- rbeta(n, 1, 1)
    if(grepl("psi", p)) x <- rgamma(n, 1, 0.1)
    if(grepl("alpha", p)) x <- rnorm(length(sd), 0, sd)

    tibble(
      value = x,
      dist = "Prior"
    )

  }

  g <- list()
  for(j in 1:nrow(params_actual)){

    p <- params_actual$node[j]
    v <- params_actual$actual[j]

    post <- tibble(
      dist = "Posterior",
      value = params_samples[[p]]
    )

    if(grepl("alpha", p)){
      if(grepl("cluster", p)) sd <- prec_2_sd(params_samples[["tau_cluster"]])
      if(grepl("project", p)) sd <- prec_2_sd(params_samples[["tau_project"]])
    } else {
      sd <- NULL
    }

    prior <- make_prior(p, sd)

    df <- bind_rows(post, prior) |>
      mutate(dist = factor(dist, levels = c("Prior", "Posterior")))

    xlim <- df |>
      filter(dist == "Posterior") |>
      pull(value) |>
      range()

    g[[j]] <- df |>
      ggplot() +
      aes(x = value,
          y = after_stat(scaled),
          fill = dist,
          color = dist) +
      geom_density(alpha = 0.5) +
      geom_vline(xintercept = v, linewidth = 1, linetype = "dashed") +
      scale_color_manual(values = c(
        "Posterior" = "black",
        "Prior" = "red"
      )) +
      scale_fill_manual(values = c(
        "Posterior" = "lightblue",
        "Prior" = "white"
      )) +
      coord_cartesian(xlim = xlim) +
      labs(x = "Value",
           y = "Density",
           color = "Distribution",
           fill = "Distribution",
           title = p) +
      theme_bw()

  }

  ggg <- do.call(ggpubr::ggarrange,
                 list(
                   plotlist = g,
                   nrow = 4,
                   ncol = 4,
                   common_legend = TRUE,
                   legend = "none"
                 )
  )

  for(j in 1:length(ggg)){
    filename <- file.path(task_dir, paste0("priorPosterior_", sprintf("%03d", j), ".pdf"))
    ggsave(filename, ggg[[j]], units = "cm", width = 18, height = 18)
  }




}

