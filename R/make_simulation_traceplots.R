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

}

