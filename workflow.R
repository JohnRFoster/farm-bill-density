# --------------------------------------------------------------------
#
# Workflow script for invasive species removal simulations
# for Farm Bill properties that have lat lon data
#
# Author:
# John Foster
#
# General workflow (conducted with run_simulation):
#
# 1. Load Farm Bill data
#
# 2. run_simulation
#   - Simulate/bootstrap data
#     - 1-method properties
#     - 2- to 4-method properties
#   - Simulate eco/take dynamics
#   - Fit MCMC
#   - Check MCMC
# 3. Summarize output
#
# --------------------------------------------------------------------

# config_name <- "default"
config_name <- "hpc_test"
# config_name <- "hpc_production"
config <- config::get(config = config_name)
interval <- 4

if(grepl("hpc", config_name)){
  Sys.setenv(RENV_CONFIG_SANDBOX_ENABLED = FALSE)
  renv::load("/home/john.foster/farm-bill-density/")

  args <- commandArgs(trailingOnly = TRUE)
  task_id <- args[1]

} else {
  task_id <- 1
}

message("Task ID: ", task_id)
set.seed(task_id)

include_project <- config$include_project
include_cluster <- config$include_cluster

message("Project random effect: ", include_project)
message("Cluster random effect: ", include_cluster)

library(nimble)
library(parallel)
library(coda)
library(readr)
library(dplyr)
library(tidyr)
library(purrr)
library(lubridate)

source("R/functions_misc.R")
source("R/make_properites.R")
source("R/functions_prep_nimble_simulation.R")
source("R/simulation.R")
source("R/fit_mcmc_simulation.R")
source("R/nimble_cluster_model.R")
source("R/check_mcmc.R")
source("R/n_method_properties.R")
source("R/one_method_properties.R")

# -----------------------------------------------------------------
# Load MIS data ----
# -----------------------------------------------------------------
message("MIS data intake")
df <- read_csv("data/clusters250km2.csv") |>
  filter(STATE == "OK") |>
  mutate(property = as.numeric(as.factor(propertyID)))

end_dates <- unique(sort(df$end.date))
min_date <- min(end_dates)
max_date <- max(end_dates)

start_dates <- seq(min_date, max_date, by = paste(interval, "week"))
end_dates <- c(start_dates[-1] - 1, max_date)

targets::tar_assert_identical(length(start_dates), length(end_dates))
# targets::tar_assert_true(min(df$start.date) >= min_date)
# targets::tar_assert_true(max(df$start.date) <= max_date)

timestep_df <- tibble(start_dates, end_dates) |>
  mutate(timestep = 1:n())
timestep_df$month <- month(timestep_df$end_dates)
timestep_df$year <- year(timestep_df$end_dates)

# for each row in the merged data, insert the integer primary period timestep
df$timestep <- NA
after_start <- before_end <- rep(NA, nrow(df))
message("Assign timesteps...")
pb <- txtProgressBar(max = nrow(df), style = 1)
for (i in 1:nrow(df)) {
  after_start[i] <- which(timestep_df$start_dates <= df$start.date[i]) |> max()
  before_end[i] <- which(timestep_df$end_dates >= df$end.date[i]) |> min()
  if (after_start[i] == before_end[i]) {
    # then the start and end date is contained within a primary period
    df$timestep[i] <- timestep_df$timestep[before_end[i]]
  } # otherwise, timestep[i] will be left as NA and filtered out later
  setTxtProgressBar(pb, i)
}
close(pb)

write_csv(timestep_df, file.path("data/timestep_df.csv"))

df_with_timesteps <- df |>
  filter(!is.na(timestep)) |>
  left_join(timestep_df) |>
  filter(year >= StartYr,
         year <= LastYr) |>
  arrange(property, timestep)

# need a lookup table for the number of properties in each cluster ----
n_property_lookup <- df_with_timesteps |>
  select(property, cluster) |>
  distinct() |>
  count(cluster, name = "n_props_in_cluster")

# get the proportion of clusters that have n properties
# for determining the number of clusters to simulate
n_rel <- n_property_lookup |>
  count(n_props_in_cluster, name = "n_clusters") |>
  mutate(rel_prop = n_clusters / sum(n_clusters),
         n_clusters_to_sim = ceiling(rel_prop * 40),
         n_props_to_sim = n_clusters_to_sim * n_props_in_cluster)

n_clusters_to_sim <- sum(n_rel$n_clusters_to_sim)
n_properties_to_sim <- sum(n_rel$n_props_to_sim)
n_props_in_cluster <- n_rel$n_props_in_cluster

start_density <- round(runif(1, 0.3, 5), 3)
message("Starting density: ", start_density)
message("n clusters: ", n_clusters_to_sim)
message("n properties: ", n_properties_to_sim)

cluster_breaks <- c(0, cumsum(n_rel$n_clusters_to_sim))
property_breaks <- c(0, cumsum(n_rel$n_props_to_sim))

all_clusters <- tibble()
for(i in 1:(length(cluster_breaks) - 1)){

  clust <- tibble(
    cluster = rep((cluster_breaks[i] + 1):(cluster_breaks[i + 1]),
                  each = n_props_in_cluster[i]),
    property = (property_breaks[i] + 1):(property_breaks[i + 1]),
  )

  all_clusters <- bind_rows(all_clusters, clust)

}

n_method <- 0
while(n_method != 5){

  properties <- list()

  for(i in 1:nrow(n_rel)){

    properties_i <- make_properites(
      df = df_with_timesteps,
      n_clusters = n_rel$n_clusters_to_sim[i],
      n_props_in_each_cluster = n_props_in_cluster[i],
      n_pp = config$n_pp)

    cluster_vec <- rep((cluster_breaks[i] + 1):(cluster_breaks[i + 1]),
                       each = n_props_in_cluster[i])

    properties_i <- properties_i |> map2(1:length(properties_i),
                                         \(x, y) assign_in(x, "cluster", cluster_vec[y]))

    properties <- c(properties, properties_i)

  }

  message("\nSimulate cluster dynamics")
  simulated_data <- simulate_cluster_dynamics(
    start_density = start_density,
    prop_ls = properties,
    n_pp = config$n_pp,
    include_project,
    include_cluster
  )

  take <- simulated_data$all_take
  abundance <-  simulated_data$all_pigs

  n_method <- length(unique(take$method))

}

if(config_name == "default"){
  library(ggplot2)

  take |>
    mutate(cluster = as.character(cluster),
           M = M / cluster_area_km2) |>
    ggplot() +
    aes(x = PPNum, y = M, color = cluster) +
    labs(y = "density") +
    geom_line()

}

take |>
  select(property, cluster) |>
  distinct() |>
  count(cluster, name = "n_properties_in_cluster") |>
  count(n_properties_in_cluster)

take |>
  select(project, cluster) |>
  distinct() |>
  arrange(project, cluster)

message("Prep data for MCMC")
nimble_ls <- prep_nimble(take, config$posterior_path) |> suppressMessages()
nimble_data <- nimble_ls$data
nimble_constants <- nimble_ls$constants

monitors_add <- c("N", "M")

params_check <- c(
  "alpha_cluster",
  "alpha_project",
  "beta1",
  "beta_p",
  "log_gamma",
  "log_nu",
  "log_rho",
  "p_mu",
  "phi_mu",
  "psi_phi"
)

if(include_project){
  monitors_add <- c(monitors_add, "alpha_project")
  params_check <- c(params_check, "tau_project")
}

if(include_cluster){
  monitors_add <- c(monitors_add, "alpha_cluster")
  params_check <- c(params_check, "tau_cluster")
}

n_chains <- config$n_chains
cl <- makeCluster(n_chains)

c_samp <- tibble(
  node = c("phi_mu", "psi_phi"),
  type = c("slice",  "slice")
)

samples <- fit_mcmc(
  cl = cl,
  task_seed = task_id,
  modelCode = modelCode,
  data = nimble_ls$data,
  constants = nimble_ls$constants,
  start_density = start_density,
  n_iter = config$n_iter,
  n_chains = n_chains,
  custom_samplers = c_samp,
  monitors_add = monitors_add,
  include_project = include_project,
  include_cluster = include_cluster
)

stopCluster(cl)

if(include_project){
  model_dir <- if_else(include_cluster, "project_cluster", "project")
} else {
  model_dir <- if_else(include_cluster, "cluster", "base")
}

out_dir <- file.path(config$project_dir, config$out_dir, config$dev_dir, model_dir, task_id)
message("\n\nWriting to: ", out_dir)

if(!dir.exists(out_dir)) dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

samples_mcmc <- as.mcmc.list(samples)

out <- check_mcmc(
  samples = samples_mcmc,
  nodes_check = params_check,
  dest = out_dir
)

all_time_ids <- make_all_pp(take)
known_abundance <- left_join(all_time_ids, abundance) |>
  # filter(!is.na(M)) |>
  mutate(propertySimID = paste(task_id, property, start_density, sep = "-"),
         clusterSimID = paste(task_id, cluster, start_density, sep = "-"),
         projectSimID = paste(task_id, project, start_density, sep = "-"))

add_ids <- function(df){
  df |>
    mutate(task_id = task_id,
           start_density = start_density)
}

out$known_abundance <- known_abundance |> add_ids()
out$data <- nimble_ls$data
out$constants <- nimble_ls$constants
out$start_density <- start_density
out$task_id <- task_id
out$take <- take |>
  mutate(propertySimID = paste(task_id, property, start_density, sep = "-"),
         clusterSimID = paste(task_id, cluster, start_density, sep = "-"),
         projectSimID = paste(task_id, project, start_density, sep = "-")) |>
  add_ids()

cluster_observation_lookup <- take |>
  select(cluster, PPNum) |>
  distinct() |>
  mutate(observation_flag = 1)

property_observation_lookup <- take |>
  select(property, PPNum) |>
  distinct() |>
  mutate(observation_flag = 1)


select_pivot_longer <- function(df, node){
  df |>
    select(contains(node)) |>
    mutate(iter = 1:n()) |>
    pivot_longer(cols = -iter,
                 names_to = "node")
}

## analysis
samples_matrix <- out$posterior_samples
out$posterior_samples <- NULL

node_pattern <- "(?<=\\[)\\d*"

M_samples <- samples_matrix |>
  select_pivot_longer("M[") |>
  mutate(m_id = as.numeric(stringr::str_extract(node, node_pattern)))

message("Writing cluster abundance samples")
write_rds(M_samples, file.path(out_dir, "clusterAbundanceSamples.rds"))

N_samples <- samples_matrix |>
  select_pivot_longer("N[") |>
  mutate(n_id = as.numeric(stringr::str_extract(node, node_pattern)))

message("Writing property abundance samples")
write_rds(N_samples, file.path(out_dir, "propertyAbundanceSamples.rds"))


M_known <- known_abundance |>
  left_join(cluster_observation_lookup) |>
  mutate(observation_flag = if_else(is.na(observation_flag), 0, observation_flag)) |>
  select(m_id, M, clusterSimID, cluster, PPNum, observation_flag) |>
  distinct()

N_known <- known_abundance |>
  left_join(property_observation_lookup) |>
  mutate(observation_flag = if_else(is.na(observation_flag), 0, observation_flag)) |>
  select(n_id, N, M, property, propertySimID, clusterSimID, cluster, PPNum, observation_flag) |>
  distinct()

cluster_areas <- take |>
  select(cluster, cluster_area_km2) |>
  distinct()

property_areas <- take |>
  select(cluster, property, cluster_area_km2, property_area_km2) |>
  distinct()

M_join <- left_join(M_samples, M_known, by = "m_id") |>
  left_join(cluster_areas) |>
  mutate(estimated_density = value / cluster_area_km2,
         known_mean_density = M / cluster_area_km2,
         bias = estimated_density - known_mean_density)

N_join <- left_join(N_samples, N_known) |>
  left_join(property_areas) |>
  mutate(estimated_density = value / property_area_km2,
         known_mean_density = N / property_area_km2,
         bias = estimated_density - known_mean_density)

# normalized metrics will be undefined when mean known density = 0
m_summary <- function(df){
  df |>
    summarise(
      low = quantile(estimated_density, 0.05),
      med = quantile(estimated_density, 0.5),
      high = quantile(estimated_density, 0.95),
      med_known = median(known_mean_density),
      mean_bias = mean(bias),
      norm_mean_bias = mean_bias / mean(known_mean_density),
      rmse = sqrt(mean(bias^2)),
      norm_rmse = rmse /  mean(known_mean_density)
    ) |>
    ungroup() |>
    mutate(recovered_flag = if_else(med_known >= low & med_known <= high, 1, 0)) |>
    add_ids()
}

# M each cluster x primary period
M_by_time <- M_join |>
  group_by(cluster, clusterSimID, PPNum, M, known_mean_density, observation_flag) |>
  m_summary()

# M each cluster
M_by_cluster <- M_join |>
  group_by(cluster, clusterSimID) |>
  m_summary()

# N each property x primary period
N_by_time <- N_join |>
  group_by(property, propertySimID, clusterSimID, PPNum, N, observation_flag) |>
  m_summary()

# N each property
N_by_property <- N_join |>
  group_by(property, clusterSimID) |>
  m_summary()

# N across clusters
N_by_cluster <- N_join |>
  group_by(clusterSimID) |>
  m_summary()

out$M_by_time <- M_by_time
out$M_by_cluster <- M_by_cluster
out$N_by_time <- N_by_time
out$N_by_property <- N_by_property
out$N_by_cluster <- N_by_cluster

message("\nCorrelation: Cluster densities for each primary period")
round(cor(M_by_time$med_known, M_by_time$med), 2)

message("\nCorrelation: Mean cluster density")
round(cor(M_by_cluster$med_known, M_by_cluster$med), 2)

message("\nCorrelation: Property densities for each primary period")
round(cor(N_by_time$med_known, N_by_time$med), 2)

message("\nCorrelation: Mean property density")
round(cor(N_by_property$med_known, N_by_property$med), 2)

message("\nCorrelation: Mean property density for each cluster")
round(cor(N_by_cluster$med_known, N_by_cluster$med), 2)

message("\nPercent recovered: cluster density")
round(sum(M_by_time$recovered_flag) / nrow(M_by_time) * 100, 1)

message("\nPercent recovered: cluster density (observed)")
M_by_time_observed <- M_by_time |> filter(observation_flag == 1)
round(sum(M_by_time_observed$recovered_flag) / nrow(M_by_time_observed) * 100, 1)

message("\nPercent recovered: property density")
round(sum(N_by_time$recovered_flag) / nrow(N_by_time) * 100, 1)

message("\nPercent recovered: property density (observed)")
N_by_time_observed <- N_by_time |> filter(observation_flag == 1)
round(sum(N_by_time_observed$recovered_flag) / nrow(N_by_time_observed) * 100, 1)

### start here with parameter recovery

known_params <- read_csv("data/posterior_95CI_range_all.csv") |>
  suppressMessages() |>
  select(node, med) |>
  rename(actual = med)

alpha_cluster <- take |>
  select(cluster, alpha_cluster) |>
  distinct() |>
  mutate(node = paste0("alpha_cluster[", cluster, "]")) |>
  rename(actual = alpha_cluster) |>
  select(node, actual)

alpha_project <- take |>
  select(project, alpha_project) |>
  distinct() |>
  mutate(node = paste0("alpha_project[", project, "]")) |>
  rename(actual = alpha_project) |>
  select(node, actual)

tau_cluster <- take |>
  select(tau_cluster) |>
  distinct() |>
  mutate(node = "tau_cluster") |>
  rename(actual = tau_cluster) |>
  select(node, actual)

tau_project <- take |>
  select(tau_project) |>
  distinct() |>
  mutate(node = "tau_project") |>
  rename(actual = tau_project) |>
  select(node, actual)

params_actual <- bind_rows(
  known_params,
  alpha_cluster,
  alpha_project,
  tau_cluster,
  tau_project
) |>
  add_ids()

out$params_actual <- params_actual
out <- purrr::compact(out)
write_rds(out, file.path(out_dir, "simulationData.rds"))









