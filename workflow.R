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

config_name <- "hpc_production"
config <- config::get(config = config_name)
interval <- 4

if(grepl("hpc", config_name)){
  Sys.setenv(RENV_CONFIG_SANDBOX_ENABLED = FALSE)
  renv::load("/home/john.foster/pigs-simulation/")

  args <- commandArgs(trailingOnly = TRUE)
  task_id <- args[1]

} else {
  task_id <- 1
}

message("Task ID: ", task_id)


library(nimble)
library(parallel)
library(coda)
library(readr)
library(dplyr)
library(tidyr)
library(purrr)
library(lubridate)

source("R/functions_misc.R")
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
df <- read_csv("data/farmBillTakeGoodStates.csv") |>
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
  arrange(propertyID, timestep)

# need a lookup table for property IDs and how may methods they employ ----
n_method_lookup <- df_with_timesteps |>
  select(property, method) |>
  distinct() |>
  count(property)

# get the proportion of properties that use n methods
# for determining the number of properties to simulate
n_rel <- n_method_lookup |>
  count(n, name = "n_sum") |>
  mutate(rel_prop = n_sum / sum(n_sum),
         n_simulate = ceiling(rel_prop * 200))

# -----------------------------------------------------------------
# 1-method properties ----
# -----------------------------------------------------------------
message("Create 1-method properties")
n_pp <- config$n_pp       # the number of primary periods to simulate
method_1 <- one_method_properties(df_with_timesteps, n_rel$n_simulate[1], n_pp)

# -----------------------------------------------------------------
# n-method properties ----
# -----------------------------------------------------------------


start_density <- round(runif(1, 0.3, 5), 3)
message("Starting density: ", start_density)

n_method <- 0
while(n_method != 5){

  message("Create 2-method properties")
  method_2 <- n_method_properties(df_with_timesteps, n_rel$n_simulate[2], 2, n_pp)
  message("Create 3-method properties")
  method_3 <- n_method_properties(df_with_timesteps, n_rel$n_simulate[3], 3, n_pp)
  message("Create 4-method properties")
  method_4 <- n_method_properties(df_with_timesteps, n_rel$n_simulate[4], 4, n_pp)

  properties <- c(
    method_1,
    method_2,
    method_3,
    method_4
  )

  # randomly assign each property a Lat Lon
  gps_info <- df_with_timesteps |>
    select(Lat, Long) |>
    slice(sample.int(nrow(df_with_timesteps), length(properties))) |>
    mutate(property = 1:n(),
           property_area_km2 = list_c(map(properties, "area")),
           STATE = "CO",
           propertyID = 1:n())

  all_clusters <- make_clusters(250, gps_info)

  n_clusters <- length(unique(all_clusters$state_cluster))

  n_per_cluster <- all_clusters |>
    group_by(state_cluster) |>
    count()

  summary(n_per_cluster$n)
  hist(n_per_cluster$n, breaks = 20)


  cluster_props <- all_clusters


  message("\nSimulate cluster dynamics")
  simulated_data <- simulate_cluster_dynamics(
    start_density = start_density,
    cluster_props = all_clusters,
    properties = properties
  )

  take <- simulated_data$all_take
  abundance <-  simulated_data$all_pigs

  n_method <- length(unique(take$method))

}

message("Prep data for MCMC")
nimble_ls <- prep_nimble(take, config$posterior_path) |> suppressMessages()
nimble_data <- nimble_ls$data
nimble_constants <- nimble_ls$constants

monitors_add <- c("N", "M")

params_check <- c(
  "beta_p",
  "beta1",
  "log_gamma",
  "log_rho",
  "phi_mu",
  "psi_phi",
  "log_nu",
  "p_mu"
)

n_chains <- config$n_chains
cl <- makeCluster(n_chains)

samples <- fit_mcmc(
  cl = cl,
  task_seed = task_id,
  modelCode = modelCode,
  data = nimble_ls$data,
  constants = nimble_ls$constants,
  start_density = start_density,
  n_iter = config$n_iter,
  n_chains = n_chains,
  custom_samplers = NULL,
  monitors_add
)

stopCluster(cl)

out_dir <- file.path(config$project_dir, config$out_dir, config$dev_dir, task_id)
message("\n\nWriting to: ", out_dir)

if(!dir.exists(out_dir)) dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

samples_mcmc <- as.mcmc.list(samples)

out <- check_mcmc(
  samples = samples_mcmc,
  nodes_check = params_check,
  dest = out_dir
)

source("R/functions_misc.R")
all_time_ids <- make_all_pp(take)
known_abundance <- left_join(all_time_ids, abundance) |>
  mutate(propertySimID = paste(task_id, property, start_density, sep = "-"),
         clustersimID = paste(task_id, cluster, start_density, sep = "-"))

out$known_abundance <- known_abundance
out$data <- nimble_ls$data
out$constants <- nimble_ls$constants
out$start_density <- start_density
out$task_id <- task_id
out$take <- take |>
  mutate(propertySimID = paste(task_id, property, start_density, sep = "-"),
         clusterSimID = paste(task_id, cluster, start_density, sep = "-"))

cluster_observation_lookup <- take |>
  select(cluster, PPNum) |>
  distinct() |>
  mutate(observation_flag = 1)

property_observation_lookup <- take |>
  select(property, propertyID, PPNum) |>
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

M_samples <- samples_matrix |>
  select_pivot_longer("M[") |>
  mutate(m_id = as.numeric(stringr::str_extract(node, "(?<=\\[)\\d*")))

message("Writing cluster abundance samples")
write_rds(M_samples, file.path(out_dir, "clusterAbundanceSamples.rds"))

N_samples <- samples_matrix |>
  select_pivot_longer("N[") |>
  mutate(n_id = as.numeric(stringr::str_extract(node, "(?<=\\[)\\d*")))

message("Writing property abundance samples")
write_rds(N_samples, file.path(out_dir, "propertyAbundanceSamples.rds"))

M_known <- known_abundance |>
  left_join(cluster_observation_lookup) |>
  mutate(observation_flag = if_else(is.na(observation_flag), 0, observation_flag)) |>
  select(m_id, M, clustersimID, cluster, PPNum, observation_flag) |>
  distinct()

N_known <- known_abundance |>
  left_join(property_observation_lookup) |>
  mutate(observation_flag = if_else(is.na(observation_flag), 0, observation_flag)) |>
  select(n_id, N, M, property, propertySimID, clustersimID, cluster, PPNum, observation_flag) |>
  distinct()

cluster_areas <- take |>
  select(cluster, cluster_area_km2) |>
  distinct()

property_areas <- take |>
  select(cluster, property, cluster_area_km2, property_area_km2) |>
  distinct()

# need to add abundance metrics

M_join <- left_join(M_samples, M_known) |>
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
    mutate(recovered_flag = if_else(med_known >= low & med_known <= high, 1, 0))
}

# M each cluster x primary period
M_by_time <- M_join |>
  group_by(cluster, clustersimID, PPNum, M, known_mean_density, observation_flag) |>
  m_summary()

# M each cluster
M_by_cluster <- M_join |>
  group_by(cluster, clustersimID) |>
  m_summary()

# N each property x primary period
N_by_time <- N_join |>
  group_by(property, propertySimID, clustersimID, PPNum, N, observation_flag) |>
  m_summary()

# N each property
N_by_property <- N_join |>
  group_by(property, clustersimID) |>
  m_summary()

# N across clusters
N_by_cluster <- N_join |>
  group_by(clustersimID) |>
  m_summary()

out$M_by_time <- M_by_time
out$M_by_cluster <- M_by_cluster
out$N_by_time <- N_by_time
out$N_by_property <- N_by_property
out$N_by_cluster <- N_by_cluster

out <- purrr::compact(out)
write_rds(out, file.path(out_dir, "simulationData.rds"))

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


