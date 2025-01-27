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
#     - 2- to 5-method properties
#   - Simulate eco/take dynamics
#   - Fit MCMC
#   - Check MCMC
# 3. Summarize output
#
# --------------------------------------------------------------------

config_name <- "dev"
config <- config::get(config = config_name)
interval <- 4


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

# -----------------------------------------------------------------
# Load MIS data ----
# -----------------------------------------------------------------
message("MIS data intake")
data_dir <- config$data_dir
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
source("R/one_method_properties.R")
n_pp <- config$n_pp       # the number of primary periods to simulate
method_1 <- one_method_properties(df_with_timesteps, n_rel$n_simulate[1], n_pp)

# -----------------------------------------------------------------
# n-method properties ----
# -----------------------------------------------------------------

source("R/n_method_properties.R")

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

start_density <- config$start_density
cluster_props <- all_clusters

simulated_data <- simulate_cluster_dynamics(
  start_density = start_density,
  cluster_props = all_clusters,
  properties = properties
)

take <- simulated_data$all_take
abundance <-  simulated_data$all_pigs

nimble_ls <- prep_nimble(take)
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

n_chains <- 3
cl <- makeCluster(n_chains)

samples <- fit_mcmc(
  cl = cl,
  modelCode = modelCode,
  data = nimble_ls$data,
  constants = nimble_ls$constants,
  start_density = start_density,
  n_iter = 10000,
  n_chains = n_chains,
  custom_samplers = NULL,
  monitors_add
)

stopCluster(cl)

samples_mcmc <- as.mcmc.list(samples)

check <- check_mcmc(
  samples = samples_mcmc,
  nodes_check = params_check
)











source("R/functions_misc.R")
all_time_ids <- make_all_pp(take)

known_abundance <- left_join(all_time_ids, abundance)

samples <- as.matrix(samples)
hist(as.matrix(samples[,"M[201]"]))
known_abundance |> filter(m_id == 201)

all_time_ids |> filter(property == 1)
abundance |> filter(property == 1)
take |> filter(property == 1)



inits <- nimble_inits(nimble_data, nimble_constants, start_density)

source("R/nimble_cluster_model.R")
Rmodel <- nimbleModel(
  code = modelCode,
  constants = nimble_constants,
  data = nimble_data,
  inits = inits
)

Rmodel$initializeInfo()
