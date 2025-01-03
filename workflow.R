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
  after_start[i] <- which(timestep_df$start_dates <= df$end.date[i]) |> max()
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
         n_simulate = ceiling(rel_prop * 110))

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

n_properties <- length(properties)
assertthat::are_equal(n_properties, sum(n_rel$n_simulate))



