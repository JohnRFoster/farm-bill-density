
library(geosphere)
library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggspatial)
library(usmap)

source("R/functions_misc.R")

data_repo <- "../pigs-statistical/data"
fb_info <- read_csv(file.path(data_repo, "FarmBill_GPS_1124_start_end_clean.csv"))

data_farm_bill <- read_csv(file.path(data_repo, "All_FB_Agreements_long_2024-05-30.csv"))
farm_bill_properties <- data_farm_bill |>
  rename(alws_agrprop_id = propertyID) |>
  select(-agreement_name, -property_name, -STATE) |>
  mutate(farm_bill = 1,
         fb_id = paste0(agreementID, "-", alws_agrprop_id))

fb_join <- left_join(fb_info, farm_bill_properties)

file <- file.path(data_repo, "insitu/MIS.Effort.Take.All.Methods.Daily.Events.01JUL2024.csv")

all_take <- read_csv(file, show_col_types = FALSE) |>
  select(-`...1`) |>
  filter(start.date >= lubridate::ymd("2019-01-01")) |>
  mutate(cnty_name = if_else(grepl("ST ", cnty_name), gsub("ST ", "ST. ", cnty_name), cnty_name),
         cnty_name = if_else(grepl("KERN", cnty_name), "KERN", cnty_name)) |>
  mutate(property_area_km2 = round(property.size * 0.00404686, 2)) |>
  mutate(effort = if_else(cmp_name %in% c("TRAPS, CAGE", "SNARE"), cmp.days, cmp.hours),
         effort_per = effort / cmp.qty,
         cmp_name = if_else(cmp_name == "TRAPS, CAGE", "TRAPS", cmp_name)) |>
  rename(method = cmp_name,
         trap_count = cmp.qty) |>
  select(-wt_work_date, -hours, -cmp.hours, -cmp.days) |>
  mutate(propertyID = paste0(agrp_prp_id, "-", alws_agrprop_id)) |>
  arrange(propertyID, start.date, end.date) |>
  rename(statefp = st_gsa_state_cd,
         countyfp = cnty_gsa_cnty_cd) |>
  mutate(countyfp = sprintf("%03d", countyfp),
         countyfp = ifelse(cnty_name == "HUMBOLDT (E)", "013", countyfp),
         FIPS = as.numeric(paste0(statefp, countyfp)),
         FIPS = sprintf("%05d", FIPS))

states <- read_csv("../pigs-statistical/data/counties/statePostalCodes.csv")
postal <- states |>
  mutate(State = toupper(State)) |>
  select(Postal, State) |>
  rename(st_name = State)

all_take <- left_join(all_take, postal)

states_model <- c("ARKANSAS", "HAWAII", "LOUISIANA", "MISSOURI", "OKLAHOMA",
                  "SOUTH CAROLINA", "TEXAS")

fb_take <- left_join(fb_join, all_take) |>
  filter(!is.na(take),
         st_name %in% states_model)

# only GA and FL are mismatched - suggests error in lat long reporting
fb_take |>
  filter(Postal != STATE) |>
  group_by(Postal, STATE) |>
  count()

max_area <- 250

all_lat_lon <- fb_take |>
  filter(!is.na(Long),
         !is.na(Lat))

write_csv(all_lat_lon, "data/farmBillTakeGoodStates.csv")

# ignoring project for now
gps_info <- all_lat_lon |>
  select(Project, propertyID, Project, STATE, FIPS, property_area_km2, Long, Lat) |>
  distinct()

# =======================
# TODO
# check there isn't too much overlap in property ID and cluster ID
#    Projects have multiple clusters
#    Most clusters belong to a single project

all_clusters <- make_clusters(max_area, gps_info)

write_csv(all_clusters, "data/clusters250km2.csv")


projs <- all_clusters |> pull(Project) |> unique()
for(i in seq_along(projs)){
  tmp <- all_clusters |>
    filter(Project == projs[i]) |>
    pull(cluster) |>
    unique()

  message("Project ", projs[i])
  print(tmp)

}


assertthat::are_equal(nrow(all_clusters), nrow(gps_info))

cluster_areas <- all_clusters |>
  group_by(cluster, state_cluster) |>
  summarise(cluster_area_km2 = sum(property_area_km2),
            n_props = n()) |>
  ungroup()

larger_than_max_area <- cluster_areas |>
  filter(cluster_area_km2 > max_area)

assertthat::assert_that(all(larger_than_max_area$n_props == 1))

n_clusters <- length(unique(all_clusters$state_cluster))

n_per_cluster <- all_clusters |>
  group_by(state_cluster) |>
  count()

summary(n_per_cluster$n)

n_per_county <- all_clusters |>
  group_by(FIPS) |>
  count()

hist(n_per_cluster$n, breaks = 25, main = "Cluster size", xlab = "Number of properties in a cluster")

all_clusters |>
  select(cluster, STATE) |>
  distinct() |>
  group_by(cluster) |>
  count(STATE) |>
  filter(n > 1)

all_clusters |>
  select(cluster, FIPS) |>
  distinct() |>
  group_by(cluster) |>
  count(FIPS) |>
  filter(n > 1)


