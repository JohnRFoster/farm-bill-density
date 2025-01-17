
library(geosphere)
library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggspatial)
library(usmap)

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
  select(propertyID, Project, STATE, FIPS, property_area_km2, Long, Lat) |>
  distinct()

# =======================
# TODO
# check there isn't too much overlap in property ID and cluster ID
#    if there is a lot of overlap wont be identifable
# get landcover information for property latlon
# =======================

area_threshold <- max_area * 1

create_clusters <- function(df, a, cmin){

  latlon <- df |>
    select(Long, Lat) |>
    as.matrix()

  dist_matrix <- distm(latlon) / 1000
  tmp <- as.dist(dist_matrix)
  hc <- hclust(tmp, method = "complete")

  d <- 2 * sqrt(a / pi)

  clust <- cutree(hc, h = d)

  df$cluster <- clust + cmin

  df

}

get_max_size <- function(df){
  df |>
    group_by(cluster) |>
    summarise(area = sum(property_area_km2)) |>
    ungroup() |>
    filter(area == max(area)) |>
    pull(area)
}

get_bad_properties <- function(dfc){

  cc <- dfc |>
    group_by(cluster) |>
    summarise(area = sum(property_area_km2)) |>
    ungroup() |>
    filter(area > max_area) |>
    pull(cluster)

  dfc |>
    filter(cluster %in% cc) |>
    pull(propertyID)

}

large_properties <- gps_info |> filter(property_area_km2 >= area_threshold)
small_properties <- gps_info |> filter(property_area_km2 < area_threshold)

clusters1 <- create_clusters(small_properties, area_threshold, 0)

bad_props <- get_bad_properties(clusters1)

all_clusters <- clusters1 |> filter(!propertyID %in% bad_props)

dfc <- clusters1

max_size <- get_max_size(clusters1)
scaler <- seq(0.9, 0, length.out = 20)
for(s in scaler){

  tmp <- small_properties |> filter(propertyID %in% bad_props)

  mc <- dfc |> pull(cluster) |> max()

  dfc <- create_clusters(tmp, max_area * s, mc)

  bad_props <- get_bad_properties(dfc)
  dfg <- dfc |> filter(!propertyID %in% bad_props)

  if(nrow(dfg) == 0){
    next
  } else {
    all_clusters <- bind_rows(all_clusters, dfg)
  }

  if(nrow(all_clusters) == nrow(small_properties)) break

  if(s == 0){
    mc <- dfc |> pull(cluster) |> max()
    tmp <- small_properties |>
      filter(propertyID %in% bad_props) |>
      mutate(cluster = seq(mc + 1, by = 1, length.out = n()))

    all_clusters <- bind_rows(all_clusters, tmp)

  }

}

mc <- all_clusters |> pull(cluster) |> max()
large_properties$cluster <- seq(mc + 1, by = 1, length.out = nrow(large_properties))

all_clusters <- bind_rows(all_clusters, large_properties) |>
  group_by(STATE) |>
  mutate(state_cluster = cluster - min(cluster) + 1) |>
  ungroup() |>
  mutate(state_cluster = paste0(STATE, "-", state_cluster))

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


