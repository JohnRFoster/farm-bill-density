
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

write_csv(fb_take, "data/farmBillTakeGoodStates.csv")


# only GA and FL are mismatched - suggests error in lat long reporting
fb_take |>
  filter(Postal != STATE) |>
  group_by(Postal, STATE) |>
  count()

# ignoring project for now
gps_info <- fb_take |>
  filter(!is.na(Long),
         !is.na(Lat)) |>
  select(propertyID, STATE, FIPS, Long, Lat) |>
  distinct()

IDs <- gps_info$propertyID

latlon <- gps_info |>
  select(Long, Lat) |>
  as.matrix()

length(unique(fb_take$propertyID))

dist_matrix <- distm(latlon) / 1000
tmp <- as.dist(dist_matrix)

hc <- hclust(tmp)

max_area <- 150
d <- 2 * sqrt(max_area / pi)
d <- d + 4

clust <- cutree(hc, h = d)

gps_info$cluster <- clust

gps_info_df <- gps_info |>
  group_by(STATE) |>
  mutate(state_cluster = cluster - min(cluster) + 1) |>
  ungroup() |>
  mutate(state_cluster = paste0(STATE, "-", state_cluster))


left_join(fb_take, gps_info)


n_clusters <- length(unique(gps_info_df$state_cluster))

n_per_cluster <- gps_info_df |>
  group_by(state_cluster) |>
  count()

n_per_county <- gps_info_df |>
  group_by(FIPS) |>
  count()

n_min2_per_cluster <- n_per_cluster |>
  filter(n >= 2) |>
  nrow()

n_1_per_cluster <- n_per_cluster |>
  filter(n == 1) |>
  nrow()

summary(n_per_cluster$n)
hist(n_per_cluster$n, breaks = 25, main = "Cluster size", xlab = "Number of properties in a cluster")

# check that clusters are within the correct distance
colnames(dist_matrix) <- IDs
dist_matrix[!upper.tri(dist_matrix)] <- -1

gps_info_df2 <- gps_info_df |>
  rename(propertyID2 = propertyID,
         STATE2 = STATE,
         FIPS2 = FIPS,
         state_cluster2 = state_cluster) |>
  select(-Long, -Lat, -cluster)

neighborhoods <- dist_matrix |>
  as_tibble() |>
  mutate(propertyID = IDs) |>
  pivot_longer(cols = -propertyID,
               names_to = "propertyID2",
               values_to = "distance") |>
  filter(distance != -1) |>
  left_join(gps_info_df) |>
  left_join(gps_info_df2) |>
  filter(state_cluster == state_cluster2) |>
  select(-state_cluster2)

neighborhoods |>
  group_by(state_cluster) |>
  filter(distance == max(distance))

# no clusters cross state lines
neighborhoods |>
  filter(STATE != STATE2)

# 69 clusters cross county lines
neighborhoods |>
  filter(FIPS != FIPS2) |>
  pull(state_cluster) |>
  unique() |>
  length()





# GeoLocations <- usmap_transform(
#   as.data.frame(gps_info_df),
#   input_names = c("Long", "Lat"))
#
# include <- c("OK")#, "OK", "LA", "MO", "AR", "MS", "AL", "GA", "NC", "SC")
#
#
#
# plot_data <- GeoLocations |>
#   filter(STATE %in% include) |>
#   filter(propertyID != "423698-420737")
#
# plot_usmap(regions = "county",
#            include = include
#            ) +
#   geom_sf(
#     data = plot_data,
#     aes(color = state_cluster)
#   ) +
#   theme(legend.position = "none")
