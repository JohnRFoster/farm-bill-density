
library(tidyverse)
library(geosphere)

data_repo <- "../pigs-statistical/data"

data_fb_start_end <- read_csv(file.path(data_repo, "FB_Properties_Start_End.csv"))
fb_start_end <- data_fb_start_end |>
  rename(alws_agrprop_id = PropID,
         agreementID = AgrID,
         STATE = State) |>
  mutate(fb_id = paste0(agreementID, "-", alws_agrprop_id)) |>
  select(-PropName, -AgrName)

n_ids <- length(unique(fb_start_end$fb_id))

data_latlon <- read_csv(file.path(data_repo, "FarmBill_GPS_1124.csv"))

data_latlon <- as.data.frame(data_latlon)
data_latlon[118, "Long"] <- -87.72853
data_latlon[906, "Long"] <- -157.8407943
data_latlon[912, "Long"] <- -157.8476406
data_latlon[3057, "Long"] <- -87.82308
data_latlon[3190, "Long"] <- -93.327563
data_latlon[3370, "Long"] <- -90.97315

fb_lat_lon <- data_latlon |>
  as_tibble() |>
  filter(!is.na(Long),
         !is.na(Lat)) |>
  select(-`Agr Name`, -`Property Name`)

which(is.na(as.numeric(fb_lat_lon$Lat)))

fb_lat_lon <- as.data.frame(fb_lat_lon)
fb_lat_lon[2132, "Lat"] <- "34.39114"
fb_lat_lon[2223, "Lat"] <- "34.23057"
fb_lat_lon[2254, "Lat"] <- "34.63701"
fb_lat_lon[2292, "Lat"] <- "34.88802"
fb_lat_lon[2305, "Lat"] <- "34.46177"

fb_lat_lon2 <- fb_lat_lon |>
  as_tibble() |>
  mutate(Lat = as.numeric(Lat)) |>
  rename(alws_agrprop_id = AgrPropId,
         agreementID = `Agreement ID`,
         STATE = State) |>
  mutate(fb_id = paste0(agreementID, "-", alws_agrprop_id),
         Lat = abs(Lat)) |>
  select(-Project) |>
  distinct() |>
  group_by(fb_id, Lat) |>
  summarise(Long = mean(Long)) |>
  ungroup()

# these properties have bad coordinates
fb_lat_lon2 |>
  filter(Lat <= 7 |
           Lat >= 85 |
           Long >= -20 |
           Long <= -179)

fb_lat_lon3 <- fb_lat_lon2 |>
  mutate(Long = if_else(Long > 0, -1 * Long, Long)) |>
  filter(Lat >= 7,
         Lat <= 85,
         Long <= -20,
         Long >= -179)

fb_join <- left_join(fb_start_end, fb_lat_lon3)

years_worked <- fb_join |>
  mutate(LastYr = if_else(is.na(LastYr), 2024, LastYr)) |>
  group_by(fb_id) |>
  filter(LastYr == min(LastYr)) |>
  distinct()

nrow(years_worked) == n_ids

n_projects <- years_worked |>
  group_by(fb_id) |>
  summarise(n_project = length(unique(Project))) |>
  ungroup()

single_proj <- n_projects |>
  filter(n_project == 1) |>
  pull(fb_id)

single_proj_data <- years_worked |>
  filter(fb_id %in% single_proj)

duplicate_proj <- n_projects |>
  filter(n_project > 1)

duplicates_removed <- years_worked |>
  filter(fb_id %in% duplicate_proj$fb_id) |>
  group_by(fb_id) |>
  mutate(p = 1:n()) |>
  ungroup() |>
  filter(p == 1) |>
  select(-p)

distinct_projects <- bind_rows(single_proj_data, duplicates_removed)

nrow(distinct_projects) == n_ids

write_csv(distinct_projects, file.path(data_repo, "FarmBill_GPS_1124_start_end_clean.csv"))


