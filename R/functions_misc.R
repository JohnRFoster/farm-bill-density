
#---------
#
# make spatial clusters given max distance (diameter of max area)
# make sure that the sum of property sizes within each cluster are less than max area
#
#---------

#'@description Group properties into spatial clusters given a max area
#'@param max_area max area that properties are grouped into
#'@param df data frame with columns property, Project, STATE, FIPS, property_area_km2, Long, Lat


make_clusters <- function(max_area, df){

  require(dplyr)

  area_threshold <- max_area * 1

  create_clusters <- function(df, a, cmin){

    latlon <- df |>
      select(Long, Lat) |>
      as.matrix()

    dist_matrix <- geosphere::distm(latlon) / 1000
    tmp <- as.dist(dist_matrix)
    hc <- hclust(tmp, method = "complete")

    d <- 2 * sqrt(a / pi)

    clust <- cutree(hc, h = d)

    df$cluster <- clust + cmin
    df$cluster_area_km2 <- a

    df

  }

  get_bad_properties <- function(dfc, max_area){

    cc <- dfc |>
      group_by(cluster) |>
      summarise(area = sum(property_area_km2)) |>
      ungroup() |>
      filter(area > max_area) |>
      pull(cluster)

    dfc |>
      filter(cluster %in% cc) |>
      pull(property)

  }

  large_properties <- df |> filter(property_area_km2 >= area_threshold)
  small_properties <- df |> filter(property_area_km2 < area_threshold)

  clusters1 <- create_clusters(small_properties, area_threshold, 0)

  bad_props <- get_bad_properties(clusters1, area_threshold)

  if(length(bad_props) == 0){
    all_clusters <- clusters1
  } else {
    all_clusters <- clusters1 |> filter(!property %in% bad_props)

    dfc <- clusters1

    scaler <- seq(0.9, 0, length.out = 20)
    for(s in scaler){

      tmp <- dfc |> filter(property %in% bad_props)

      mc <- dfc |> pull(cluster) |> max()

      if(s > 0){
        area_threshold <- max_area * s
        dfc <- create_clusters(tmp, area_threshold, mc)

        bad_props <- get_bad_properties(dfc, area_threshold)
        dfg <- dfc |> filter(!property %in% bad_props)

        if(nrow(dfg) == 0){
          next
        } else {
          all_clusters <- bind_rows(all_clusters, dfg)
        }
      }

      if(nrow(all_clusters) == nrow(small_properties)) break

      if(s == 0){
        mc <- dfc |> pull(cluster) |> max()
        tmp <- small_properties |>
          filter(property %in% bad_props)

        tmp$cluster <- seq(mc + 1, by = 1, length.out = nrow(tmp))
        tmp$cluster_area_km2 <- tmp$property_area_km2

        all_clusters <- bind_rows(all_clusters, tmp)

      }

    }
  }

  mc <- all_clusters |> pull(cluster) |> max()
  large_properties$cluster <- seq(mc + 1, by = 1, length.out = nrow(large_properties))
  large_properties$cluster_area_km2 <- large_properties$property_area_km2

  all_clusters <- bind_rows(all_clusters, large_properties)

  cluster_size <- all_clusters |>
    group_by(cluster) |>
    summarise(n_props = n()) |>
    ungroup()

  assertthat::assert_that(sum(cluster_size$n_props) == length(unique(df$property)),
                        msg = "Lost a property!")

  cluster_area_1 <- left_join(all_clusters, cluster_size) |>
    filter(n_props == 1) |>
    mutate(cluster_area_km2 = property_area_km2)

  cluster_area_n <- left_join(all_clusters, cluster_size) |>
    filter(n_props > 1)

  cluster_areas <- bind_rows(cluster_area_1, cluster_area_n) |>
    arrange(property) |>
    select(-n_props)

  check_size <- cluster_areas |>
    group_by(cluster) |>
    summarise(area = sum(property_area_km2)) |>
    left_join(cluster_areas) |>
    mutate(d = cluster_area_km2 - area)

  assertthat::assert_that(all(check_size$d >= 0))
  assertthat::assert_that(all(!is.na(cluster_areas$cluster_area_km2)))
  assertthat::assert_that(all(cluster_areas$cluster_area_km2 > 0))

  bad_clusters <- left_join(cluster_areas, cluster_size) |>
    mutate(drop_flag = if_else(n_props == 1 & cluster_area_km2 < 1.8, 1, 0))

  cluster_return <- all_clusters |>
    select(-cluster_area_km2) |>
    left_join(bad_clusters) |>
    mutate(cluster = as.numeric(as.factor(cluster)))

  larger_than_max_area <- left_join(cluster_areas, cluster_size) |>
    filter(cluster_area_km2 > max_area)

  assertthat::assert_that(all(larger_than_max_area$n_props == 1))

  return(cluster_return)
}


# ecological process
process_model <- function(N, zeta, a_phi, b_phi){
  n <- length(N)
  phi <- rbeta(n, a_phi, b_phi)
  lambda <- (N / 2) * zeta + (N * phi)
  rpois(n, lambda)
}

make_all_pp <- function(df){

  cluster_vec <- df |>
    pull(cluster) |>
    unique()

  all_timesteps <- tibble()
  for(i in seq_along(cluster_vec)){

    tmp <- df |>
      filter(cluster == cluster_vec[i])

    cluster_props <- tmp |>
      pull(property) |>
      unique()

    cluster_pps <- tmp |>
      pull(PPNum) |>
      range()

    # message("i ", i)
    # print(cluster_pps)

    for(j in seq_along(cluster_props)){
      tmpp <- tibble(
        cluster = cluster_vec[i],
        property = cluster_props[j],
        PPNum = cluster_pps[1]:cluster_pps[2]
      )

      all_timesteps <- bind_rows(all_timesteps, tmpp)

    }

  }

  all_cluster_pp <- all_timesteps |>
    select(cluster, PPNum) |>
    distinct() |>
    mutate(m_id = 1:n())

  all_property_pp <- all_timesteps |>
    arrange(property, PPNum) |>
    mutate(n_id = 1:n())

  all_time_ids <- left_join(all_property_pp, all_cluster_pp)

}

check_input <- function(x){
  assertthat::assert_that(
    all(x > 0),
    msg = "All input values must be positive!")
}

prec_2_sd <- function(prec){
  check_input(prec)
  return(sqrt(1 / prec))
}

prec_2_var <- function(prec){
  check_input(prec)
  return(1 / prec)
}

sd_2_prec <- function(sd){
  check_input(sd)
  return(sd ^ -2)
}

sd_2_var <- function(sd){
  check_input(sd)
  return(sd * sd)
}

var_2_prec <- function(v){
  check_input(v)
  return(1 / v)
}

var_2_sd <- function(v){
  check_input(v)
  return(sqrt(v))
}


