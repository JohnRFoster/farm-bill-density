
#---------
#
# make spatial clusters given max distance (diameter of max area)
# make sure that the sum of property sizes within each cluster are less than max area
#
#---------

#'@description Group properties into spatial clusters given a max area
#'@param max_area max area that properties are grouped into
#'@param df data frame with columns propertyID, Project, STATE, FIPS, property_area_km2, Long, Lat


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
      pull(propertyID)

  }

  large_properties <- df |> filter(property_area_km2 >= area_threshold)
  small_properties <- df |> filter(property_area_km2 < area_threshold)

  clusters1 <- create_clusters(small_properties, area_threshold, 0)

  bad_props <- get_bad_properties(clusters1, area_threshold)

  if(length(bad_props) == 0){
    all_clusters <- clusters1
  } else {
    all_clusters <- clusters1 |> filter(!propertyID %in% bad_props)

    dfc <- clusters1

    scaler <- seq(0.9, 0, length.out = 20)
    for(s in scaler){

      tmp <- dfc |> filter(propertyID %in% bad_props)

      mc <- dfc |> pull(cluster) |> max()

      area_threshold <- max_area * s
      dfc <- create_clusters(tmp, area_threshold, mc)

      bad_props <- get_bad_properties(dfc, max_area)
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

    all_clusters <- bind_rows(all_clusters, large_properties)
  }

  all_clusters <- all_clusters |>
    mutate(cluster = as.numeric(as.factor(cluster))) |>
    group_by(STATE) |>
    mutate(state_cluster = cluster - min(cluster) + 1) |>
    ungroup() |>
    mutate(state_cluster = paste0(STATE, "-", state_cluster))

  cluster_areas <- all_clusters |>
    group_by(cluster) |>
    summarise(cluster_area_km2 = sum(property_area_km2),
              n_props = n()) |>
    ungroup()

  larger_than_max_area <- cluster_areas |>
    filter(cluster_area_km2 > max_area)

  assertthat::assert_that(all(larger_than_max_area$n_props == 1))

  return(all_clusters)
}


# ecological process
process_model <- function(N, zeta, a_phi, b_phi){
  n <- length(N)
  phi <- rbeta(n, a_phi, b_phi)
  lambda <- (N / 2) * zeta + (N * phi)
  rpois(n, lambda)
}
