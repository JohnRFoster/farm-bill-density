
make_properites <- function(df, n_clusters, n_props_in_each_cluster, n_pp){

  cluster_sub <- df |>
    select(property, cluster) |>
    distinct() |>
    count(cluster) |>
    filter(n == n_props_in_each_cluster) |>
    pull(cluster)

  property_info <- df |>
    filter(cluster %in% cluster_sub) |>
    select(Project, propertyID, method, property_area_km2, cluster_area_km2) |>
    distinct()

  prop <- property_info |>
    group_by(propertyID) |>
    count(propertyID, name = "n_method") |>
    ungroup() |>
    count(n_method) |>
    mutate(p = n / sum(n))

  areas <- df |>
    filter(cluster %in% cluster_sub) |>
    pull(cluster_area_km2)

  p <- prop$p
  n_method <- prop$n_method
  properties <- list()

  for(i in seq_len(n_clusters)){

    n_method_sample <- sample(
      n_method,
      n_props_in_each_cluster,
      replace = TRUE,
      prob = p)

    draws <- sample(nrow(property_info), n_props_in_each_cluster)

    property_areas <- property_info$property_area_km2[draws]
    cluster_areas <- property_info$cluster_area_km2[draws]
    project <- property_info$Project[draws]

    for(j in seq_along(n_method_sample)){

      if(n_method_sample[j] == 1){

        sim <- one_method_properties(df, 1, n_pp)

      } else {

        sim <- n_method_properties(df, 1, n_method_sample[j], n_pp)

      }

      sim[[1]]$property_area <- property_areas[j]
      sim[[1]]$cluster_area <- cluster_areas[j]
      sim[[1]]$project <- project[j]

      properties <- c(properties, sim)

    }
  }

  return(properties)

}


