# --------------------------------------------------------------------
#
# Generate psuedo-data for properties that use a single method ----
# for removal model simulations
#
# John Foster
#
# --------------------------------------------------------------------

#' @param df data frame of MIS take values
#' @param n_properties the number of properties to simulate
#' @param n_pp the maximum number of primary periods to simulate
#' @return list of properties, their methods, and sample frequencies

one_method_properties <- function(df, n_properties, n_pp){
  require(dplyr)

  # function to get property IDs given number of methods used ----
  get_properties <- function(n_method){
    df |>
      select(property, method) |>
      distinct() |>
      count(property) |>
      filter(n == n_method) |>
      pull(property)
  }

  one_method_props <- get_properties(1)

  # list of things we need to keep track of
  property_attributes <- list(
    num = NULL,
    method_1 = NULL,
    area = NULL,
    effort = NULL
  )

  # initiate list with the number of properties we need ----
  properties_one <- rep(list(property_attributes), n_properties)

  ## function to place vectors into each property attributes' list ----
  place_vec <- function(ls, n, where, what){
    require(purrr)
    ls |> map2(1:n, \(x, y) assign_in(x, where, what[y]))
  }

  properties_one <- place_vec(properties_one, n_properties, "num", 1:n_properties)

  ## relative proportion of properties that use a specific method ----
  ## for sampling from to assign 1-method properties a method
  one_method_prob <- df |>
    filter(property %in% one_method_props) |>
    select(property, method) |>
    distinct() |>
    count(method) |>
    mutate(prob = n / sum(n))

  ## assign methods to 1-method properties ----
  sample_one_method <- sample(one_method_prob$method,
                              n_properties,
                              one_method_prob$prob,
                              replace = TRUE)

  properties_one <- place_vec(properties_one, n_properties, "method_1", sample_one_method)

  ## function to sample property area from a given set of properties and methods ----
  get_area <- function(dat, m){
    dat |>
      filter(method == m) |>
      pull(property.size) |>
      sample(1)
  }

  one_method_areas <- df |>
    filter(property %in% one_method_props) |>
    select(property, method, property.size) |>
    distinct()

  ## assign areas to 1-method properties ----
  areas <- sample_one_method |>
    map_vec(\(x) get_area(one_method_areas, x))

  properties_one <- place_vec(properties_one, n_properties, "area", round(areas, 2))

  ## now we need to determine which PPs are sampled for each property ----
  ### start with determining the first PP ----

  ### function to assign a stating PP given method ----
  start_pp <- function() sample.int(round(n_pp*0.8), n_properties, replace = TRUE)
  start <- start_pp()

  ### lookup table for return intervals by method for single method-properties ----
  one_method_return <- df |>
    filter(property %in% one_method_props) |>
    select(property, timestep, method) |>
    distinct() |>
    group_by(property, method) |>
    mutate(delta = c(0, diff(timestep)),
           single_occurance = if_else(n() == 1, 1, 0),
           first_occurance = if_else(delta == 0 & single_occurance == 0, 1, 0))

  #### function to assign sampling occasions ----
  get_sample_occasions <- function(return_df, m, s, max_pp){
    sample_occasions <- s
    end_pp <- s
    while(end_pp <= max_pp){
      interval <- return_df |>
        filter(method %in% m,
               first_occurance == 0) |>
        pull(delta) |>
        sample(1)
      end_pp <- end_pp + interval
      sample_occasions <- c(sample_occasions, end_pp)

      # need to make sure we get at least two sample occasions
      if(length(sample_occasions) == 2 & end_pp > max_pp){
        sample_occasions <- sample_occasions[1]
        end_pp <- sample_occasions
      }
    }
    sample_occasions[sample_occasions <= max_pp] |> as.vector()
  }

  ### lookup table for the number of repeat events in a PP by method for single-method properties ----
  n_reps_single_method <- df |>
    filter(property %in% one_method_props) |>
    group_by(property, timestep, method) |>
    count() |>
    ungroup()

  #### function to sample number of removal events in a PP given method
  get_reps <- function(reps_df, m, size){
    reps_df |>
      filter(method == m) |>
      pull(n) |>
      sample(size, replace = TRUE)
  }

  ### sample the number of observations and reps, place in properties list
  pb <- txtProgressBar(min = 1, max = n_properties, style = 1)
  for(i in seq_len(n_properties)){
    sample_occasions <- get_sample_occasions(one_method_return, sample_one_method[i], start[i], n_pp)
    n_reps <- get_reps(n_reps_single_method, sample_one_method[i], length(sample_occasions))

    effort <- tibble(sample_occasions = sample_occasions,
                     n_reps = n_reps,
                     method_1 = sample_one_method[i])

    testthat::expect(effort$sample_occasions[nrow(effort)] <= n_pp,
                     paste("Last PP exceeds boundary for 1-method property", i))

    properties_one <- assign_in(properties_one, list(i, "effort"), effort)
    setTxtProgressBar(pb, i)
  }
  close(pb)

  return(properties_one)

}

