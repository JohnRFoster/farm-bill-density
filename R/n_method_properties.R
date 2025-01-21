# --------------------------------------------------------------------
#
# Generate psuedo-data for properties that use more than one method ----
# for removal model simulations
#
# John Foster
#
# --------------------------------------------------------------------

#' @param df data frame of MIS take values
#' @param n_properties the number of properties to simulate
#' @param n the number of methods to use at a property
#' @param n_pp the maximum number of primary periods to simulate
#' @return list of properties, their methods, and sample frequencies

n_method_properties <- function(df, n_properties, n, n_pp){
  require(dplyr)
  require(tidyr)
  require(purrr)

  # function to get property IDs given number of methods used ----
  get_properties <- function(n_method){
    df |>
      select(property, method) |>
      distinct() |>
      count(property) |>
      filter(n == n_method) |>
      pull(property)
  }

  n_method_props <- get_properties(n)

  # initiate list with the number of properties we need ----
  # list of things we need to keep track of
  property_attributes <- list(
    num = NULL,
    method_1 = NULL,
    method_2 = NULL,
    method_3 = NULL,
    method_4 = NULL,
    method_5 = NULL,
    area = NULL,
    effort = NULL
  )
  properties_n <- rep(list(property_attributes), n_properties)

  ## function to place vectors into each property attributes' list ----
  place_vec <- function(ls, n, where, what){
    require(purrr)
    ls |> map2(1:n, \(x, y) assign_in(x, where, what[y]))
  }

  properties_n <- place_vec(properties_n, n_properties, "num", 1:n_properties)

  ## count of method pairs by first method used then second ----
  ## for sampling from to assign 2-method properties their methods
  m_vec <- c("Traps", "Snares", "Firearms", "Fixed wing", "Helicopter")

  get_method_rel_freq <- function(prop_vec){
    temp <- df |>
      filter(property %in% prop_vec) |>
      select(property, method) |>
      distinct() |>
      group_by(property) |>
      mutate(m = paste0("m_", 1:n())) |>
      ungroup() |>
      pivot_wider(names_from = m,
                  values_from = method) |>
      select(starts_with("m_"))

    m <- colnames(temp)

    temp |> count(!!!syms(m)) |>
      mutate(prob = n / sum(n))
  }

  ### the probability of using a given method first for 2-method properties ----
  two_method_prob <- get_method_rel_freq(n_method_props)

  ### sample the index of rows in probability table given probability column
  pairs_sample <- sample(1:nrow(two_method_prob), n_properties, prob = two_method_prob$prob, replace = TRUE)

  ### use index samples to select method pairs
  sample_n_method <- two_method_prob |>
    slice(pairs_sample) |>
    select(starts_with("m_"))

  ## place the pair of methods in the 2-method data list ----
  for(i in 1:n){
    where <- paste0("method_", i)
    what <- sample_n_method |> pull(paste0("m_", i))
    properties_n <- place_vec(properties_n, n_properties, where, what)
  }

  ## function to sample property area from a given set of properties and methods ----
  subset_area_n_method <- function(prop_vec){
    df |>
      filter(property %in% prop_vec) |>
      select(property, method, property_area_km2) |>
      distinct() |>
      group_by(property) |>
      mutate(n = paste0("n_", 1:n())) |>
      ungroup() |>
      pivot_wider(names_from = n,
                  values_from = method)
  }

  get_area_n <- function(area_df, m){

    area_df |>
      pivot_longer(cols = starts_with("n_"),
                   names_to = "n",
                   values_to = "method") |>
      group_by(property) |>
      filter(all(method %in% m)) |>
      ungroup() |>
      select(property, property_area_km2) |>
      distinct() |>
      pull(property_area_km2) |>
      sample(1)

  }

  ## assign areas to properties ----
  n_method_areas <- subset_area_n_method(n_method_props)
  areas <- 1:n_properties |>
    map(\(x) get_area_n(n_method_areas,
                        unname(unlist(slice(sample_n_method, x)))
    )) |>
    list_c()

  properties_n <- place_vec(properties_n, n_properties, "area", round(areas, 2))

  ## joint and single return intervals ----
  n_return_intervals <- function(prop_vec, n){

    temp <- df |>
      filter(property %in% prop_vec) |>
      select(property, timestep, method) |>
      distinct() |>
      pivot_wider(names_from = method,
                  values_from = method) |>
      unite(method, -c(property, timestep), sep = "&", na.rm = TRUE) |>
      group_by(property, method) |>
      mutate(delta = c(0, diff(timestep))) |>
      # filter(delta > 0) |>
      ungroup() |>
      separate(method, paste0("method_", 1:n), sep = "&") |>
      suppressWarnings() # warns about not enough pieces, fill with NA is what we want

    n <- apply(temp[,paste0("method_", 1:n)], 1, function(x) sum(!is.na(x)))
    temp |> mutate(n = n)
  }

  all_return_intervals <- n_return_intervals(n_method_props, n)

  ## function to create a return interval df given two methods ----
  get_joint_return <- function(return_df, m){

    nm <- length(m)

    p <- return_df |>
      pivot_longer(cols = starts_with("method"),
                   names_to = "Mn",
                   values_to = "method",
                   values_drop_na = TRUE) |>
      select(property, method) |>
      distinct() |>
      group_by(property) |>
      filter(all(method %in% m)) |>
      ungroup() |>
      pull(property) |>
      unique()

    rdf <- return_df |>
      filter(property %in% p,
             delta > 0)

    missing_m <- rep(FALSE, 2)
    if(nm >= 2) missing_m[1] <- ifelse(all(is.na(rdf$method_2)), TRUE, FALSE)
    if(nm >= 3) missing_m[2] <- ifelse(all(is.na(rdf$method_3)), TRUE, FALSE)
    if(nm >= 4) missing_m[3] <- ifelse(all(is.na(rdf$method_3)), TRUE, FALSE)
    if(nm >= 5) missing_m[4] <- ifelse(all(is.na(rdf$method_3)), TRUE, FALSE)

    create_deltas <- function(df, p){
      temp <- df |>
        filter(property %in% p)

      z <- which(temp$delta == 0)
      temp$delta[z] <- sample.int(20, length(z), replace = TRUE)
      temp |> filter(delta > 0)
    }

    if(any(missing_m)){
      return(create_deltas(return_df, p))
    } else {
      return(rdf)
    }

  }

  ## function to generate sample occasions of single or joint use of methods
  get_sample_occasions_n <- function(joint_return_df, m, max_pp){
    ### function to assign a stating PP given method ----
    start <- sample.int(round(max_pp*0.8), 1)

    obs <- joint_return_df |>
      slice(sample.int(nrow(joint_return_df), 1))
    effort <- obs |>
      select(starts_with("method")) |>
      mutate(sample_occasions = start)

    create_effort_df <- function(start_effort, start_pp){
      effort <- start_effort
      end_pp <- start_pp
      for(xx in start_pp:max_pp){
        obs <- joint_return_df |>
          slice(sample.int(nrow(joint_return_df), 1))

        interval <- obs |> pull(delta)
        end_pp <- end_pp + interval

        m_new <- obs |>
          select(starts_with("method")) |>
          mutate(sample_occasions = end_pp)
        effort <- bind_rows(effort, m_new)

        # need to make sure we get at least two sample occasions
        if(nrow(effort) == 2 & end_pp >= max_pp){
          effort <- effort |> slice(1)
          end_pp <- effort |> pull(sample_occasions)
        }

        if(end_pp > max_pp){
          effort <- effort |> slice(-nrow(effort))
        }

      }
      effort
    }

    effort_sample <- create_effort_df(effort, start)

    # need to make sure all methods are used at least once
    n_method <- effort_sample |>
      select(starts_with("method_")) |>
      unlist() |>
      unique()
    not_all_used <- length(n_method[!is.na(n_method)]) < length(m)
    # not_all_used

    while(not_all_used){

      m_need <- m[which(!m %in% n_method)]
      obs_need <- joint_return_df |>
        pivot_longer(cols = starts_with("method_"),
                     names_to = "M",
                     values_to = "method",
                     values_drop_na = TRUE) |>
        filter(method %in% m_need) |>
        group_by(n, delta, timestep, property) |>
        mutate(M = paste0("method_", 1:n())) |>
        ungroup() |>
        pivot_wider(names_from = "M",
                    values_from = "method")

      obs <- obs_need |>
        slice(sample.int(nrow(obs_need), 1))

      occasions_possible <- 1:max_pp
      occasions_have <- effort_sample$sample_occasions
      occasions_pool <- occasions_possible[which(!occasions_possible %in% occasions_have)]

      if(length(occasions_pool) == 0){
        effort_sample <- effort_sample |> slice(seq(2, nrow(effort_sample), by = 2))
        occasions_have <- effort_sample$sample_occasions
        occasions_pool <- occasions_possible[which(!occasions_possible %in% occasions_have)]
      }

      new_occasion <- sample(occasions_pool, 1)

      m_new <- obs |>
        select(starts_with("method")) |>
        mutate(sample_occasions = new_occasion)
      effort_sample <- bind_rows(effort_sample, m_new) |>
        arrange(sample_occasions)
      n_method <- effort_sample |>
        select(starts_with("method_")) |>
        unlist() |>
        unique()
      not_all_used <- length(n_method[!is.na(n_method)]) < length(m)

    }
    effort_sample
  }

  get_n_reps_joint <- function(prop_vec){
    df |>
      filter(property %in% prop_vec) |>
      group_by(property, timestep, method) |>
      count() |>
      ungroup() |>
      pivot_wider(names_from = method,
                  values_from = n)
  }

  n_reps_method <- get_n_reps_joint(n_method_props)

  #### function to sample number of removal events in a PP given a pair of methods ----
  get_reps_n <- function(m){

    m <- m[complete.cases(m)]

    df <- n_reps_method |>
      select(all_of(m)) |>
      filter(if_all(all_of(m), ~ !is.na(.)))

    if(ncol(df) == 1){
      df <- df |> filter(if_all(all_of(m), ~ . >= 1))
      reps <- c(apply(df, 2, sample, 1), NA)
    } else {
      reps <- apply(df, 2, sample, 1)
    }
    reps_r <- t(as.data.frame(reps))
    colnames(reps_r) <- paste0("n_reps_", 1:ncol(reps_r))
    as.data.frame(reps_r)
  }


  ### sample the number of observations and reps, place in properties list
  pb <- txtProgressBar(max = n_properties, style = 1)
  for(i in seq_len(n_properties)){
    m <- unname(unlist(slice(sample_n_method, i)))

    sample_occasions <- get_joint_return(all_return_intervals, m) |>
      get_sample_occasions_n(m, n_pp)

    m_sample <- sample_occasions |>
      select(starts_with("method_"))

    n_reps <- 1:nrow(m_sample) |>
      map(\(x) get_reps_n(unname(unlist(slice(m_sample, x))))) |>
      list_rbind()

    effort <- bind_cols(sample_occasions, n_reps)

    testthat::expect(effort$sample_occasions[nrow(effort)] <= n_pp,
                     paste("Last PP exceeds boundary for 2-method property", i))

    properties_n <- assign_in(properties_n, list(i, "effort"), effort)
    setTxtProgressBar(pb, i)
  }
  close(pb)

  return(properties_n)

}
