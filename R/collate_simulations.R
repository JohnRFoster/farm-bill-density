library(dplyr)
library(readr)
config_name <- "hpc_production"
config <- config::get(config = config_name)

# simulations are stored here
out_dir <- file.path(config$project_dir, config$out_dir, config$dev_dir)

sim_files <- list.files(out_dir)

converged <- rep(NA, length(sim_files))

all_take <- tibble()
all_M_by_time <- tibble()
all_M_by_cluster <- tibble()
all_N_by_time <- tibble()
all_N_by_property <- tibble()
all_N_by_cluster <- tibble()

pb <- txtProgressBar(max = length(sim_files), style = 1)

for(i in seq_along(sim_files)){

  task_id <- sim_files[i]

  out_file <- file.path(out_dir, task_id, "simulationData.rds")

  if(file.exists(out_file)){
    out_data <- read_rds(out_file)


    psrf <- out_data$psrf
    psrf <- psrf[row.names(psrf) != "psi_phi",]
    psrf <- psrf[row.names(psrf) != "phi_mu",]

    good <- all(psrf[,1] <= 1.1)

    if(good){
      start_density <- out$start_density

      take <- out_data$take |> mutate(start_density = start_density, simulation = task_id)
      all_take <- bind_rows(all_take, take)

      M_by_time <- out_data$M_by_time |> mutate(start_density = start_density, simulation = task_id)
      all_M_by_time <- bind_rows(all_M_by_time, M_by_time)

      M_by_cluster <- out_data$M_by_cluster |> mutate(start_density = start_density, simulation = task_id)
      all_M_by_cluster <- bind_rows(all_M_by_cluster, M_by_cluster)

      N_by_time <- out_data$N_by_time |> mutate(start_density = start_density, simulation = task_id)
      all_N_by_time <- bind_rows(all_N_by_time, N_by_time)

      N_by_property <- out_data$N_by_property |> mutate(start_density = start_density, simulation = task_id)
      all_N_by_property <- bind_rows(all_N_by_property, N_by_property)

      N_by_cluster <- out_data$N_by_cluster |> mutate(start_density = start_density, simulation = task_id)
      all_N_by_cluster <- bind_rows(all_N_by_cluster, N_by_cluster)

    }
  }
  setTxtProgressBar(pb, i)
}
close(pb)

dest <- file.path(config$project_dir, config$analysis_dir)

write_rds(all_take, file.path(dest, "all_take.rds"))
write_rds(all_M_by_time, file.path(dest, "all_M_by_time.rds"))
write_rds(all_M_by_cluster, file.path(dest, "all_M_by_cluster.rds"))
write_rds(all_N_by_time, file.path(dest, "all_N_by_time.rds"))
write_rds(all_N_by_property, file.path(dest, "all_N_by_property.rds"))
write_rds(all_N_by_cluster, file.path(dest, "all_N_by_cluster.rds"))



