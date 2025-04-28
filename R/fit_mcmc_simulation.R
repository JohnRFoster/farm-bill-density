fit_mcmc <- function(cl, task_seed, modelCode, data, constants, start_density,
                     n_iter, n_chains, custom_samplers = NULL, monitors_add = NULL,
                     include_project, include_cluster){

  require(nimble)
  require(parallel)
  source("R/functions_prep_nimble_simulation.R")

  base <- if_else(!include_project & !include_cluster, TRUE, FALSE)
  include_project_cluster <- if_else(include_project & include_cluster, TRUE, FALSE)
  project_prior <- if_else(include_project, TRUE, FALSE)
  cluster_prior <- if_else(include_cluster, TRUE, FALSE)

  export <- c(
    "modelCode",
    "data",
    "constants",
    "n_iter",
    "monitors_add",
    "custom_samplers",
    "include_project",
    "include_cluster",
    "include_project_cluster",
    "base",
    "project_prior",
    "cluster_prior"
  )

  clusterExport(cl, export, envir = environment())

  for(i in seq_along(cl)){
    set.seed(as.numeric(task_seed) + i)
    init <- nimble_inits(data, constants, start_density, buffer = 100, include_cluster, include_project)
    clusterExport(cl[i], "init", envir = environment())
  }

  message("Compiling model and initial parallel sampling...")
  start <- Sys.time()

  out <- clusterEvalQ(cl, {

    library(nimble)
    source("R/calc_log_potential_area.R")

    if(include_project_cluster){
      include_project <- FALSE
      include_cluster <- FALSE
    }

    Rmodel <- nimbleModel(
      code = modelCode,
      constants = constants,
      data = data,
      inits = init,
      calculate = TRUE
    )

    for(i in 1:constants$n_clusters){
      for(t in 1:constants$n_time_clust[i]){
        M_model <- Rmodel$M[constants$mH[i, t]]
        Z <- M_model - constants$rem[i, t]

        if(Z <= 0){
          message("cluster ", i, "pp ", t)
          n <- ifelse(Z == 0, 2, Z)
          Rmodel$M[constants$mH[i, t]] <- M_model + n^2
        }

      }

      n <- round(M_model - constants$y_sum[i])

    }

    # Rmodel$initializeInfo()
    # Rmodel$simulate()
    # Rmodel$calculate()

    # default MCMC configuration
    mcmcConf <- configureMCMC(Rmodel, useConjugacy = TRUE)

    if(!is.null(monitors_add)){
      mcmcConf$addMonitors(monitors_add)
    }

    if(!is.null(custom_samplers)){
      for(i in seq_len(nrow(custom_samplers))){
        node <- custom_samplers$node[i]
        type <- custom_samplers$type[i]
        mcmcConf$removeSampler(node)
        mcmcConf$addSampler(node, type)
      }
    }

    for(i in 1:5){
      node <- paste0("beta_p[", i, ", ", 1:constants$m_p, "]")
      node <- c(paste0("beta1[", i, "]"), node)
      mcmcConf$removeSampler(node)
      mcmcConf$addSampler(node, "AF_slice")
    }

    # mcmcConf$addSampler(target = c("alpha_project", "tau_project"), type = "crossLevel")
    # mcmcConf$addSampler(target = c("alpha_cluster", "tau_cluster"), type = "crossLevel")

    mcmcConf$printSamplers(byType = TRUE)

    Rmcmc <- buildMCMC(mcmcConf)
    Cmodel <- compileNimble(Rmodel)
    Cmcmc <- compileNimble(Rmcmc)

    runMCMC(
      Cmcmc,
      niter = n_iter,
      nchains = 1,
      thin = 10,
      samplesAsCodaMCMC = TRUE
    )

  })

  message("Run time for ", n_iter, " iterations across ", n_chains, " chains in parallel:")
  print(Sys.time() - start)

  return(out)

}
