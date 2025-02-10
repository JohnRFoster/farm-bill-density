fit_mcmc <- function(cl, task_seed, modelCode, data, constants, start_density,
                     n_iter, n_chains, custom_samplers = NULL, monitors_add = NULL){

  require(nimble)
  require(parallel)
  source("R/functions_prep_nimble_simulation.R")

  export <- c(
    "modelCode",
    "data",
    "constants",
    "n_iter",
    "monitors_add",
    "custom_samplers"
  )

  clusterExport(cl, export, envir = environment())

  for(i in seq_along(cl)){
    set.seed(task_seed + i)
    init <- nimble_inits(data, constants, start_density)
    clusterExport(cl[i], "init", envir = environment())
  }

  message("Compiling model and initial parallel sampling...")
  start <- Sys.time()

  out <- clusterEvalQ(cl, {

    library(nimble)
    source("R/calc_log_potential_area.R")

    Rmodel <- nimbleModel(
      code = modelCode,
      constants = constants,
      data = data,
      inits = init,
      calculate = TRUE
    )

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

    Rmcmc <- buildMCMC(mcmcConf)
    Cmodel <- compileNimble(Rmodel)
    Cmcmc <- compileNimble(Rmcmc)

    runMCMC(
      Cmcmc,
      niter = n_iter,
      nchains = 1,
      nburnin = round(n_iter / 2),
      samplesAsCodaMCMC = TRUE
    )

  })

  message("Run time for ", n_iter, " iterations across ", n_chains, " chains in parallel:")
  print(Sys.time() - start)

  return(out)

}
