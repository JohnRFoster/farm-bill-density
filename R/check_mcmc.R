check_mcmc <- function(samples, nodes_check, n_mcmc = 1000, dest = NULL){

  require(coda)

  all_nodes <- colnames(samples[[1]])
  n_iter <- nrow(samples[[1]])

  j <- unlist(lapply(nodes_check, function(x) grep(x, all_nodes)))

  params <- samples[,j]

  total_iter <- nrow(params[[1]])
  GBR <- gelman.plot(params)
  burnin <- GBR$last.iter[tail(which(apply(GBR$shrink[, , 2] > 1.1, 1, any)), 1) + 1]
  message("Burnin: ", burnin)

  if(is.na(burnin)) burnin <- round(total_iter / 2)
  params_burnin <- window(params, start = burnin)

  message("Calculating PSRF...")
  psrf <- gelman.diag(params_burnin, multivariate = FALSE)
  print(psrf)

  message("Calculating effective samples...")
  effective_samples <- effectiveSize(params)
  print(effective_samples)

  if(!is.null(dest)){
    write_rds(list(params_burnin), file.path(dest, "parameterSamples.rds"))
  }

  samples_burn_mcmc <- window(samples, start = burnin)
  samples_burn_mat <- as.matrix(samples_burn_mcmc)
  draws <- sample.int(nrow(samples_burn_mat), n_mcmc, replace = TRUE)

  samples_draw <- as.data.frame(samples_burn_mat[draws, ])

  list(
    posterior_samples = samples_draw,
    psrf = psrf$psrf,
    effective_samples = effective_samples,
    burnin = burnin,
    converged = max(psrf$psrf[, "Upper C.I."]) <= 1.1,
    bad_mcmc = any(is.na(psrf$psrf)) | any(psrf$psrf > 3)
  )

}
