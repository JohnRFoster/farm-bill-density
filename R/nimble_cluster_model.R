
library(nimble)

source("R/calc_log_potential_area.R")

modelCode <- nimbleCode({

  # priors
  for(i in 1:n_method){
    log_rho[i] ~ dnorm(log_rho_mu[i], tau = log_rho_tau[i])
  }

  for(i in 1:2){
    p_mu[i] ~ dnorm(p_mu_mu[i], tau = p_mu_tau[i])
    logit(p_unique[i]) <- p_mu[i]

    log_gamma[i] ~ dnorm(log_gamma_mu[i], tau = log_gamma_tau[i])
  }

  # non time varying coefficients - observation model
  for(i in 1:n_method){
    beta1[i] ~ dnorm(beta1_mu[i], tau = beta1_tau[i])
  }

  for(i in 1:n_betaP){
    beta_p[beta_p_row[i], beta_p_col[i]] ~ dnorm(beta_p_mu[i], tau = beta_p_tau[i])
  }

  # project random effect - in data model
  if(project_prior){

    tau_project ~ dgamma(1, 1)

    for(i in 1:n_projects){
      alpha_project[i] ~ dnorm(0, tau = tau_project)
    }

  }

  # cluster random effect - in data model
  if(cluster_prior){

    tau_cluster ~ dgamma(1, 1)

    for(i in 1:n_clusters){
      alpha_cluster[i] ~ dnorm(0, tau = tau_cluster)
    }

  }

  # estimate apparent survival
  phi_mu ~ dbeta(phi_mu_a, phi_mu_b)
  psi_phi ~ dgamma(psi_shape, psi_rate)
  a_phi <- phi_mu * psi_phi
  b_phi <- (1 - phi_mu) * psi_phi

  log_nu ~ dnorm(log_nu_mu, tau = log_nu_tau)  # mean litter size
  log(nu) <- log_nu

  ## convert to expected number of pigs per primary period
  log_zeta <- log(pp_len) + log_nu - log(365)
  log(zeta) <- log_zeta
  for(i in 1:n_ls){
    J[i] ~ dpois(nu)
  }

  for(i in 1:n_survey){

    log_potential_area[i] <- calc_log_potential_area(
      log_rho = log_rho[1:n_method],
      log_gamma = log_gamma[1:2],
      p_unique = p_unique[1:2],
      log_effort_per = log_effort_per[i],
      effort_per = effort_per[i],
      n_trap_m1 = n_trap_m1[i],
      log_pi = log_pi,
      method = method[i]
    )

    # base probability of capture,
    if(base){
      logit(theta[i]) <- beta1[method[i]] +
        inprod(X_p[i, 1:m_p], beta_p[method[i], 1:m_p])
    }

    # base + project random effect
    if(include_project){
      logit(theta[i]) <- beta1[method[i]] +
        inprod(X_p[i, 1:m_p], beta_p[method[i], 1:m_p]) +
        alpha_project[project[i]]
    }

    # base + cluster random effect
    if(include_cluster){
      logit(theta[i]) <- beta1[method[i]] +
        inprod(X_p[i, 1:m_p], beta_p[method[i], 1:m_p]) +
        alpha_cluster[cluster[i]]
    }

    # base + project + cluster random effects
    if(include_project_cluster){
      logit(theta[i]) <- beta1[method[i]] +
        inprod(X_p[i, 1:m_p], beta_p[method[i], 1:m_p]) +
        alpha_project[project[i]] +
        alpha_cluster[cluster[i]]
    }

    # capture given that an individual is in the surveyed area
    log_theta[i] <- log(theta[i]) + min(0, log_potential_area[i] - log_survey_area_km2[i])

    # likelihood
    y[i] ~ dpois(p[i] * (N[nH_p[i]] - y_sum[i]))

  }

  # the probability an individual is captured on the first survey
  for(i in 1:n_first_survey){
    log(p[first_survey[i]]) <- log_theta[first_survey[i]]
  }

  # the probability an individual is captured after the first survey
  for(i in 1:n_not_first_survey){
    log(p[not_first_survey[i]]) <- log_theta[start[not_first_survey[i]]] +
      sum(log(1 - exp(log_theta[start[not_first_survey[i]]:end[not_first_survey[i]]])))
  }

  for(i in 1:n_clusters){

    log_lambda_1[i] ~ dunif(0, 10)
    log(M[mH[i, 1]]) <- log_lambda_1[i]

    # population growth across time steps
    for(j in 2:n_time_clust[i]){ # loop through every PP, including missing ones

      lambda[mH[i, j-1]] <- (M[mH[i, j-1]] - rem[i, j-1]) * zeta / 2 +
        (M[mH[i, j-1]] - rem[i, j-1]) * phi[mH[i, j-1]]

      M[mH[i, j]] ~ dpois(lambda[mH[i, j-1]])
      phi[mH[i, j-1]] ~ dbeta(a_phi, b_phi)

    }

  }

  # single property clusters
  # property abundance = cluster abundance
  for(i in 1:n_single_property_clusters){
    for(t in 1:n_time_single_property_clusters[i]){
      N[nH_single[i, t]] <- M[nmH_single[i, t]]
    }
  }

  # multiple property clusters
  # distribute based on average density
  for(i in 1:n_multi_property_clusters){
    for(t in 1:n_time_multi_property_clusters[i]){
      N[nH_multi[i, t]] ~ dpois(d[i, t])
      log(d[i, t]) <- log(M[nmH_multi[i, t]]) - log_cluster_area[i] + log_property_area[i]
    }
  }

})
