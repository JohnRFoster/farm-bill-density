
library(dplyr)
library(readr)

source("R/functions_removal.R")

effort_data <- read_csv("../multi-method-removal-model/data/effort_data.csv") |>
  mutate(method_name = if_else(method_name == "TRAPS, CAGE", "TRAPS", method_name))

method_lookup <- effort_data |>
  select(method, method_name) |>
  mutate(idx = method) |>
  distinct()

start_density <- start_density
hierarchy <- "Poisson"

params <- read_csv("../pigs-property/data/posterior_95CI_range_all.csv")

draw_value <- function(x){
  params |>
    filter(grepl(x, node)) |>
    pull(med) |>
    round(2)
}

phi_mu <- draw_value("phi_mu")
psi_phi <- draw_value("psi_phi")
nu <- draw_value("nu")
beta_p <- matrix(draw_value("beta_p"), 5, 3, byrow = TRUE)
beta1 <- matrix(draw_value("beta1"), 5, 1)
beta <- cbind(beta1, beta_p)
omega <- draw_value("omega")
gamma <- draw_value("gamma")
rho <- draw_value("rho")

a_phi <- phi_mu * psi_phi
b_phi <- (1 - phi_mu) * psi_phi
zeta <- 28 * nu / 365


x <- 1
cluster_props <- spatial_clusters |>
  filter(cluster == x)

property_nums <- cluster_props$property

n_props <- nrow(cluster_props)
property_data <- properties[cluster_props$property]
survey_area <- list_c(sapply(property_data, function(x) x["area"]))

log_survey_area <- log(survey_area)
log_rho <- log(rho)
log_gamma <- log(gamma)
p_unique <- omega

sigma <- runif(1, 0, 5)

# storage
# initial population (cluster) size

# =======================
# TODO
# decide on what to do about properties with area larger than cluster size
# =======================
area <- sum(survey_area) # km2
Mspin <- rpois(1, area * start_density)

# ecological process
process_model <- function(N, zeta, a_phi, b_phi){
  n <- length(N)
  phi <- rbeta(n, a_phi, b_phi)
  lambda <- (N / 2) * zeta + (N * phi)
  rpois(n, lambda)
}

for(i in 1:6){
  Mspin <- process_model(Mspin, zeta, a_phi, b_phi)
}

# distribute based on density
d <- Mspin / area * survey_area

N <- matrix(NA, n_props, n_pp)

if(hierarchy == "Poisson"){
  N[,1] <- rpois(n_props, d)
} else if(hierarchy == "Gaussian"){
  alpha <- log(d)
  lambda <- rnorm(n_props, alpha, sigma)
  N[,1] <- round(exp(lambda))
}

M <- rep(NA, n_pp)
M[1] <- Mspin

take <- tibble()
for(t in 2:n_pp){

  # reset catch
  Y <- numeric(n_props)

  for(i in seq_len(n_props)){
    tmp_data <- property_data[[i]]

    removal_effort <- tmp_data$effort
    sample_occasions <- removal_effort$sample_occasions

    remove_pigs <- t %in% sample_occasions

    if(remove_pigs){
      sample_effort <- removal_effort |> filter(sample_occasions == t)

      # determine order of removal events
      removal_order <- determine_removal_order(sample_effort)

      # conduct removals
      X <- land_cover[cluster_props$county[i], ]

      nn <- N[i, t-1]
      take_t <- conduct_removals(nn, removal_order, effort_data, log_survey_area[i],
                                 X, beta, t,
                                 log_rho, log_gamma, p_unique, method_lookup) |>
        mutate(property = property_nums[i],
               cluster = cluster,
               county = cluster_props$county[i],
               N = nn,
               M = M[t-1])

      take <- bind_rows(take, take_t)
      Y[i] <- sum(take_t$take)

      # how many pigs are left?
      zi <- N[i, t-1] - Y[i]
      assertthat::assert_that(zi >= 0,
                              msg = paste("More pigs removed than are alive! Time =", t,
                                          "Property =", i))

    }

  }

  Z <- M[t-1] - sum(Y)

  assertthat::assert_that(Z >= 0,
                          msg = paste("More pigs removed from cluster than are alive! Time =", t))

  M[t] <- process_model(Z, zeta, a_phi, b_phi)

  # distribute based on density
  d <- M[t] / area * survey_area

  if(hierarchy == "Poisson"){
    N[,t] <- rpois(n_props, d)
  } else if(hierarchy == "Gaussian"){
    alpha <- log(d)
    lambda <- rnorm(n_props, alpha, sigma)
    N[,t] <- round(exp(lambda))
  }


}



M <- matrix(NA, n_clusters, n_time)
y <- array(NA, dim = c(n_clusters, max(n_prop_in_cluster), n_time))




# spin up

# =========== start with spin up

# initial catch
for(i in 1:n_clusters){
  for(j in 1:n_prop_in_cluster[i]){
    for(t in 1:n_time){
      g <- rbinom(1, 1, 0.3)
      lambda <- ifelse(g == 0, g, exp(rnorm(1, -1, 0.1)))
      y[i, j, t] <- rpois(1, lambda)
    }
  }
}



for(t in 2:n_time){
  C <- apply(y[,,t-1], 1, sum, na.rm = TRUE)
  # C <- 0
  Z <- pmax(0, M[, t-1] - C)
  M[, t] <- process_model(Z, zeta, a_phi, b_phi)
}
M

# mean population size for the cluster
beta_cluster <- 5
# beta_cluster <- rnorm(1, 0, 1 / sqrt(0.1))
lambda_cluster <- exp(beta_cluster)


# (1) mean at each site in a cluster
# precision across properties in each cluster modeled explicitly
sd_cluster <- runif(1, 0.001, 10)
b1 <- rnorm(n_prop_in_cluster, beta_cluster, sd_cluster)
Z1 <- round(exp(b1))  # population for each site before removals

# (2) mean at each site in a cluster
# precesion across properties in each cluster = variance of Poisson
b2 <- rpois(n_prop_in_cluster, lambda_cluster)
Z2 <- b2

rem <- 2

R1 <- Z1 - rem # population for each site after removals
R2 <- Z2 - rem # population for each site after removals





N_site <- process_model(R1, zeta, a_phi, b_phi)
N_site <- process_model(R2, zeta, a_phi, b_phi)

p <- 0.3
y_site <- rpois(length(N_site), N_site * p)






hist(rbeta(1e6, a_phi, b_phi), breaks = 50)



