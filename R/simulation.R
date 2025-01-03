
library(dplyr)
library(readr)

params <- read_csv("../pigs-property/data/posterior_95CI_range_all.csv")

draw_value <- function(x){
  params |> filter(node == x) |> pull(med)
}

phi_mu <- draw_value("phi_mu")
psi_phi <- draw_value("psi_phi")
nu <- draw_value("nu")

a_phi <- phi_mu * psi_phi
b_phi <- (1 - phi_mu) * psi_phi
zeta <- 28 * nu / 365

n_clusters <- 10
n_prop_in_cluster <- round(runif(n_clusters, 0.5, 13))
n_time <- 40

# storage
M <- matrix(NA, n_clusters, n_time)
y <- array(NA, dim = c(n_clusters, max(n_prop_in_cluster), n_time))

# initial population (cluster) size
area <- 250 # km2
initial_density <- 5

M[, 1] <- rpois(n_clusters, area * initial_density)

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

# ecological process
process_model <- function(N, zeta, a_phi, b_phi){
  n <- length(N)
  phi <- rbeta(n, a_phi, b_phi)
  lambda <- (N / 2) * zeta + (N * phi)
  rpois(n, lambda)
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



