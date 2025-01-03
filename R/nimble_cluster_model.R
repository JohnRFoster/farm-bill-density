
# prior for mean cluster population
for(i in 1:n_cluster){
  for(t in 1:n_time){
    alpha[i, t] ~ dnorm(0, 0.01)
    M[i, t] <- exp(alpha[i, t])
  }
}

# property populations from cluster population
for(i in 1:n_property){
  for(t in 1:n_time){
    lambda[i, t] ~ dnorm(alpha[cluster[i], t], tau_cluster)
    N[i, t] <- exp(lambda[i, t])
  }
}


# population dynamics
for(i in 1:n_property){
  for(t in 2:n_time){

    lambda[i, t] ~ dnorm(alpha[cluster[i], t], tau_cluster)
    N[i, t-1] <- exp(lambda[i, t])

    Z[i, t] <- N[i, t-1] - rem[i, t-1]
    N[i, t] <- Z[i, t] * phi[i, t] + Z[i, t] / 2 * zeta
  }
}




