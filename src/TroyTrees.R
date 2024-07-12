library(deSolve)
library(Matrix)
library(mosaicCalc)

get_bridge_pmf <- function(times,beta, gamma, initial, final, N){
  # The function does not consider lineages yet. 
  # Initial and final are given in this form:
  # list(I=i, time=t)
  
  ### Calculate the Kolmogorov Forward Equation (consider initial state)
  
  lambda <- function(I) {
    if (0 < I & I < N) return(I*(N-I)*beta/N)
    else return(0) }
  mu <- function(I){
    if (0 < I & I <= N) return(I*gamma)
    else return(0) }
  
  initial_P <- sapply(0:N, function(x) 0 + (x == initial$I)) 
  
  A <- Matrix(data = 0, nrow = N+1, ncol = N+1)
  for (j in 1:(N+1)){
    # index is infected + 1
    A[j,j] <- -lambda(j-1) - mu(j-1)
    if (j <= N) A[j,j+1] <- lambda(j-1)
    if (j > 1) A[j,j-1] <- mu(j-1) }
  
  forward_P <- function(t){
    # returns a vector with the probability of each state as an entry, 
    # at some time t.
    return(initial_P %*% expm(A * (t-initial$time)))
  }
  
  ### Now do it backwards in time to form a stochastic bridge (consider final state)
  
  up <- function(t, I) {
    # index still is infected + 1
    if (0 <= I & I < N){
      P <- forward_P(t)
      return(mu(I+1) * P[I + 2] / P[I + 1])
    }
    else return(0)
  }
  
  down <- function(t, I) {
    # index still is infected + 1
    if (0 < I & I <= N){
      P <- forward_P(t)
      return(lambda(I-1) * P[I] / P[I+1])
    }
    else return(0)
  } 

  final_P <- sapply(0:N, function(x) 0 + (x == final$I))
  
  fundamental_matrix <- function(t) {
    out <- Matrix(data = 0, nrow = N+1, ncol = N+1)
    for (j in 1:(N+1)){
      # print(up(t, j-1))
      # print(down(t, j-1))
      out[j,j] <- integrate(Vectorize(function(x) - up(x, j-1) - down(x, j-1)), lower = final$time, upper = t)[[1]]
      if (j <= N) out[j, j+1] <- integrate(Vectorize(function(x) down(x, j-1)), lower = final$time, upper = t)[[1]]
      if (j > 1) out[j, j-1] <- integrate(Vectorize(function(x) up(x, j-1)), lower = final$time, upper = t)[[1]]
    }
    return(out)
  }
  
  return( sapply(times, function(t) final_P %*% expm(fundamental_matrix(t))) )
  
}

prob = get_bridge_pmf(c(9),0.7, 0.3, list(I=1, time=0), list(I=30, time=10), N=60)[[1]]
prob <- as.numeric(prob)
prob
sum(prob)