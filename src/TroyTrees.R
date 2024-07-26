library(deSolve)
library(Matrix)
library(mosaicCalc)
library(expm)

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
    return(initial_P %*% expm(A * (t-initial$time), do.sparseMsg = F))
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
      out[j,j] <- integrate(Vectorize(function(x)  up(x, j-1) + down(x, j-1)), lower = final$time, upper = final$time-t,subdivisions = 2000, rel.tol = 0.01)[[1]]
      if (j > 1) out[j, j-1] <- integrate(Vectorize(function(x) - down(x, j-1)), lower = final$time, upper = final$time-t,subdivisions = 2000, rel.tol=0.01)[[1]]
      if (j <= N) out[j, j+1] <- integrate(Vectorize(function(x) - up(x, j-1)), lower = final$time, upper = final$time-t,subdivisions = 2000, rel.tol=0.01)[[1]]
    }
    return(out)
  }
  return( sapply(times, function(s) final_P %*% expm(fundamental_matrix(final$time - s), do.sparseMsg = F)) )
}


#prob = get_bridge_pmf(c(0.0001),0.7, 0.3, list(I=1, time=0), list(I=11, time=10), N=60)#[[1,2]]


get_indexing_table <- function(N){
  df <- data.frame(index=c(),state=c())
  k <- 1
  l <- 1
  for (i in 1:((N)*(N+1)/2)){
    # add state
    df <- rbind(df, list(index=i, I=k, L=l))
    # setup next state
    if (l < k){
      l <- l + 1
    }
    else{
      k <- k + 1
      l <- 1
    }
  }
  return(df)
}
# prob <- as.numeric(prob)
# prob
# sum(prob)
get_partial_info_bridge_pmf <- function(time, beta, gamma, initial, final, N){
  ### NOTE HARDCODED TIME STARTS AT 1, TODO CHANGE THIS!!
  # The function does not consider lineages yet. 
  # Initial and final are given in this form:
  # list(I=i, L=l, time=t)
  indices <- get_indexing_table(N)
  # Helper functions for index table
  get_index <- function(i, l) subset(indices, I == i & L == l)$index
  get_state <- function(index){
    wanted_row <- indices[index,]
    return(list(I=wanted_row$I, L=wanted_row$L))
  }
  # Going forwards:
  lambda <- function(I) {
    if (0 <= I & I < N) return(I*(N-I)*beta/N)
    else return(0) }
  mu <- function(I){ 
    if (0 < I & I <= N) return(I*gamma)
    else return(0) }
  
  initial_P <- sapply(0:N, function(x) 0 + (x == initial$I)) 
  A <- matrix(data=0, nrow=N+1, ncol=N+1)
  for(i in 0:N){ # rememeber, the index is infected + 1
    A[i+1,i+1] <- - lambda(i) - mu(i)
    if (i > 0) A[i+1, i] <- mu(i-1)
    if (i < N) A[i+1, i+2] <- lambda(i+1)
  }
  forward_P <- function(t) return(expAtv(A, initial_P, t=t)$eAtv) 
  #use eAtv as more efficient.
  
  # Going backwards:
  chi <- function(I, L) choose(L, 2) / choose(I, 2)
  up <- function(t, I) {
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
  
  STATES_COUNT <- (N)*(N+1)/2
  final_P <- numeric(length = STATES_COUNT)
  final_P[get_index(final$I, final$L)] <- 1
  B <- matrix(0, nrow = STATES_COUNT, ncol=STATES_COUNT)
  for (k in 1:STATES_COUNT){
    curr <- get_state(k)
#    if(k %% 100==1){
#      print(curr)
#    }
    B[k,k] <- integrate(Vectorize(function(s) - up(t=final$time-s,I=curr$I) - down(t=final$time-s,I=curr$I)), lower = 0, upper = final$time - time)$value
    if (curr$L < curr$I){
      B[k, get_index(curr$I-1, curr$L)] <- 
        integrate(Vectorize(function(s) up(t=final$time-s,I=curr$I-1)), lower = 0, upper = final$time - time)$value
    }
    if (curr$I < N){
      B[k, get_index(curr$I+1, curr$L)] <- 
        integrate(Vectorize(function(s) down(t=final$time-s, I=curr$I+1)*(1-chi(I=curr$I+1,L=curr$L))), lower=0, upper= final$time - time)$value
      B[k, get_index(curr$I+1, curr$L+1)] <- 
        integrate(Vectorize(function(s) down(t=final$time-s, I=curr$I+1)*chi(I=curr$I+1,L=curr$L+1)), lower=0, upper = final$time - time)$value
    }
  } 
  # print(B)
  output_vector <- expm(B) %*% final_P
  output_df <- data.frame(row.names = c("I", "L", "P"))
  for (k in 1:STATES_COUNT){
    curr <- get_state(k)
    infected <- curr$I
    lineages <- curr$L
    probability <- output_vector[[k,1]]
    output_df <- rbind(output_df, list(I=infected, L=lineages, P=probability))
  }
  return(output_df)
}

#dff <- get_partial_info_bridge_pmf(5, 0.7, 0.3, initial=list(I=1, L=1, time=0), final=list(I=5, L=3, time=10), 10)