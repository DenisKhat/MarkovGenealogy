library(yaml)     # Uncomment when local
# library(yaml, lib.loc = "../packages")  # Uncomment when on slurm

params = read_yaml("slurm_simulation/params.yaml")  # Uncomment when local
# params = read_yaml("/params.yaml")  # Uncomment when on slurm

N <- params$num
I <- params$I

markov_virus <- function(beta, I=params$I, t=0, end_t=params$T_f){ #SI model
  I_list <- c(I)
  t_list <- c(t)
  while (t < end_t & I < N){
    S = N - I
    next_time <- rexp(1, beta * S *I / N)
    t <- t + next_time
    I <- I + 1
    I_list <- c(I_list, I)
    t_list <- c(t_list, t)
  }
  table <- data.frame(time=t_list, I=I_list)
  if (t > end_t){
    table <- head(table, -1)
  }
  return(table)
}

markov_probability <- function(times, beta, gamma=0, N=params$num, initial_I=params$I){
  # Returns a matrix of probabilities, given time (past initial condition) and I.
  # output in form of matrix[time, I_value]
  
  # Code for diff. eqtn solving
  lambda <- sapply(seq(0, N), function(i) beta * (N-i) * i / N)
  mu <- sapply(seq(0, N), function(i) gamma * i)
  
  # Construct matrix to throw into exponent.
  A <- matrix(0, nrow=N+1, ncol=N+1)
  for (i in seq(N+1)){
    A[i,i] <- -lambda[i] - mu[i]
    if (i <= N) A[i,i+1] <- lambda[i]
    if (i > 1) A[i,i-1] <- mu[i]
  }
  
  initial_p = sapply(seq(0,N), function(i) 0 + (i == initial_I))
  mass_vectors <- sapply(times, function(t) initial_p %*% expm(A*t))
  mass_matrix <- do.call(rbind, mass_vectors)
  rownames(mass_matrix) <- times
  colnames(mass_matrix) <- seq(0,N)
  return(as.matrix(mass_matrix))
}