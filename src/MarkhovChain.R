library(deSolve)
library(ggplot2)
library(gridExtra)  # install.packages("gridExtra")
library(Matrix)


markhov_virus <- function(end_time, beta, gamma, S, I, R = 0){
  #output value as: "time_of_event: infected/0, infected/recovered, 
  N <- S + I + R
  population <- 1:N
  p_0 = sample(population, 1)
  I_list <- list(p_0)
  S_list <- population[-p_0]
  R_list <- c()
  table <- matrix(nrow=0, ncol=7)
  colnames(table) <- c("time", "infector", "affected", "S", "I", "R", "I List")
  current_time <- 0
  infecter <- 0
  affected <- 0
  table <- rbind(table, c(0,0,0,S,I,R, toString(I_list)))
  while (current_time < end_time & I > 0 & (gamma > 0 || I < N)){
    next_time <- rexp(1, beta * S *I / N + gamma *I)
    current_time <- current_time + next_time
    U <- runif(1)
    if (U < (beta * S / N) / (beta * S / N + gamma)) {
      I = I + 1
      S = S - 1
      infecter <- sample(I_list, 1)
      affected <- sample(S_list, 1)
      S_list <- setdiff(S_list, c(affected))
      I_list <- c(I_list, affected)
    }
    else {
      # S = S + 1
      I = I - 1
      R = R + 1
      infecter <- sample(I_list, 1)
      affected <- 0
      I_list <- setdiff(I_list, c(infecter))
      # S_list <- c(S_list, affected)
      R_list <- c(R_list, infecter)
    }
    table <- rbind(table, c(current_time, infecter, affected, S, I, R, toString(I_list)))
  }
  table <- table[-nrow(table), ]
  table <- data.frame(table)
  table$time <- as.numeric(table$time)
  return(table)
}

table <- markhov_virus(end_time=10,beta=0.7, gamma=0, S=59, I=1)
p_0 <- table$infector[2]

markhov_probability <- function(times, beta, gamma, initial_S, initial_I, initial_R=0){
  # Returns a matrix of probabilities, given time and I.
  # output in form of matrix[time, I_value]
  #temporarily for this function, assume recovery back into susceptible, may change that later.
  N <- initial_S + initial_I + initial_R
  # Code for diff. eqtn solving
  lambda <- sapply(seq(0, N), function(i) beta * (N-i) * i / N)
  mu <- sapply(seq(0, N), function(i) gamma * i)
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


