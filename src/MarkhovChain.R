library(deSolve)
library(ggplot2)
library(gridExtra)  # install.packages("gridExtra")
library(Matrix)


markhov_virus <- function(end_time, beta, gamma, S, I, R = 0, curr_time=0, S_list=NULL, I_list=NULL, R_list=NULL, SIS=FALSE){
  #output value as: "time_of_event: infected/0, infected/recovered,
  N <- S + I + R
  pop <- 1:N
  
  # initialize I_list, S_list, R_list if none given
  # store as lists so they can be stored in matrix
  if (is.null(I_list)){
    I_list <- list(sample(pop, 1))
  }
  else{
    I_list <- I_list
  }
  if (is.null(S_list)){
    S_list <- list(setdiff(1:N, I_list[[1]]))
  }
  else{
    S_list <- S_list
  }
  if (is.null(R_list)){
    R_list <- list(c())
  }
  else{
    R_list <- R_list
  }
  p_0 = I_list[1]  # get patient 0 
  
  # initialize matrix
  table <- matrix(nrow=0, ncol=9)
  colnames(table) <- c("time", "infector", "affected", "S", "I", "R", "S_list", "I_list", "R_list")
  
  # initialize other params
  current_time <- curr_time
  infecter <- 0
  affected <- 0
  
  R <- length(unlist(R_list))  # get number of recovered
  table <- rbind(table, c(curr_time, 0, 0, S ,I, R, S_list, I_list, R_list))
  
  while (current_time < end_time & I > 0 & (gamma > 0 || I < N)){
    # unlist to access as vector
    I_list <- unlist(I_list)
    S_list <- unlist(S_list)
    R_list <- unlist(R_list)
    
    # get next time
    next_time <- rexp(1, beta * S *I / N + gamma *I)
    current_time <- current_time + next_time
    U <- runif(1)
    if (U < (beta * S / N) / (beta * S / N + gamma)) {
      I = I + 1
      S = S - 1
      if (length(I_list) == 1){
        infecter <- I_list[1]
      }else{
        infecter <- sample(I_list, 1)
      }
      affected <- sample(S_list, 1)
      S_list <- setdiff(S_list, c(affected))
      I_list <- c(I_list, affected)
    }
    else {
      # S = S + 1
      if (length(I_list) == 1){
        infecter <- I_list[1]
      }else{
        infecter <- sample(I_list, 1)
      }
      affected <- 0
      if (!SIS) {
        R = R + 1
        R_list <- c(R_list, infecter)
        }
      else {
        S = S + 1
        S_list <- c(S_list, infecter)
      }
      I = I - 1
      I_list <- setdiff(I_list, c(infecter))
      # S_list <- c(S_list, affected)
    }
    # put in list to store in matrix
    I_list <- list(I_list)
    S_list <- list(S_list)
    R_list <- list(R_list)
    table <- rbind(table, c(current_time, infecter, affected, S, I, R, S_list, I_list, R_list))
  }
  table <- data.frame(table)
  if (current_time > end_time){
    table <- head(table, -1)
  }
  table$time <- as.numeric(table$time)
  return(table)
}
# 
# table <- markhov_virus(end_time=10,beta=0.1, gamma=0, S=59, I=1, curr_time=5)
# print(table[, 1:6])

markhov_probability <- function(times, beta, gamma, initial_S, initial_I, initial_R=0, immunity=TRUE){
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
    if (i > 1 & !immunity) A[i,i-1] <- mu[i]
    }
  initial_p = sapply(seq(0,N), function(i) 0 + (i == initial_I))
  mass_vectors <- sapply(times, function(t) initial_p %*% expm(A*t))
  mass_matrix <- do.call(rbind, mass_vectors)
  rownames(mass_matrix) <- times
  colnames(mass_matrix) <- seq(0,N)
  return(as.matrix(mass_matrix))
}

markhov_probability(0, beta = 0.4, gamma = 0, initial_S = 59, initial_I = 1)[[1]]
