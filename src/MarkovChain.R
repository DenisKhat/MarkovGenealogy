library(deSolve)
# library(ggplot2)
# library(gridExtra)  # install.packages("gridExtra")
library(Matrix)
source("Volz.R")


markov_virus <- function(end_time, beta, gamma, S, I, R = 0, curr_time=0, S_list=NULL, I_list=NULL, R_list=NULL, SIS=FALSE){
  #output value as: "time_of_event: infected/0, infected/recovered,
  # when reading the table, an infector of 0 implies a recovery event.
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

# table <- markhov_virus(end_time=10,beta=0.1, gamma=0, S=59, I=1, curr_time=5)
# print(table[, 1:6])

markov_probability_SIS <- function(times, beta, gamma, initial_S, initial_I){
  # Returns a matrix of probabilities, given some times, for all I.
  # output in form of matrix[time, I_value]
  N <- initial_S + initial_I
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
  return(mass_vectors)
}


markov_probability_SIR <- function(times, beta, gamma, initial_S, initial_I, initial_R=0){
  N = initial_I + initial_S + initial_R
  MATRIX_SIZE = (N+1) * (N+1)
  ### some functions to help create the matrix eventually
  get_SI_from_index <- function(index){
    return(c(S=floor((index-1)/(N+1)), I=((index - 1) %% (N+1))))
  }
  
  get_index_from_SI <- function(SI){
    return(SI[1]*(N+1) + SI[2] + 1)
}
  
  get_roi_from_SI <- function(SI){  
    # rate of infection, lambda. takes in SI from the get_SI_from_index function.
    return(beta*SI[1]*SI[2]/N)
  }
  
  get_ror_from_SI <- function(SI){ 
    # rate of recovery, mu. takes in SI from the get_SI_from_index function.
    return(gamma*SI[2])
  }
  ###
  
  initial_P <- numeric(MATRIX_SIZE)
  initial_P[get_index_from_SI(c(S=initial_S, I=initial_I))] = 1
  # print(initial_P)
  
  ### construct the matrix
  A <- matrix(data=0,nrow=MATRIX_SIZE, ncol=MATRIX_SIZE)
  for (j in seq(1,MATRIX_SIZE)){
   SI <- get_SI_from_index(j)
   A[j, j] <- - get_ror_from_SI(SI) - get_roi_from_SI(SI)
   if (SI[2] > 1 && SI[1] < N){
     k <- get_index_from_SI(SI + c(S=1, I=-1))
     A[j, k] <- get_roi_from_SI(SI+c(S=1, I=-1))
   }
   if (SI[2] < N){
     k <- get_index_from_SI(SI + c(S=0, I=1))
     A[j, k] <- get_ror_from_SI(SI + c(S=0, I=1))
   }
  }
  
  mass_vectors <- sapply(times, function(t) expm(A*t) %*% initial_P)
  return(mass_vectors)
}


expected_I <- function(beta, gamma, time, initial_S = 59, initial_I = 1){
  probabilities <- markov_probability_SIS(time, beta, gamma, initial_S, initial_I)
  N <- initial_I + initial_S
  return(sum(sapply(1:N, function(i) i * probabilities[i+1] )))
}

expected_I(0.7, 0.3, 10)


get_indices_for_I <- function(I, N){
  return(sapply(seq(0,N), function(s) s*(N+1) + I + 1))
}


get_event_likelihood_SI <- function(data, N=60, initial_I=1, end_time=10){
  # Data is a list of times. Log likelihood function input of form (beta).
  # Using pop size of 60.
  M <- length(data)
  L <- function(parameter){
    lambda = function(i) parameter * (N-i) * i / N
    if (M >= 60) out <- log(lambda(1)) - data[i] * lambda(1)
    else if (M == 0) return (end_time * lambda(1))
    else out <- (end_time - tail(data,1)) * lambda(M + 1) + log(lambda(1)) - data[i] * lambda(1)
    for (i in 2:M){
      out <- out + log(lambda(i)) - (data[i] - data[i-1])*lambda(i)
    }
    return(out)
  }
}


get_sample_times_likelihood_SIS <- function(data){
  # Data is a df of of the form (time,I), remember to fix times when using this.
  # representing the observed data. Returns a log likelihood function which can
  # then be optimized. log likelihood function input of form (beta, gamma).
  L <- function(parameter){
    probabilities <- markov_probability_SIS(data$time, parameter[1], parameter[2], 59, 1)
    likelihood <- 0
    for (i in length(probabilities)){
      likelihood <- likelihood + log((probabilities[[i]])[1,data$I[[i]]+1])
    }
    return(likelihood)
  }
  return(L)
}


get_exact_likelihood_SIS <- function(data, end_time=10){
  # data is a dataframe from markov_virus.
  # make sure data includes the initial state, i.e. time = 0
  # time is a real number, type is +1 if infection, -1 if recovery.
  N <- data$I[[1]] + data$S[[1]]
  initial_I <- data$I[[1]]
  M <- nrow(data)
  L <- function(parameter){ # parameter is c(beta,gamma)
    I <- initial_I
    lambda <- function(i) parameter[[1]] * (N-i) * i / N
    mu <- function(i) parameter[[2]] * i
    combined_rate <- function(i) lambda(i) + mu(i)
    out <- 0
    for(k in 2:M){
      out <- out + log(combined_rate(I)) - combined_rate(I) * (data$time[[k]] - data$time[[k-1]])
      if(data$affected[[k]] == 0){
        out <- out + log(mu(I)/combined_rate(I))
        I <- I - 1
      }
      else{
        out <- out + log(lambda(I)/combined_rate(I)) 
        I <- I + 1
      }
      # print(out)
      # print(I)
    }
    if(I != 0){
      # print("using thingy")
      # print(out)
      out <- out - combined_rate(I)*(end_time - data$time[[M]])
    }
    return(out)  
  }
  return(L)
}

# sample_dta <- markov_virus(end_time = 10,beta = 0.7,gamma=0.3,S = 59,I=1,R=0,SIS = T)
# sample_dta$I
# L <- get_exact_likelihood_SIS(sample_dta)
# print(L(c(0.7,0.3)))
# print(L(c(0.7,0.0)))

get_volz_likelihood_SIS <- function(data, N=60, initial_I=1){
  return(0)
}