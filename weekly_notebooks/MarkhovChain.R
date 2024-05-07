library(deSolve)
library(ggplot2)
library(gridExtra)  # install.packages("gridExtra")

markhov_virus <- function(end_time, beta, gamma, S, I, R = 0){
  #output value as: "time_of_event: infected/0, infected/recovered, 
  N = S + I + R
  population = seq(1, N, by=1)
  I_list = population[1:I]
  S_list = population[I:S+I]
  R_list = population[S+I:N]
  table = matrix(nrow=0, ncol=6)
  colnames(table) <- c("time", "infector", "affected", "S", "I", "R")
  current_time = 0
  infecter = 0
  affected = 0
  table <- rbind(table, c(0,0,0,S,I,R))
  while (current_time < end_time & I > 0 & (gamma > 0 | I < N)){
    next_time <- rexp(1, beta * S *I / N + gamma *I)
    current_time <- current_time + next_time
    U <- runif(1)
    if (U < (beta * S / N) / (beta * S / N + gamma)) {
      I = I + 1
      S = S - 1
      infecter <- sample(I_list, 1)
      affected <- sample(S_list, 1)
      S_list <- S_list[S_list != affected]
      I_list <- c(I_list, affected)
    }
    else {
      # S = S + 1
      I = I - 1
      R = R + 1
      infecter <- 0
      affected <- sample(I_list, 1)
      I_list <- I_list[I_list != affected]
      # S_list <- c(S_list, affected)
      R_list <- c(R_list, affected)
    }
    table <- rbind(table, c(current_time, infecter, affected, S, I, R))
  }
  return(table)
}

data.frame(markhov_virus(end_time=60,beta=0.88, gamma=0.4, S=59, I=1))
