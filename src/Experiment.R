library(ggplot2)

#source("src/MarkhovChain.R") # copy path from MarkhovChain File


# Function to find the number closest to a certain number as long as it is lower
closest_lower <- function(target, numbers) {
  return(max(numbers[numbers < target]))
}

experiment <- function(n, time, N, beta, gamma, S, I, SIS=FALSE){
  I_at_time <- c(0:N) # x axis
  counts <- rep(0, N+1)
  for (i in 1:n){
    table <- markhov_virus(end_time=time+1,beta=beta, gamma=gamma, S=S, I=I, R=0, SIS=SIS) # produce table
    index <- which(unlist(table$time == time))  # get index of desired time
    
    if (length(index) == 0){ # if no event at exact time
      index <- which(unlist(table$time == closest_lower(time, unlist(table$time)))) # look at closest earlier event
    }
    
    # get the number of infected at time
    num_of_inf <- as.numeric(table$I[index])
    # increment count of I at time 
    counts[num_of_inf+1] <- counts[num_of_inf+1] + 1
  }
  
  # get probabilities
  probs <- counts / n
  
  
  combined <- list(I_at_time = I_at_time, probs = probs)
  df = as.data.frame(do.call(cbind, combined))
  return(df)
}


# df <- experiment(10, 10, 60, 0.7, 0, 59, 1)
# print(df)



