library(ggplot2)

#source("src/MarkhovChain.R") # copy path from MarkhovChain File


# Function to find the number closest to a certain number as long as it is lower
closest_lower <- function(target, numbers) {
  max(numbers[numbers < target])
}

experiment <- function(n, time, N, beta, gamma, S, I){
  ids <- c(0:N)
  id_counts <- rep(0, N+1)
  for (i in 1:n){
    table <- markhov_virus(end_time=time+1,beta=beta, gamma=gamma, S=S, I=I) # produce table
    indices <- which(table$time == time)  # get indice of desired time
    if (length(indices) == 0){ # if no exact event at exact time
      # print(closest_lower(time, table$time))
      indices <- which(table$time == closest_lower(time, table$time)) # look at closest earlier event
    }
    
    # print(paste("infected at time:", table$I.List[indices]))
    # print(table[-7])
    # print(indices)
    # infected_at_t <- c()  # set of those infected at desired time
    
    # infected <- table$I.List[indices] # look at infected at desired time
    # print(infected)
    # infected <- as.integer(unlist(strsplit(as.character(infected), ", "))) # split into list of ints
    # infected_at_t <- c(infected_at_t, infected) # append infected
    # print(infected_at_t)
    
    # print(to_incr)
    # print(table$I[indices])
    num_of_inf <- as.numeric(table$I[indices])
    # print(num_of_inf)
    id_counts[num_of_inf+1] <- id_counts[num_of_inf+1] + 1
    # print(id_counts)
    
    # print(infected_at_t)
    # print(id_probs)
  }
  
  id_probs <- id_counts / n
  
  
  combined <- list(ids = ids, prob = id_probs)
  df = as.data.frame(do.call(cbind, combined))
  return(df)
}






