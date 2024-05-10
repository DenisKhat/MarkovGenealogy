library(ggplot2)

source("C:/Users/Simon/MarkhovGenealogy/weekly_notebooks/MarkhovChain.R") # copy path from MarkhovChain File

time <- 5
N <- 100

ids <- c(1:60)
id_counts <- rep(0, 60)

# Function to find the number closest to a certain number as long as it is lower
closest_lower <- function(target, numbers) {
  max(numbers[numbers < target])
}

for (i in 1:N){
  table <- markhov_virus(end_time=10,beta=0.7, gamma=0, S=59, I=1) # produce table
  indices <- which(table$time == time)  # get indice of desired time
  if (length(indices) == 0){ # if no exact event at exact time
    indices <- which(table$time == closest_lower(time, table$time)) # look at closest earlier event
  }
  
  # print(paste("infected at time:", table$I.List[indices]))
  # print(table)
  # print(indices)
  infected_at_t <- c()  # set of those infected at desired time
  
  # infected <- table$I.List[indices] # look at infected at desired time
  # print(infected)
  # infected <- as.integer(unlist(strsplit(as.character(infected), ", "))) # split into list of ints
  # infected_at_t <- c(infected_at_t, infected) # append infected
  # print(infected_at_t)
  
  # print(to_incr)
  # print(table$I[indices])
  num_of_inf <- as.numeric(table$I[indices])
  id_counts[num_of_inf] <- id_counts[num_of_inf] + 1
  # print(id_counts)
  
  # print(infected_at_t)
  # print(id_probs)
}

id_probs <- id_counts / N


combined <- list(ids = ids, prob = id_probs)
df = as.data.frame(do.call(cbind, combined))

ggplot(df, aes(x = ids, y = id_probs)) + 
  geom_col() +
  labs(x = "Ids", y = "Probabilities") + 
  scale_y_continuous(breaks = seq(0, 1, by = 0.1)) + 
  scale_x_continuous(breaks = ids)




