library(ggplot2)

source("C:/Users/Simon/MarkhovGenealogy/weekly_notebooks/MarkhovChain.R") # copy path from MarkhovChain File

time <- 2
N <- 1

ids <- c(1:60)
id_counts <- rep(0, 60)

# Function to find the number closest to a certain number as long as it is lower
closest_lower <- function(target, numbers) {
  max(numbers[numbers < target])
}

for (i in 1:N){
  result <- data.frame(markhov_virus(end_time=10,beta=0.7, gamma=0, S=59, I=1)) # produce table
  table <- result$table
  
  if (gamma > 0){ # if recovery only check times that start with desired time, ie 2.1, 2.3, etc
    indices <- which(substr(as.character(table$time), 1, 1) == time)  # get indices of desired time
    if (length(indices) == 0){ # if no events at time, use latest info, ie no events at time = 2 look at 1.9
      indices <- which(as.numeric(table$time) == closest_lower(time, as.numeric(table$time)))
    }
  }
  else{
    indices <- which(substr(as.character(table$time), 1, 1) <= time)  # get indices of desired time
  }
  print(table[-7])
  print(indices)
  infected_at_t <- c()  # set of those infected at desired time
  for (j in indices){
    infected <- as.integer(strsplit(table$I.List[j], ", ")[[1]])
    infected_at_t <- c(infected_at_t, infected)
    infected_at_t <- unique(infected_at_t)
  }
  to_incr <- match(infected_at_t, 1:60)
  print(to_incr)
  id_counts[to_incr] <- id_counts[to_incr] + 1
  print(id_counts)
  
  print(infected_at_t)
  print(id_probs)
}

id_probs <- id_counts / N


combined <- list(ids = ids, prob = id_probs)
df = as.data.frame(do.call(cbind, combined))

ggplot(df, aes(x = ids, y = id_probs)) + 
  geom_col() +
  labs(x = "Ids", y = "Probabilities") + 
  scale_y_continuous(breaks = seq(0, 1, by = 0.1)) + 
  scale_x_continuous(breaks = ids)




