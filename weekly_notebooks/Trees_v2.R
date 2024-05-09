plot.new()

source("C:/Users/Simon/MarkhovGenealogy/weekly_notebooks/MarkhovChain.R")  # copy path from MarkhovChain File
result <- markhov_virus(end_time=10,beta=0.7, gamma=0, S=59, I=1)
table <- result$table
p_0 <- result$p_0

end_time = 10


height = 1  # need as global variable as is tracked outside of recursion
# print(table)

# begin with patient 0 and infection time 0, recovery time == end time, 
phylog <- function(table, recent, end_time, time_of_infection = 0, counter = 0){
  indices <- rev(which(as.numeric(table$infector) == as.numeric(recent)))
  # print("_____")
  # print(indices)
  # print(paste("recent", recent))
  # print(paste("infects", table$affected[indices[1]]))
  
  segment_count <- 0  # Initialize segment count for this node
  
  # if recovery
  if (length(indices) > 0 && table$affected[indices[1]] == 0){ 
    recovery_time = table$time[indices[1]]
    # print(parent_height)
    # print("drawn if recovery")
    # print(time_of_infection)
    segments(time_of_infection, height, recovery_time, height, lwd=2)
    segments(time_of_infection, height, time_of_infection,  height + 0.016 * counter, lwd=2)
    text(recovery_time + 0.1, height + 0.003, labels = recent, cex = 0.6)
    indices <- indices[-1] #remove recovery
    segment_count <- segment_count + 1
  }
  # if patient does not recover
  else {
    
    # print(parent_height)
    # print("drawn")
    
    # print(time_of_infection)
    segments(time_of_infection, height, end_time, height, lwd=2)
    segments(time_of_infection, height, time_of_infection,  height + 0.016 * counter, lwd=2)
    text(end_time + 0.1, height + 0.003, labels = recent, cex = 0.6)
    segment_count <- segment_count + 1
  }
  
  
  for (i in indices){
    height <<- height - 0.016
    counter <- segment_count
    # print(paste("counter", counter))
    # print(table$time[i])
    segment_count <- segment_count + phylog(table, table$affected[i], end_time, table$time[i], segment_count)
  }
  return(segment_count)
}

# AXIS
ticks <- seq(0, end_time, by = 1)
labels <- as.character(ticks)
plot(1, type = "n", xlim = c(0, end_time + 0.1), ylim = c(0, 1), xlab = "Time", ylab = "", axes = FALSE)
axis(1, at = ticks, labels = labels)

phylog(table, p_0, end_time)



print(table[-7])