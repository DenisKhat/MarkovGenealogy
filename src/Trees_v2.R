plot.new()

source("src/MarkhovChain.R")  # copy path from MarkhovChain File
table <- markhov_virus(end_time=10,beta=0.7, gamma=0.5, S=59, I=1)
p_0 <- table$I_list[1]

end_time = 10


height = 1  # need as global variable as is tracked outside of recursion
# print(table)

# begin with patient 0 and infection time 0, recovery time == end time, 
phylog <- function(table, recent, end_time, time_of_infection = 0, counter = 0){
  # get indices of events involving recent id
  indices <- rev(which(as.numeric(unlist(table$infector)) == as.numeric(unlist(recent)))) 
  
  segment_count <- 0  # Initialize segment count for this id
  
  # if recovery
  if (length(indices) > 0 && table$affected[indices[1]] == 0){ 
    recovery_time <- table$time[indices[1]] # get recovery time
    # draw line
    segments(time_of_infection, height, recovery_time, height, lwd=2)
    segments(time_of_infection, height, time_of_infection,  height + 0.016 * counter, lwd=2)
    text(recovery_time + 0.1, height + 0.003, labels = recent, cex = 0.6)
    indices <- indices[-1] #remove recovery
    segment_count <- segment_count + 1
  }
  # if patient does not recover
  else {
    # draw lines till end time since no recovery
    segments(time_of_infection, height, end_time, height, lwd=2)
    segments(time_of_infection, height, time_of_infection,  height + 0.016 * counter, lwd=2)
    text(end_time + 0.1, height + 0.003, labels = recent, cex = 0.6)
    segment_count <- segment_count + 1
  }
  
  
  for (i in indices){
    height <<- height - 0.016 # decrement height for next line to be drawn
    segment_count <- segment_count + phylog(table, table$affected[i], end_time, table$time[i], segment_count)
  }
  return(segment_count)
}

# draw axis and set plot
ticks <- seq(0, end_time, by = 1)
labels <- as.character(ticks)
plot(1, type = "n", xlim = c(0, end_time + 0.1), ylim = c(0, 1), xlab = "Time", ylab = "", axes = FALSE)
axis(1, at = ticks, labels = labels)



# TO DISPLAY:
print(table[, 1:6])
phylog(table, p_0, end_time)
