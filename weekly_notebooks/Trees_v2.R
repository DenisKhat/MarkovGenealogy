plot.new()

source("C:/Users/Simon/MarkhovGenealogy/weekly_notebooks/MarkhovChain.R")
table <- data.frame(markhov_virus(end_time=10,beta=0.7, gamma=0.3, S=59, I=1))


patient_0 = 1
end_time = 10
height = 1


# begin with patient 0 and infection time 0, recovery time == end time, 
phylog <- function(table, recent, time_of_infection, recovery_time, end_time, counter){
  indeces = rev(which(table$infector == recent))
  print("_____")
  print(indeces)
  print(paste("recent", recent))
  print(paste("infects", table$affected[indeces[1]]))
  
  segment_count <- 0  # Initialize segment count for this node
  
  # if recovery
  if (length(indeces) > 0 && table$affected[indeces[1]] == 0){ 
    recovery_time = table$time[indeces[1]]
    # print(parent_height)
    print("drawn if recovery")
    segments(time_of_infection, height, recovery_time, height, lwd=2)
    segments(time_of_infection, height, time_of_infection,  height + 0.015 * counter, lwd=2)
    text(recovery_time + 0.1, height, labels = recent, cex = 0.6)
    indeces <- indeces[-1] #remove recovery
    segment_count <- segment_count + 1
  }
  # if patient does not recover
  else {
    
    # print(parent_height)
    print("drawn")
    
    
    segments(time_of_infection, height, recovery_time , height, lwd=2)
    segments(time_of_infection, height, time_of_infection,  height + 0.015 * counter, lwd=2)
    text(recovery_time + 0.1, height, labels = recent, cex = 0.6)
    segment_count <- segment_count + 1
  }
  
  
  for (i in indeces){
    height <<- height - 0.015
    counter <- segment_count
    # print(paste("counter", counter))
    segment_count <- segment_count + phylog(table, table$affected[i], table$time[i], end_time, end_time, segment_count)
  }
  return(segment_count)
}

# Define tick positions and labels
ticks <- seq(0, end_time, by = 1)
labels <- as.character(ticks)
plot(1, type = "n", xlim = c(0, end_time + 0.1), ylim = c(0, 1), xlab = "Time", ylab = "", axes = FALSE)
axis(1, at = ticks, labels = labels)

phylog(table, 1, 0, end_time, end_time, 0)



print(table)