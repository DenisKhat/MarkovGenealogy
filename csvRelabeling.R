# If a simulation breaks, or id's are not distributed properly, this tool should relabel every run with a unique id.
# Simulation with a unique id.

CSV_PATH <- "2betasimulation.csv"

data <- read.csv(CSV_PATH)

new_id <- 1
new_ids <- c(1)
ids <- data$id

for(i in 2:nrow(data)) {
  if (ids[i] != ids[i-1]){
    new_id <- new_id + 1
  }
  new_ids <- c(new_ids, new_id)
}

out_data <- data.frame(id=new_ids, time=data$time ,I=data$I)
write.csv(out_data,file="2beta_fixed_sim.csv")