source("slurm_simulation/virusSimulation.R")  # Uncomment locally
# source("virusSimulation") # Uncomment on slurm
dir <- getwd()
# job_id = 1  # Uncomment locally
# args = commandArgs(TRUE)  # Uncomment on slurm
# job_id = as.numeric(args[1])  # Uncomment on slurm
m <- ceiling(params$nsim/params$ncores)
sim_numbers <- ((job_id-1)*m + 1):(min(job_id*m,params$nsim))

run_simulation <- function(id){
  first <- markov_virus(beta = params$beta1, I = params$I, t = 0, end_t = params$T_1)
  new_I <- tail(first$I, n = 1)
  last <- markov_virus(beta = params$beta2, I = new_I, t = params$T_1, end_t = params$T_f)
  data <- rbind(first, last[-1,])
  data <- cbind(id, data)
  return(data)
}

data <- run_simulation(sim_numbers[1])
for (i in sim_numbers[-1]){
  data <- rbind(data, run_simulation(i))
}


path = file.path(dir, "Results")
if (!file.exists(path)){
  dir.create(path)
}


write.csv(data,file=paste0("Results/two_beta_simulation",job_id,".csv"),row.names=F)
