source("src/Volz.R")
source("slurm_simulation/virusSimulation.R")
# source("virusSimulation.R") # Uncomment on slurm
# source("Volz.R")

LIKELIHOOD_SIMULATION_NUM = 800

# original_simulated <- readRDS("SIS_virus_infections.RDS")
original_simulated <- readRDS("src/SIS_virus_infections.RDS")
original_simulated <- subset(original_simulated, sim_num == LIKELIHOOD_SIMULATION_NUM)
original_simulated <- subset(original_simulated, time < 10)
original_simulated <- subset(original_simulated, time > 0)



dir <- getwd()
job_id <- 99  # Uncomment locally
# args <- commandArgs(TRUE)  # Uncomment on slurm
# job_id <- as.numeric(args[1])  # Uncomment on slurm
m <- 10
# sim_numbers <- ((job_id-1)*m + 1):(job_id*m)

sim_beta <- (floor(job_id / 10) + 1) * 0.1
sim_gamma <- ((job_id %% 10) + 1) * 0.1
# print(sim_beta)
# print(sim_gamma)

simulations <- markov_virus_full(10, beta=sim_beta, gamma=sim_gamma, S = 59, I = 1, SIS=T)
simulations <- cbind(sim_num = 1, simulations)

for (i in 2:10000){
  sim <- markov_virus_full(10, beta=sim_beta, gamma=sim_gamma, S = 59, I = 1, SIS=T)
  sim <- cbind(sim_num = i, sim)
  simulations <- rbind(simulations, sim)
}


path = file.path(dir, "Simulation")
if (!file.exists(path)){
  dir.create(path)
}


saveRDS(simulations, file=paste("Simulation/simulations_",job_id,".RDS"))

times <- seq(0, 10, by = 0.05)

experimental_A <- data.frame()
for (i in 1:10000){
# for (i in 1:10){
  run_data <- subset(simulations, sim_num == i)
  # print(run_data)
  experimental_A <- rbind(experimental_A, A_from_data(times=times, sample_simulation = run_data, final_time = 10))
}

A <- colMeans(experimental_A,na.rm = T)
A <- as.vector(A)

A_dt <- A_dt_from_data(times=times, A_values=A)
print(A_dt)
A_dt <- approxfun(times, A_dt)

log_likelihood <- -log(A[1] + tail(A, n=1))
for (t in as.vector(original_simulated$time)){
  log_likelihood <- log(A_dt(t)) + log_likelihood
}
# data <- run_simulation(sim_numbers[1])
# for (i in sim_numbers[-1]){
#   data <- rbind(data, run_simulation(i))
# }


output <- data.frame(beta=sim_beta, gamma=sim_gamma, likelihood=log_likelihood)
write_csv(output, file=paste("Simulation/likelihood_",job_id,".csv"))

