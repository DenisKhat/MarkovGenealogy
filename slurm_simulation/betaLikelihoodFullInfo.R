source("slurm_simulation/virusSimulation.R")  # Uncomment locally
# source("virusSimulation") # Uncomment on slurm

job_id = 1  # Uncomment locally
# args = commandArgs(TRUE)  # Uncomment on slurm
# job_id = as.numeric(args[1])  # Uncomment on slurm

rate <- function(i) params$beta * (N-i) * i / N

