source('MarkovChain.R')
# library(deSolve, lib.loc = "../packages") #Uncomment on slurm
# library(Matrix, lib.loc = "../packages")
# library(yaml, lib.loc = "../packages")

# set up cluster
args <- commandArgs(TRUE)
# job_id <- as.numeric(args[1])
job_id <- 2
dir <- getwd()

num_of_sims <- params$nsim
ncores <- params$ncores  # will be from sim file

m <- ceiling(num_of_sims / ncores)

df <- data.frame()  

for (i in ((job_id-1)*m+1):min(job_id*m, num_of_sims)){
  data <- markov_virus(10, 0.7, 0.3, S=59, I=1,SIS=TRUE)
  data <- data[,c("time", "I")]
  sample_states <- data_frame() # states sampled at time 5, 10
  for(i in c(5.0,10.0)){
    sample_states <- rbind(sample_states, tail(data[which(data$time < i),],1))
  }
  LS <- get_sample_times_likelihood_SIS(data)
  MLE <- optim(c(0.5,0.5), function(theta) -LS(theta), lower = c(0,0), upper = c(1,1))
  df <- rbind(df, list(id=i, beta=MLE$par[[1]], gamma=MLE$par[[2]]))
}

path <- file.path(dir,"Results")
if (!file.exists(path)) {
  dir.create(path)
}

output_table <- data.frame(id, beta, gamma)
setwd(path)
write.csv(df,file=paste0("sampled_likelihood",job_id,".csv"),row.names=FALSE, col.names = FALSE)
setwd(dir)

