rm(list=ls())

source('virusSimulation.R')
# source('slurm_simulation/virusSimulation.R')
library(deSolve, lib.loc = "../packages")
library(Matrix, lib.loc = "../packages")
library(yaml, lib.loc = "../packages")

N = params$num
# N = 60
beta1= 0.6
beta2 = 0.3

# set up cluster
args <- commandArgs(TRUE)
job_id <- as.numeric(args[1])
# job_id <- 1
dir <- getwd()

num_of_sims <- params$nsim
ncores <- params$ncores  # will be from sim file

# num_of_sims <- 1
# ncores <- 1
m <- ceiling(num_of_sims / ncores)

# helper functions
I_at_t <- function(data, time){
  index <- which(data$time == time)  # get index of desired time
  if (length(index) == 0){ # if no event at exact time
    index <- which(data$time == closest_lower(time, data$time)) # look at closest earlier event
  }
  I <- as.numeric(data$I[index])
  return(I)
}

closest_lower <- function(target, numbers) {
  return(max(numbers[numbers < target]))
}

# METHOD B for sample times: 5, 10
loglike_B <- function(betas, data, T1 = 5, TF = 10, N = 60){
  b1 <- betas[1]
  b2 <- betas[2]
  
  first_prob = log(markov_probability(data$time[1], b1, 0, N, 1)[data$I[1] + 1])
  second_prob = log(markov_probability(data$time[2] - data$time[1], b2, 0, N, data$I[1])[data$I[2] + 1])
  final_prob = first_prob + second_prob
  return(final_prob)
}

get_loglike_prof_1B <- function(data, T_1 = 5, T_f = 10, N = 60, precision=0.05){
  sampled_points = seq(0+precision, 1, by=precision)
  # Returns a fitted approx curve of profile likelihood, on domain of 0 to 1. 
  # Doing nested optim is bad, so this "pre bakes" the function.
  # Yes, max_at_b1 is just the likelihood function, but using it straight crashes R.
  max_at_b1 <- function(beta1){
    LB <- function(beta2){-loglike_B(c(beta1, beta2), data, T_1, T_f, N)}
    mle <- optim(0.5, LB, method="L-BFGS-B", lower=0.01, upper= 0.99)
    return(-mle$value)}
  sampled_maxes = sapply(sampled_points, max_at_b1)
  return(approxfun(x=sampled_points, y=sampled_maxes, method="linear", rule=1))
}


get_loglike_prof_2B <- function(data, T_1 = 5, T_f = 10, N = 60, precision=0.05){
  sampled_points = seq(0+precision, 1, by=precision)
  # Returns a fitted approx curve of profile likelihood, on domain of 0 to 1. 
  # Doing nested optim is bad, so this "pre bakes" the function.
  # Yes, max_at_b1 is just the likelihood function, but using it straight crashes R.
  max_at_b2 <- function(beta2){
    LB <- function(beta1){-loglike_B(c(beta1, beta2), data, T_1, T_f, N)}
    mle <- optim(0.5, LB, method="L-BFGS-B", lower=0.01, upper= 0.99)
    return(-mle$value)}
  sampled_maxes = sapply(sampled_points, max_at_b2)
  return(approxfun(x=sampled_points, y=sampled_maxes, method="linear", rule=1))
}

conf_interval <- function(betahat, data, T1, T2, N){
  chi <- qchisq(p = 0.95, df = 1)/2
  
  likelihood_P1 <- get_loglike_prof_1B(data, T1, T2, N=N, precision = 0.001)
  likelihood_P2 <- get_loglike_prof_2B(data, T1, T2, N=N, precision = 0.001)
  
  #fix beta1, find max beta 2 (profile likelihood for beta 1)
  beta1_hat = betahat[1]
  
  beta2_hat = betahat[2]
  
  shifted_L1 <- function(x) likelihood_P1(x) + beta1_hat + chi
  shifted_L2 <- function(x) likelihood_P2(x) + beta2_hat + chi
  
  b1_inCI = NA
  b2_inCI = NA
  b1_CI_width = NA
  b2_CI_width = NA
  
  if (shifted_L1(0.01) > 0 || shifted_L1(1) > 0){
    b1_inCI = NA
    b1_CI_width = NA
  }
  else{
    # confidence interval b1 (Wilks' estimate)
    upper_b1 <- 1
    lower_b1 <- 0
    lower_b1 <- uniroot(shifted_L1, lower = 0.001, upper=beta1_hat)$root 
    upper_b1 = uniroot(shifted_L1, lower = beta1_hat, upper = 0.99)$root
    
    if (beta1_hat > lower_b1  && beta1_hat < upper_b1){
      b1_inCI = 1
    } else{
      b1_inCI = 0
    }
    b1_CI_width = upper_b1 - lower_b1
  }
  
  if (shifted_L2(0.01) > 0 || shifted_L2(1) > 0){
    b2_inCI = NA
    b2_CI_width = NA
  }else{
    #confidence interval b2 (Wilks' estimate)
    upper_b2 <- 1
    lower_b2 <- 0
    lower_b2 = uniroot(shifted_L2, lower = 0.001, upper=beta2_hat)$root
    upper_b2 = uniroot(shifted_L2, lower = beta2_hat, upper = 0.99)$root
    
    if (beta2_hat > lower_b2  && beta2_hat < upper_b2){
      b2_inCI = 1
    } else{
      b2_inCI = 0
    }
    b2_CI_width = upper_b2 - lower_b2
  }
  
  return(list(b1_inCI, b1_CI_width, b2_inCI, b2_CI_width))
}


id <- c()
b1_estimate <- c()
b2_estimate <- c()
b1_inCI <- c()
b2_inCI <- c()
b1_CI_width <- c()
b2_CI_width <- c()
table <- read.csv("simulation2beta.csv")

for (i in ((job_id-1)*m+1):min(job_id*m, num_of_sims)){
  
  if (length(table$id == i) == 0){
    next
  }
  
  id <- c(id, i)
  
  sim_table <- table[table$id == i, c("time", "I")]
  sampled_times=c(5, 10)
  
  times <- sim_table$time
  sampled_I <- c()
  
  for (j in sampled_times){
    index <- which(sim_table$time == j) # get index of desired time
    
    if (length(index) == 0){ # if no event at exact time
      index <- which(sim_table$time == closest_lower(j, sim_table$time)) # look at closest earlier event
    }
    # get the number of infected at time
    num_of_inf <- as.numeric(sim_table$I[index])
    
    sampled_I <- c(sampled_I, num_of_inf)
  }
  #data frame fit for method b
  #times, I
  dfb = data.frame(time=as.vector(sampled_times), I=as.vector(sampled_I))
  print(dfb)
  
  likelihood_B <- function(beta) -loglike_B(beta, dfb, T1, TF, N)
  
  betas_hat = optim(c(0.5,0.5), likelihood_B, method = "L-BFGS-B", lower=c(0.01, 0.01), upper=c(0.99, 0.99))
  betas_hat <- betas_hat$par
  b1_estimate <- c(b1_estimate, betas_hat[1])
  b2_estimate <- c(b2_estimate, betas_hat[2])
  
  CI <- conf_interval(betas_hat, dfb, sampled_times[1], sampled_times[2], N)
  
  b1_inCI <- c(b1_inCI, CI[1])
  b1_CI_width <- c(b1_CI_width, CI[2])
  
  b2_inCI <- c(b2_inCI, CI[3])
  b2_CI_width <- c(b2_CI_width, CI[4])
}

path <- file.path(dir,"Results")
if (!file.exists(path)) {
  dir.create(path)
}

output_table <- data.frame(b1_estimate, b1_inCI, b1_CI_width, b2_estimate, b2_inCI, b2_CI_width)
setwd(path)
write.csv(output_table,file=paste0("2betaMethodB",job_id,".csv"),row.names=FALSE)
setwd(dir)




