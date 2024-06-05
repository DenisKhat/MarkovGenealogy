source("slurm_simulation/virusSimulation.R")  # Uncomment locally
# source("virusSimulation") # Uncomment on slurm

dir <- getwd()
job_id = 2  # Uncomment locally
data <- read.csv("simulation2beta.csv")

# args = commandArgs(TRUE)  # Uncomment on slurm
# job_id = as.numeric(args[1])  # Uncomment on slurm

m <- ceiling(params$nsim/params$ncores)
sim_numbers <- ((job_id-1)*m + 1):(min(job_id*m,params$nsim)) 
  
get_rate <- function(beta) return(function(i) beta * (params$num-i) * i / params$num)


loglike_C <- function(betas, times, T_1 = params$T_1, T_f = params$T_f, N = params$num){
  # times is a list of [1, 2.3, 5, 10, ...] times when infections occurred.
  # T_1 is the time at which beta changes
  # T_f is the time at which our sampling ends.
  rate1 <- get_rate(betas[1])
  rate2 <- get_rate(betas[2])
  # print(log(rate2(N-4)))
  # Condition that NOTHING happens at anytime
  # if (length(times) == 0){
  #   return( -rate1(1)*T_1 - rate2(1)*(T_f-T_1))
  # }
  
  times <- sort(times)
  
  times1 <- times[times <= T_1] # times while rate 1 active
  times2 <- times[times > T_1]  # times while rate 2 active
  out <- 0
  
  K = length(times1)
  M = length(times2) + K
  #K will be the last time a beta1 event happens (not K-1), a little different to board notation.
  #Similarly, with M. M'th is last event that happened, not M-1th. (M+1 first time after T_f)
  
  # Condition that NOTHING happens while beta1 is active
  if (K == 0){
    out <- -rate1(1)*T_1
  }
  #Something does happen while beta1 active
  else {
    out <- log(rate1(1)) - rate1(1) * (times1[1])
    if (K > 1){
      for (i in 2:K){
        out <- out + log(rate1(i)) - rate1(i) * (times[i]-times[i-1])
      }
    }
    if (K < N){
      out <- out - rate1(K+1) * (T_1 - times[K])
    }
    }
  # Condition that NOTHING happens while beta2 is active 
  if (M==K){
    out <- out - rate2(M+1)*(T_f-T_1)
  }
  # Condition that SOMETHING happens while beta2 active.
  else {
    out <- out + log(rate2(K+1)) - rate2(K+1)*(times[K+1]-T_1)
    # print(M)
    if(M > K+1){
      for (i in seq(K+2,M)){
        # print(i)
        out <- out + log(rate2(i)) - rate2(i)*(times[i]-times[i-1])
      }
    }
    if (M < N){
      out <- out - rate2(M+1)*(T_f - times[M])
    }
  }
  return(out)
}


get_loglike_prof_1 <- function(times, T_1 = params$T_1, T_f = params$T_f, N = params$num, precision=0.01){
  sampled_points = seq(0+precision, 1, by=precision)
  # Returns a fitted approx curve of profile likelihood, on domain of 0 to 1. 
  # Doing nested optim is bad, so this "pre bakes" the function.
  # Yes, max_at_b1 is just the likelihood function, but using it straight crashes R.
  max_at_b1 <- function(beta1){
    LC <- function(beta2){-loglike_C(c(beta1, beta2), times, T_1, T_f, N)}
    mle <- optim(0.5, LC, method="L-BFGS-B", lower=0.01, upper= 0.99)
    return(-mle$value)}
  sampled_maxes = sapply(sampled_points, max_at_b1)
  return(approxfun(x=sampled_points, y=sampled_maxes, method="linear", rule=1))
}


get_loglike_prof_2 <- function(times, T_1 = params$T_1, T_f = params$T_f, N = params$num, precision=0.01){
  sampled_points = seq(0+precision, 1, by=precision)
  # Returns a fitted approx curve of profile likelihood, on domain of 0 to 1. 
  # Doing nested optim is bad, so this "pre bakes" the function.
  # Yes, max_at_b2 is just the likelihood function, but using it straight crashes R.
  max_at_b2 <- function(beta2){
    LC <- function(beta1){-loglike_C(c(beta1, beta2), times, T_1, T_f, N)}
    mle <- optim(0.5, LC, method="L-BFGS-B", lower=0.01, upper= 0.99)
    return(-mle$value)}
  sampled_maxes = sapply(sampled_points, max_at_b2)
  return(approxfun(x=sampled_points, y=sampled_maxes, method="linear", rule=1))
}



full_knowledge_estimation <- function(data, alpha=0.05){
  LC <- function(betas) -loglike_C(betas, data)
  likelihood_P1 <- get_loglike_prof_1(times=data)
  likelihood_P2 <- get_loglike_prof_2(times=data)
  MLE <- optim(c(0.5,0.5), LC, method = "L-BFGS-B", lower = 0.01, upper= 0.99, hessian = FALSE)
  in_1 <- 1
  in_2 <- 1
  chi <- qchisq(p = 1-alpha, df = 1)
  wilks_cutoff <-  -MLE$val - chi/2
  shifted_L1 <- function(x) likelihood_P1(x) - wilks_cutoff
  shifted_L2 <- function(x) likelihood_P2(x) - wilks_cutoff
  if (shifted_L1(0.02) > 0 || shifted_L1(0.98) > 0){
    # confidence interval not found in range.
    width_1 <- NA
    in_1 <- NA
  }
  else{
    lower_wilks = uniroot(shifted_L1, lower = 0.01, upper=MLE$par[1])$root
    upper_wilks = uniroot(shifted_L1, lower = MLE$par[1], upper = 1)$root
    width_1 <- upper_wilks - lower_wilks
    if (params$beta1 < lower_wilks || params$beta1 > upper_wilks){
      in_1 <- 0
    }
  }
    
  if (shifted_L2(0.02) > 0 || shifted_L2(0.98) > 0){
    # confidence interval not found in range.
    width_2 <- NA
    in_2 <- NA
  }
  
  else{
    lower_wilks = uniroot(shifted_L2, lower = 0.01, upper=MLE$par[2])$root
    upper_wilks = uniroot(shifted_L2, lower = MLE$par[2], upper = 1)$root
    width_2 <- upper_wilks - lower_wilks
    if (params$beta2 < lower_wilks || params$beta2 > upper_wilks){
      in_2 <- 0
    }
  }
  
  return(list(estimate_1=MLE$par[1], in_1=in_1, width_1=width_1, estimate_2=MLE$par[2], in_2=in_2, width_2=width_2))
}


ids <- c()
output_table <- data.frame(estimate_1 <- c(), in_1 <- c(), width_1 <- c(), estimate_2 <- c(), in_2 <- c(), width_2 <- c())

for (i in sim_numbers){
  times <- data[which(data$id == i), "time"]
  ids <- c(ids, i)
  output_table <- rbind(output_table, full_knowledge_estimation(times))
}

output_table <- cbind(ids, output_table)
output_table

