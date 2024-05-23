library(yaml)

rm(list=ls())

source('src/MarkhovChain.R')
# library(deSolve, lib.loc = "/global/home/hpc5441/packages")
# library(Matrix, lib.loc = "/global/home/hpc5441/packages")

params <- read_yaml("slurm_simulation/params.yaml")
# source("slurm_simulation/simulation.R")  # source simulations

# set up cluster
args <- commandArgs(TRUE)
job_id <- as.numeric(args[1])
job_id <- 2
dir <- getwd()

num_of_sims <- 100
ncores <- params$ncores  # will be from sim file

# helper functions
I_at_t <- function(data, time){
  index <- which(table$time == time)  # get index of desired time
  if (length(index) == 0){ # if no event at exact time
    index <- which(table$time == closest_lower(time, table$time)) # look at closest earlier event
  }
  I <- as.numeric(table$I[index])
  return(I)
}

closest_lower <- function(target, numbers) {
  return(max(numbers[numbers < target]))
}

get_MLE <- function(func, alpha){
  MLE <- optim(0.5, func, method = "L-BFGS-B", lower = 0.2, upper= 0.8, hessian = TRUE)
  
  mle <- MLE$par
  var <- solve(MLE$hessian[[1,1]])
  # z <- qnorm(alpha / 2, lower.tail = FALSE)
  chi <- qchisq(p = 1-alpha, df = 1)
  # lower <- (mle - z * sqrt(var))
  # upper <- (mle + z * sqrt(var))
  wilks_cutoff <-  -func(mle) - chi/2
  x_points <- seq(0, 10, by=0.01)
  y_points <- sapply(x_points, function(x) -func(x) - wilks_cutoff)
  curve <- approxfun(x_points, y_points)
  # plot(curve)
  # roots <- uniroot.all(f=function(x) -func(x) - wilks_cutoff, interval = c(0,1), n=2)
  lower_wilks = uniroot(curve, lower = 0, upper=mle)$root
  upper_wilks = uniroot(curve, lower = mle, upper = 10)$root
  # lower = uniroot(curve, lower = 0, upper=mle)$root
  # upper = uniroot(curve, lower = mle, upper = 10)$root
  return(list(mle, lower_wilks, upper_wilks))
}

# vars
N <-params$Num  # will be from sim file
I <- params$I   # will be from sim file
end_time <- 10
beta <- 0.7
gamma <- 0

table <- markhov_virus(end_time, beta, gamma, 59, I, 0)
sampled_times=c(5, 10)

times <- table$time
sampled_I <- c()

for (i in sampled_times){
  index <- which(table$time == i) # get index of desired time
  
  if (length(index) == 0){ # if no event at exact time
    index <- which(table$time == closest_lower(i, table$time)) # look at closest earlier event
  }
  # get the number of infected at time
  num_of_inf <- as.numeric(table$I[index])
  
  sampled_I <- c(sampled_I, num_of_inf)
}
#data frame fit for method b
#times, I
dfb = data.frame(times=as.vector(sampled_times), I=as.vector(sampled_I))


methodB <- function(data, N){  # METHODS A AND B
  x <- seq(0, 10, by = 0.01)
  y <- c()
  
  data <- data[order(data$times),]
  row.names(data) <- NULL
  LB <- function(beta){
    val = log(markhov_probability(data[[1,1]], beta, 0, N-1, 1)[ data[[1,2]] + 1 ])
    if (length(data) > 2){
      for (i in 2:nrow(data)){
        val = val + log(markhov_probability(data[[i,1]] - data[[i-1,1]], beta, 0, N-data[[i-1,2]], data[[i-1,2]])[[data[[i,2]] + 1 ]])
      }
    }
    return(val)
  }
  
  for (b in x){
    y <- c(y, LB(b))
  }
  
  mle <- get_MLE(function(x) -LB(x), 0.05)
  max <- mle[[1]]
  lower <- mle[[2]]
  upper <- mle[[3]]
  
  return(mle)
}

results <- data.frame()
results <- results[((job_id-1)*num_of_sims+1):(job_id*num_of_sims),]
results$estimate <- 0
results$inCI <- 0
print(results)

for (i in 1:num_of_sims){
  mle <- methodB(dfb, N)
  estimate <- mle[1]
  results$estimate[i] <- estimate
  if (mle[2] < estimate < mle[3]){
    results$inCI[i] <- 1
  }
}
  
  
  