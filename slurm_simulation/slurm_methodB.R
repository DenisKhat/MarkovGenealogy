source('virusSimulation.R')
source('MarkovChain.R')
library(deSolve, lib.loc = "../packages")
library(Matrix, lib.loc = "../packages")
library(yaml, lib.loc = "../packages")

beta = 0.5

# set up cluster
args <- commandArgs(TRUE)
job_id <- as.numeric(args[1])
# job_id <- 2
dir <- getwd()

num_of_sims <- params$nsim
ncores <- params$ncores  # will be from sim file

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

get_MLE <- function(func, alpha){
  MLE <- optim(0.5, func, method = "L-BFGS-B", lower = 0.01, upper= 1, hessian = TRUE)
  
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
  inCI = 0
  width = NA
  #if error
  if (curve(0) > 0 || curve(10) > 0){
    inCI = NA
  }
  else{
    lower_wilks = uniroot(curve, lower = 0, upper=mle)$root
    upper_wilks = uniroot(curve, lower = mle, upper = 10)$root
    #if in CI
    if (mle[1] > lower_wilks && mle[1] < upper_wilks){
      inCI = 1
      width = upper_wilks - lower_wilks
    } else{
      inCI = 0
      width = upper_wilks - lower_wilks
    }
  }
  return(list(mle, inCI, width))
  
  # plot(curve)
  # roots <- uniroot.all(f=function(x) -func(x) - wilks_cutoff, interval = c(0,1), n=2)
  
  # lower = uniroot(curve, lower = 0, upper=mle)$root
  # upper = uniroot(curve, lower = mle, upper = 10)$root
}

# vars
# I <- params$I   # will be from sim file
# end_time <- 10
# gamma <- 0

methodB <- function(data, N){  # METHODS A AND B
  x <- seq(0, 10, by = 0.01)
  y <- c()
  
  data <- data[order(data$times),]
  row.names(data) <- NULL
  LB <- function(beta){
    val = log(markov_probability(data[[1,1]], beta, 0, N, 1)[ data[[1,2]] + 1 ])
    if (length(data) > 1){
      for (i in 2:nrow(data)){
        val = val + log(markov_probability(data[[i,1]] - data[[i-1,1]], beta, 0, N, data[[i-1,2]])[[data[[i,2]] + 1 ]])
      }
    }
    return(val)
  }
  
  for (b in x){
    y <- c(y, LB(b))
  }
  
  mle <- get_MLE(function(x) -LB(x), 0.05)
  
  return(mle)
}

# results <- data.frame()
# results <- results[((job_id-1)*m+1):min(job_id*m, num_of_sims),]
# results$estimate <- 0
# results$inCI <- 0
# row.names(results) <- NULL
# print(results)

id <- c()
estimate <- c()
inCI <- c()
width <- c()
table <- read.csv("simulation1beta.csv")


for (i in ((job_id-1)*m+1):min(job_id*m, num_of_sims)){
  if (length(table$id == i) == 0){
    next
  }
  
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
  dfb = data.frame(times=as.vector(sampled_times), I=as.vector(sampled_I))
  # print(dfb)
  mle <- methodB(dfb, N)
  mle <- as.numeric(mle)
  if (mle[1] == 0.01){
    mle[1] = 0
  }
  id <- c(id, i)
  estimate <- c(estimate, mle[1])
  inCI <- c(inCI, mle[2])
  width <- c(width, mle[3])
  
  
  
  # results$estimate[i] <- estimate
  # results$inCI[i] <- inCI
}

# print(results)
path <- file.path(dir,"Results")
if (!file.exists(path)) {
  dir.create(path)
}

output_table <- data.frame(id, estimate, inCI, width)
setwd(path)
write.csv(output_table,file=paste0("MethodB",job_id,".csv"),row.names=FALSE)
setwd(dir)

