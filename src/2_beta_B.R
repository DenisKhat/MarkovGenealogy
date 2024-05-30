source('src/MarkhovChain.R')
source('src/Experiment.R')
source('src/Trees_v2.R')
source('src/PiecewiseBeta.R')

N = 60
b1 = 0.6
b2 = 0.3
gamma = 0
change_time = 5
end_time = 10

first <- markhov_virus(change_time, b1, gamma, 59, 1)
# print(first[, 1:6])
first_R <- as.numeric(first$R[nrow(first)])
first_S <- as.numeric(first$S[nrow(first)])
first_I <- as.numeric(first$I[nrow(first)])
# print(first$S_list)


last <- markhov_virus(end_time, b2, gamma, first_S, first_I, first_R, curr_time=change_time, S_list=first$S_list[nrow(first)], I_list = first$I_list[nrow(first)], R_list=first$R_list[nrow(first)])

# first <- first[-1,]
last <- last[-1,] # to include/not include the 5
# print(last[, 1:6])

# PRINT GRAPH
full_table <- rbind(first, last)

I_at_t <- function(data, time){
  index <- which(table$time == time)  # get index of desired time
  if (length(index) == 0){ # if no event at exact time
    index <- which(table$time == closest_lower(time, table$time)) # look at closest earlier event
  }
  I <- as.numeric(table$I[index])
  return(I)
}

# METHOD B
loglike_B <- function(betas, data, T1 = 5, TF = 10, N = 60){
  b1 <- betas[1]
  b2 <- betas[2]

  first_prob = log(markhov_probability(data$time[1], b1, 0, N-1, 1)[data$I[1] + 1])
  second_prob = log(markhov_probability(data$time[2] - data$time[1], b2, 0, N-data$I[1], data$I[1])[data$I[2] + 1])
  final_prob = first_prob + second_prob
  return(final_prob)
}


sampled_times <- c(5, 10)
times <- full_table$time
sampled_I <- c()

for (i in sampled_times){
  index <- which(full_table$time == i) # get index of desired time
  
  if (length(index) == 0){ # if no event at exact time
    index <- which(full_table$time == closest_lower(i, full_table$time)) # look at closest earlier event
  }
  # get the number of infected at time
  num_of_inf <- as.numeric(full_table$I[index])
  
  sampled_I <- c(sampled_I, num_of_inf)
}

dfb = data.frame(times=as.vector(sampled_times), I=as.vector(sampled_I))
print(dfb)

likelihood_B <- function(beta) -loglike_B(beta, dfb, T1, TF, N)
betas_hat = optim(c(0.5,0.5), likelihood_B, method = "L-BFGS-B", lower=c(0.01, 0.01), upper=c(0.99, 0.99))
betas_hat$par
# 

heatmap_B <- function(data, T1, TF, N){
  b1 <- seq(0.01, 1, by = 0.01)
  b2 <- seq(0.01, 1, by = 0.01)
  
  LF <- function(b1, b2){
    betas <- c(b1, b2)
    loglike_B(betas, data, T1, TF, N)
  } 
  matrix_result <- outer(b1, b2, Vectorize(LF))
  print(matrix_result)
  
  # NORMALIZE
  # min_val <- min(matrix_result)
  # max_val <- max(matrix_result)
  # normalized_result <- (matrix_result - min_val) / (max_val - min_val)
  
  data <- expand.grid(b1 = b1, b2 = b2)
  # data$Z <- as.vector(normalized_result)
  data$Z <- as.vector(matrix_result)
  
  p <- ggplot(data, aes(b1, b2, fill= Z)) +
    geom_tile()+
    scale_fill_gradient(low="blue", high="yellow", limits = c(-25, 0))+
    geom_point(aes(x = betas_hat$par[1], y = betas_hat$par[2]), color = "black", size = 3)
  print(p)
}


heatmap_B(dfb, T1, TF, N)

# get_loglike_prof_1B <- function(data, T_1 = 5, T_f = 10, N = 60, precision=0.05){
#   sampled_points = seq(0+precision, 1, by=precision)
#   # Returns a fitted approx curve of profile likelihood, on domain of 0 to 1. 
#   # Doing nested optim is bad, so this "pre bakes" the function.
#   # Yes, max_at_b1 is just the likelihood function, but using it straight crashes R.
#   max_at_b1 <- function(beta1){
#     LB <- function(beta2){-loglike_B(c(beta1, beta2), data, T_1, T_f, N)}
#     mle <- optim(0.5, LB, method="L-BFGS-B", lower=0.01, upper= 0.99)
#     return(-mle$value)}
#   sampled_maxes = sapply(sampled_points, max_at_b1)
#   return(approxfun(x=sampled_points, y=sampled_maxes, method="linear", rule=1))
# }
# 
# 
# get_loglike_prof_2B <- function(data, T_1 = 5, T_f = 10, N = 60, precision=0.05){
#   sampled_points = seq(0+precision, 1, by=precision)
#   # Returns a fitted approx curve of profile likelihood, on domain of 0 to 1. 
#   # Doing nested optim is bad, so this "pre bakes" the function.
#   # Yes, max_at_b1 is just the likelihood function, but using it straight crashes R.
#   max_at_b2 <- function(beta2){
#     LB <- function(beta1){-loglike_B(c(beta1, beta2), data, T_1, T_f, N)}
#     mle <- optim(0.5, LB, method="L-BFGS-B", lower=0.01, upper= 0.99)
#     return(-mle$value)}
#   sampled_maxes = sapply(sampled_points, max_at_b2)
#   return(approxfun(x=sampled_points, y=sampled_maxes, method="linear", rule=1))
# }
# 
# likelihood_P1 <- get_loglike_prof_1B(data = dfb, T_1=change_time, T_f=end_time, N=N, precision = 0.001)
# likelihood_P2 <- get_loglike_prof_2B(data = dfb, T_1=change_time, T_f=end_time, N=N, precision = 0.001)
# 
# #fix beta1, find max beta 2 (profile likelihood for beta 1)
# beta1_hat = optim(0.5, function(x) -likelihood_P1(x), method = "L-BFGS-B", lower=0.005, upper=0.99)
# beta1_hat$par
# 
# beta2_hat = optim(0.5, function(x) -likelihood_P2(x), method = "L-BFGS-B", lower=0.005, upper=0.99)
# beta2_hat$par
# 
# chi <- qchisq(p = 0.95, df = 1)/2
# shifted_L1 <- function(x) likelihood_P1(x) + beta1_hat$value + chi
# shifted_L2 <- function(x) likelihood_P2(x) + beta2_hat$value + chi
# 
# # confidence interval b1 (Wilks' estimate)
# upper_b1 <- 1
# lower_b1 <- 0
# if (shifted_L1(0.01) < 0){
#   lower_b1 <- uniroot(shifted_L1, lower = 0.001, upper=beta1_hat$par)$root }
# if (shifted_L1(1) <= 0){
#   upper_b1 = uniroot(shifted_L1, lower = beta1_hat$par, upper = 0.99)$root }
# 
# #confidence interval b2 (Wilks' estimate)
# upper_b2 <- 1
# lower_b2 <- 0
# if (shifted_L2(0.01) < 0){
#   lower_b2 = uniroot(shifted_L2, lower = 0.001, upper=beta2_hat$par)$root}
# if (shifted_L2(1) <= 0){
#   upper_b2 = uniroot(shifted_L2, lower = beta2_hat$par, upper = 0.99)$root
# }
# 
# 
# print(lower_b1)
# print(upper_b1)
# 
# print(lower_b2)
# print(upper_b2)