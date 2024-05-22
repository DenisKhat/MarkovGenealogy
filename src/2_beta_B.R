source('src/MarkhovChain.R')
source('src/Experiment.R')
source('src/Trees_v2.R')
source('src/PiecewiseBeta.R')

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
    scale_fill_gradient(low="blue", high="red")
  print(p)
}

b1 = 0.7
b2 = 0.3
T1= 5
TF = 10
N = 60

# heatmap_B(dfb, T1, TF, N)
likelihood_B <- function(beta) -loglike_B(beta, dfb, T1, TF, N)
betas_hat = optim(c(0.5,0.5), likelihood_B, method = "L-BFGS-B", lower=c(0.01, 0.01), upper=c(0.99, 0.99))
betas_hat$par