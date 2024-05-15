source('src/MarkhovChain.R')
source('src/Experiment.R')

time = 10
N = 60
beta = 0.4

table <- markhov_virus(10, beta, 0, 59, 1)
combined <- list(times = table$time, I=table$I)
df = as.data.frame(do.call(cbind, combined))


I_at_t <- function(data, time){
  index <- which(table$time == time)  # get index of desired time
  if (length(index) == 0){ # if no event at exact time
    index <- which(table$time == closest_lower(time, table$time)) # look at closest earlier event
  }
  I <- as.numeric(table$I[index])
  return(I)
}

get_MLE <- function(func, alpha){
  MLE <- optim(0.5, func, method = "Brent", lower = 0, upper= 1, hessian = TRUE)
  
  mle <- MLE$par
  var <- solve(MLE$hessian[[1,1]])
  
  z <- qnorm(alpha / 2, lower.tail = FALSE)
  lower <- (mle - z * sqrt(var))
  upper <- (mle + z * sqrt(var))
    
  return(list(mle, lower, upper))
}

# methodA <- function(I, t, N){
#   x <- seq(0, 1, by = 0.01)
#   y <- c()
#   LA <- function(beta) markhov_probability(t, beta, 0, N-1, 1)[[I]]
#   
#   log_LA <- function(beta){
#     -(log(LA(beta)))
#   }
# 
#   for (b in x){
#     y <- c(y, -log_LA(b))
#   }
#   
#   mle <- get_MLE(log_LA, 0.05)
#   max <- mle[[1]]
#   lower <- mle[[2]]
#   upper <- mle[[3]]
#   
#   
#   y_max <- -log_LA(max)
#   y_lower <- -log_LA(lower)
#   y_upper <- -log_LA(upper)
#   y_true <- -log_LA(beta)
#   
#   
#   df <- data.frame(beta = x, LA = y)
#   p <- ggplot(df, aes(x = beta, y = LA)) +
#     geom_point() +
#     geom_line() +
#     annotate("segment", x = lower, y = -Inf, xend = lower, yend = y_lower, linetype = "dashed", color = "red") +
#     annotate("segment", x = upper, y = -Inf, xend = upper, yend = y_upper, linetype = "dashed", color = "red") +
#     annotate("segment", x = max, y = -Inf, xend = max, yend = y_max, color = "red") + 
#     annotate("segment", x = beta, y = -Inf, xend = beta, yend = y_true, color = "blue") + 
#     ggtitle("Method A") + 
#     xlab('Beta') +
#     ylab('Log-Likelihood')
#   # print(p)
# 
#   return(list(mle, p))
# }
# methodA(I, time, N)

methodB <- function(data, N){  # METHODS A AND B
  x <- seq(0, 1, by = 0.01)
  y <- c()
  
  data <- data[order(data$times),]
  row.names(data) <- NULL
  LB <- function(beta){
    val = markhov_probability(data[[1,1]], beta, 0, N-1, 1)[ data[[1,2]] + 1 ]
    if (length(data) > 2){
      for (i in 2:nrow(data)){
        val = val * markhov_probability(data[[i,1]] - data[[i-1,1]], beta, 0, N-data[[i-1,2]], data[[i-1,2]])[[data[[i,2]] + 1 ]]
      }
    }
    return(val)
  }
  
  log_LB <- function(beta){
    -(log(LB(beta)))
  }
  
  for (b in x){
    y <- c(y, -log_LB(b))
  }
  
  mle <- get_MLE(log_LB, 0.05)
  max <- mle[[1]]
  lower <- mle[[2]]
  upper <- mle[[3]]
  
  
  y_max <- -log_LB(max)
  y_lower <- -log_LB(lower)
  y_upper <- -log_LB(upper)
  y_true <- -log_LB(beta)
  
  
  df <- data.frame(beta = x, LA = y)
  p <- ggplot(df, aes(x = beta, y = LA)) +
    # geom_point() +
    geom_line() +
    annotate("segment", x = lower, y = -Inf, xend = lower, yend = y_lower, linetype = "dashed", color = "red") +
    annotate("segment", x = upper, y = -Inf, xend = upper, yend = y_upper, linetype = "dashed", color = "red") +
    annotate("segment", x = max, y = -Inf, xend = max, yend = y_max, color = "red") +
    annotate("segment", x = beta, y = -Inf, xend = beta, yend = y_true, color = "blue") + 
    ggtitle("Method A/B") + 
    xlab('Beta') +
    ylab('Log-Likelihood') +
    labs(caption = paste("lower:", lower, "max:", max, "upper:", upper, "beta:", beta))
  # print(p)
  
  return(list(mle, p))
}


methodC <- function(data, t_final,N){
  x <- seq(0, 1, by = 0.01)
  y <- c()
  #data must be incremented times. (corresponding to time of infections)
  s_obs = c(0,diff(data))
  M = length(data)
  LC <- function(beta) {
    lambda = sapply(seq(N), function(i) beta * (N-i) * i / N)
    pdfs = sapply(seq(M), function(i) dexp(s_obs[i],lambda[i]))
    p_M_last = 1 - pexp(t_final-data[M], lambda[M+1])
    return(p_M_last * prod(as.vector(pdfs)))
  }
  
  # print(y)
  df <- data.frame(beta = x, LA = y)
  p <- ggplot(df, aes(x = beta, y = LA)) +
      geom_point() +
      geom_line() +
      xlab('B') +
      ylab('L(B)')
  print(p)

  return(optimize(LA, interval=c(0,1), maximum = TRUE)[["maximum"]])
}
# methodA(I, time, N)


methodB <- function(data, N){
  x <- seq(0, 1, by = 0.01)
  y <- c()
  
  log_LC <- function(beta){
    -(log(LC(beta)))
  }
  
  for (b in x){
    y <- c(y, -log_LC(b))
  }
  
  mle <- get_MLE(log_LC, 0.05)
  max <- mle[[1]]
  lower <- mle[[2]]
  upper <- mle[[3]]
  
  
  y_max <- -log_LC(max)
  y_lower <- -log_LC(lower)
  y_upper <- -log_LC(upper)
  y_true <- -log_LC(beta)
  
  
  df <- data.frame(beta = x, LA = y)
  p <- ggplot(df, aes(x = beta, y = LA)) +
    # geom_point() + 
    geom_line() +
    annotate("segment", x = lower, y = -Inf, xend = lower, yend = y_lower, linetype = "dashed", color = "red") +
    annotate("segment", x = upper, y = -Inf, xend = upper, yend = y_upper, linetype = "dashed", color = "red") +
    annotate("segment", x = max, y = -Inf, xend = max, yend = y_max, color = "red") + 
    annotate("segment", x = beta, y = -Inf, xend = beta, yend = y_true, color = "blue") + 
    ggtitle("Method C") + 
    xlab('Beta') +
    ylab('Log-Likelihood') + 
    labs(caption = paste("lower:", lower, "max:", max, "upper:", upper, "beta:", beta))
  # print(p)
  
  
  return(list(mle, p))
}

# dfb = df[sample(nrow(df), 8),]


# just some code to fake sampling data from one simulation
I <- I_at_t(df, time)
sampled_times = c(10) # FOR METHOD A/B
times = as.vector(df[["times"]])
times = do.call(rbind, times)

sampled_I = sapply(sampled_times, function(i) df[[which(times == closest_lower(i,times)), 2]])

dfb = data.frame(times=as.vector(sampled_times), I=as.vector(sampled_I))
dfc = as.vector(do.call(rbind, df[["times"]]))

# A <- methodA(I, time, N)
B <- methodB(dfb, 60)
C <- methodC(dfc, 10, 60)

# A_plot <- A[[2]]
B_plot <- B[[2]]
C_plot <- C[[2]]

grid.arrange( B_plot, C_plot, nrow = 2)

