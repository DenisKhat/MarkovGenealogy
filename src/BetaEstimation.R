library(rootSolve)
source('src/MarkhovChain.R')
source('src/Experiment.R')

time = 10
N = 60
beta = 0.6

table <- markhov_virus(10, beta, 0, 59, 1)
# table$time <- as.numeric(table$time)
# table$I <- as.numeric(table$I)
# 
# print(class(time))
# print(class(table$time))
# print(range(table$time))
# print(time)


combined <- list(times = table$time, I=table$I)
df = as.data.frame(do.call(cbind, combined))
# print(df)

I_at_t <- function(data, time){
  index <- which(table$time == time)  # get index of desired time
  if (length(index) == 0){ # if no event at exact time
    index <- which(table$time == closest_lower(time, table$time)) # look at closest earlier event
  }
  I <- as.numeric(table$I[index])
  return(I)
}

get_MLE <- function(func, alpha){
  MLE <- optim(0.5, func, method = "L-BFGS-B", lower = 0.2, upper= 0.8, hessian = TRUE)
  
  mle <- MLE$par
  var <- solve(MLE$hessian[[1,1]])
  z <- qnorm(alpha / 2, lower.tail = FALSE)
  
  chi <- qchisq(p = 1-alpha, df = 1)
  lower <- (mle - z * sqrt(var))
  upper <- (mle + z * sqrt(var))
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
  return(list(mle, lower, upper, lower_wilks, upper_wilks))
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
  w_lower <- mle[[4]]
  w_upper <- mle[[5]]
  
  
  y_max <- LB(max)
  y_lower <- LB(lower)
  y_upper <- LB(upper)
  y_true <- LB(beta)
  y_w_lower <- LB(w_lower)
  y_w_upper <- LB(w_upper)
  
  
  df <- data.frame(beta = x, LB = y)
  p <- ggplot(df, aes(x = beta, y = LB)) +
    # geom_point() +
    geom_line() +
    annotate("segment", x = w_lower, y = -Inf, xend = w_lower, yend = y_w_lower, linetype = "dashed", color = "olivedrab") +
    annotate("segment", x = w_upper, y = -Inf, xend = w_upper, yend = y_w_upper, linetype = "dashed", color = "olivedrab") +
    annotate("segment", x = lower, y = -Inf, xend = lower, yend = y_lower, linetype = "dashed", color = "red") +
    annotate("segment", x = upper, y = -Inf, xend = upper, yend = y_upper, linetype = "dashed", color = "red") +
    annotate("segment", x = max, y = -Inf, xend = max, yend = y_max, color = "red") +
    annotate("segment", x = beta, y = -Inf, xend = beta, yend = y_true, color = "blue") + 
    coord_cartesian(xlim=c(0.01, 1), ylim=c(-100, y_max)) +
    ggtitle("Method A/B") + 
    xlab('Beta') +
    ylab('Log-Likelihood') +
    # labs(caption = paste("lower:", lower, "max:", max, "upper:", upper, "beta:", beta))
    labs(caption = paste("lower:",lower, "w_lower:",w_lower,"upper:",upper,"w_upper:",w_upper,"estimated:", max, "true:", beta))
  # print(p)
  
  return(list(mle, p))
}


methodC <- function(data, t_final, N){
  x <- seq(0, 10, by = 0.01)
  y <- c()
  #data must be incremented times. (corresponding to time of infections)
  s_obs = c(0,diff(data))
  M = length(data)
  LC <- function(beta) {
    if (M == 0){
      return(log(1-pexp(t_final, lambda[1])))
    }
    lambda = sapply(seq(N), function(i) beta * (N-i) * i / N)
    pdfs = sapply(seq(M), function(i) dexp(s_obs[i],lambda[i]))
    log_pdfs = sapply(pdfs[1:M-1], log)
    p_M_last = log(1 - pexp(t_final-data[M], lambda[M+1]))
    # if (is.na(p_M_last)) p_M_last <- 0
    # print(p_M_last)
    # print(p_M_last)
    # print(p_M_last * prod(as.vector(pdfs)))
    return(p_M_last + sum(as.vector(log_pdfs)))
  }
  
  for (b in x){
    y <- c(y, LC(b))
  }
  print(LC(0.5))
  mle <- get_MLE(function(x) -LC(x), 0.05)
  max <- mle[[1]]
  lower <- mle[[2]]
  upper <- mle[[3]]
  w_lower <- mle[[4]]
  w_upper <- mle[[5]]
  
  y_max <- LC(max)
  y_lower <- LC(lower)
  y_upper <- LC(upper)
  y_true <- LC(beta)
  y_w_lower <- LC(w_lower)
  y_w_upper <- LC(w_upper)
  
  df <- data.frame(beta = x, LC = y)
  p <- ggplot(df, aes(x = beta, y = LC)) +
    # geom_point() + 
    geom_line() +
    annotate("segment", x = w_lower, y = -Inf, xend = w_lower, yend = y_w_lower, linetype = "dashed", color = "olivedrab") +
    annotate("segment", x = w_upper, y = -Inf, xend = w_upper, yend = y_w_upper, linetype = "dashed", color = "olivedrab") +
    annotate("segment", x = lower, y = -Inf, xend = lower, yend = y_lower, linetype = "dashed", color = "red") +
    annotate("segment", x = upper, y = -Inf, xend = upper, yend = y_upper, linetype = "dashed", color = "red") +
    annotate("segment", x = max, y = -Inf, xend = max, yend = y_max, color = "red") + 
    annotate("segment", x = beta, y = -Inf, xend = beta, yend = y_true, color = "blue") + 
    coord_cartesian(xlim=c(0.01, 1), ylim=c(-100, y_max)) +
    ggtitle("Method C") + 
    xlab('Beta') +
    ylab('Log-Likelihood') + 
    # labs(caption = paste("lower:", lower, "max:", max, "upper:", upper, "beta:", beta))
    # labs(caption = paste("estimated:", max, "true:", beta))
    labs(caption = paste("lower:",lower, "w_lower:",w_lower,"upper:",upper,"w_upper:",w_upper,"estimated:", max, "true:", beta))
  # print(p)
  
  
  return(list(mle, p))
}

# dfb = df[sample(nrow(df), 8),]


# just some code to fake sampling data from one simulation
I <- I_at_t(df, time)

sampled_times <- c(5, 10) # FOR METHOD A/B
times <- as.list(df[["times"]])
times <- do.call(rbind, times)

sampled_I = sapply(sampled_times, function(i) df[[which(times == closest_lower(i,times)), 2]])

dfb = data.frame(times=as.vector(sampled_times), I=as.vector(sampled_I))
# print(dfb)
dfc = as.vector(do.call(rbind, df[["times"]]))


# A <- methodA(I, time, N)
B <- methodB(dfb, 60)
C <- methodC(dfc, 10, 60)

# A_plot <- A[[2]]
B_plot <- B[[2]]
C_plot <- C[[2]]

grid.arrange( B_plot, C_plot, nrow = 2)

