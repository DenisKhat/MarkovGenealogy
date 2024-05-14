source('src/MarkhovChain.R')
source('src/Experiment.R')

time = 10
N = 60

table <- markhov_virus(10, 0.4, 0, 59, 1)
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
  
  for (b in x){
    y <- c(y, LC(b))
  }
  spline_fun <- splinefun(x, y)
  spline_interp <- data.frame(beta = x, LA = spline_fun(seq(0, 1, by = 0.01)))
  p <- ggplot(spline_interp, aes(x = beta, y = LA)) +
    geom_point() +
    geom_line() +
    xlab('B') +
    ylab('L(B)')
  print(p)
  
  return(optimize(LC, interval=c(0,1), maximum = TRUE)[["maximum"]])
}
dfc = as.vector(do.call(rbind, df[["times"]]))
methodC(dfc, 10, 60)

I <- I_at_t(df, time)

methodA <- function(I, t, N){
  x <- seq(0, 1, by = 0.01)
  y <- c()
  LA <- function(beta) markhov_probability(t, beta, 0, N-1, 1)[[I]]
  
  for (b in x){
    y <- c(y, LA(b))
  }
  # print(y)
  spline_fun <- splinefun(x, y)
  spline_interp <- data.frame(beta = x, LA = spline_fun(seq(0, 1, by = 0.01)))
  # list_A <- list(beta = x, LA = y)
  # df_A <- as.data.frame(do.call(cbind, list_A))
  p <- ggplot(spline_interp, aes(x = beta, y = LA)) +
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
  
  data <- data[order(data$times),]
  row.names(data) <- NULL
  # print(data)
  LB <- function(beta){
    val = markhov_probability(data[[1,1]], beta, 0, N-1, 1)[[ data[[1,2]] + 1 ]]
    for (i in 2:nrow(data)){
      val = val * markhov_probability(data[[i,1]] - data[[i-1,1]], beta, 0, N-data[[i-1,2]], data[[i-1,2]])[[data[[i,2]] + 1 ]]
    }
    return(val)
  }
  
  for (b in x){
    y <- c(y, LB(b))
  }
  spline_fun <- splinefun(x, y)
  spline_interp <- data.frame(beta = x, LA = spline_fun(seq(0, 1, by = 0.01)))
  p <- ggplot(spline_interp, aes(x = beta, y = LA)) +
    geom_point() +
    geom_line() +
    xlab('B') +
    ylab('L(B)')
  print(p)
  
  return(optimize(LB, interval=c(0,2), maximum = TRUE)$maximum)
}
# methodB(df, 60)

# dfb = df[sample(nrow(df), 8),]

# just some code to fake sampling data from one simulation
sampled_times = c(3, 5, 7, 9, 11)
times = as.vector(df[["times"]])
times = do.call(rbind, times)
# print(times)
# print(df)
sampled_I = sapply(sampled_times, function(i) df[[which(times == closest_lower(i,times)), 2]])
# print(sampled_I)
dfb = data.frame(times=as.vector(sampled_times), I=as.vector(sampled_I))
print(dfb$times)
methodB(dfb, 60)

# data <- data[order(do.call(rbind, data$times)),]
# print(df)
# df[[4,1]]

