source('src/MarkhovChain.R')

time = 10
N = 60

table <- markhov_virus(60, 0.7, 0, 59, 1)
combined <- list(times = table$time, I=table$I)
df = as.data.frame(do.call(cbind, combined))


I_at_t <- function(data, time){
  index <- which(table$time == time)  # get index of desired time
  if (length(index) == 0){ # if no event at exact time
    index <- which(table$time == closest_lower(time, table$time)) # look at closest earlier event
  }
  # get the number of infected at time
  I <- as.numeric(table$I[index])
  return(I)
}

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
      geom_line() +
      xlab('B') +
      ylab('L(B)')
  print(p)

  return(optimize(LA, interval=c(0,1), maximum = TRUE)[["maximum"]])
}
methodA(I, time, N)


# methodB <- function(data, N){
#   data <- data[order(do.call(rbind, data$times)),]
#   row.names(data) <- NULL
#   print(data)
#   LB <- function(beta){
#     val = markhov_probability(data[[1,1]], beta, 0, N-1, 1)[ data[[1,2]] ]
#     for (i in 2:nrow(data)){
#       val = val * markhov_probability(data[[i,1]] - data[[i-1,1]], beta, 0, N-data[[i-1,2]], data[[i-1,2]])[data[[i,2]]]
#     }
#     return(val)
#   }
#   return(optimize(LB, interval=c(0,2), maximum = TRUE)$maximum)
# }
# 
# df = df[sample(nrow(df), 30), ]
# # data <- data[order(do.call(rbind, data$times)),]
# # print(df)
# # df[[4,1]]
# methodB(df, 60)

#df$times[i]