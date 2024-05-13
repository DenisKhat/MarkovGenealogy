source('src/MarkhovChain.R')


methodA <- function(I, t, N){
  LA <- function(beta) markhov_probability(t, beta, 0, N-1, 1)[[I]]
  return(optimize(LA, interval=c(0,1), maximum = TRUE)[["maximum"]])
}


table <- markhov_virus(10, 0.6, 0, 59, 1)[]
combined <- list(times = table$time, I=table$I)
data = as.data.frame(do.call(cbind, combined))


methodB <- function(data, N){
  # data <- data[order(data$times),]
  LB <- function(beta) markhov_probability(data$times[[1]], beta, 0, N-1, 1)[[data$I[[1]]]]
  for (i in 2:nrow(data)){
    LB <- LB(beta) * function(beta) markhov_probability(data$times[[i]] - data$times[[i-1]], beta, 0, N-data$I[[i-1]], data$I[[i-1]])[[data$I[[i]]]]
  }
  return(optimize(LB, interval=c(0,1), maximum = TRUE))[["maximum"]]
}

# print(data)
methodB(data, 60)

#df$times[i]