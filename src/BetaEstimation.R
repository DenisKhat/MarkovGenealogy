source('src/MarkhovChain.R')


methodA <- function(I, t, N){
  LA <- function(beta) markhov_probability(t, beta, 0, N-1, 1)[[I]]
  return(optimize(LA, interval=c(0,1), maximum = TRUE)[["maximum"]])
}


table <- markhov_virus(60, 1, 0, 59, 1)
combined <- list(times = table$time, I=table$I)
df = as.data.frame(do.call(cbind, combined))


methodB <- function(data, N){
  data <- data[order(do.call(rbind, data$times)),]
  row.names(data) <- NULL
  print(data)
  LB <- function(beta){
    val = markhov_probability(data[[1,1]], beta, 0, N-1, 1)[ data[[1,2]] ]
    for (i in 2:nrow(data)){
      val = val * markhov_probability(data[[i,1]] - data[[i-1,1]], beta, 0, N-data[[i-1,2]], data[[i-1,2]])[data[[i,2]]]
    }
    return(val)
  }
  return(optimize(LB, interval=c(0,2), maximum = TRUE)$maximum)
}

df = df[sample(nrow(df), 30), ]
# data <- data[order(do.call(rbind, data$times)),]
# print(df)
# df[[4,1]]
methodB(df, 60)

#df$times[i]