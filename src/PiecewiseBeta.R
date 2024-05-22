
rate <- function(beta, N) return( function(i) beta * (N-i) * i / N )


loglike_C <- function(betas, times, T_1 = 5, T_f = 10, N = 60){
  # times is a list of [1, 2.3, 5, 10, ...] times when infections occurred.
  # T_1 is the time at which beta changes
  # T_f is the time at which our sampling ends.
  rate1 <- rate(betas[1], N)
  rate2 <- rate(betas[2], N)
  # print(log(rate2(N-4)))
  # Condition that NOTHING happens at anytime
  # if (length(times) == 0){
  #   return( -rate1(1)*T_1 - rate2(1)*(T_f-T_1))
  # }
  
  times <- sort(times)
  
  times1 <- times[times <= T_1] # times while rate 1 active
  times2 <- times[times > T_1]  # times while rate 2 active
  out <- 0
  
  K = length(times1)
  M = length(times2) + K
  #K will be the last time a beta1 event happens (not K-1), a little different to board notation.
  #Similarly, with M. M'th is last event that happened, not M-1th. (M+1 first time after T_f)
  # print("TIMES 1")
  # print(times1)
  # print("TIMES 2")
  # print(times2)
  # Condition that NOTHING happens while beta1 is active
  if (K == 0){
    out <- -rate1(1)*T_1
  }
  #Something does happen while beta1 active
  else {
    out <- log(rate1(1)) - rate1(1) * (times1[1])
    if (K > 1){
      for (i in 2:K){
        out <- out + log(rate1(i)) - rate1(i) * (times[i]-times[i-1])
      }
    }
    if (K < N){
      out <- out - rate1(K+1) * (T_1 - times[K])
    }
    }
  
  # Condition that NOTHING happens while beta2 is active 
  if (M==K){
    out <- out - rate2(M+1)*(T_f-T_1)
  }
  # Condition that SOMETHING happens while beta2 active.
  else {
    out <- out + log(rate2(K+1)) - rate2(K+1)*(times[K+1]-T_1)
    # print(M)
    if(M > K+1){
      for (i in seq(K+2,M)){
        # print(i)
        out <- out + log(rate2(i)) - rate2(i)*(times[i]-times[i-1])
      }
    }
    if (M < N){
      out <- out - rate2(M+1)*(T_f - times[M])
    }
    return(out)
  }
}
# a <- matrix(0, 10, 10)

get_loglike_prof_1 <- function(times, T_1 = 5, T_f = 10, N = 60, precision=0.05){
  sampled_points = seq(0+precision, 1, by=precision)
  # Returns a fitted approx curve of profile likelihood, on domain of 0 to 1. 
  # Doing nested optim is bad, so this "pre bakes" the function.
  # Yes, max_at_b1 is just the likelihood function, but using it straight crashes R.
  max_at_b1 <- function(beta1){
    LC <- function(beta2){-loglike_C(c(beta1, beta2), times, T_1, T_f, N)}
    mle <- optim(0.5, LC, method="L-BFGS-B", lower=0.01, upper= 0.99)
    return(-mle$value)}
  sampled_maxes = sapply(sampled_points, max_at_b1)
  return(approxfun(x=sampled_points, y=sampled_maxes, method="linear", rule=1))
}


get_loglike_prof_2 <- function(times, T_1 = 5, T_f = 10, N = 60, precision=0.05){
  sampled_points = seq(0+precision, 1, by=precision)
  # Returns a fitted approx curve of profile likelihood, on domain of 0 to 1. 
  # Doing nested optim is bad, so this "pre bakes" the function.
  # Yes, max_at_b1 is just the likelihood function, but using it straight crashes R.
  max_at_b2 <- function(beta2){
    LC <- function(beta1){-loglike_C(c(beta1, beta2), times, T_1, T_f, N)}
    mle <- optim(0.5, LC, method="L-BFGS-B", lower=0.01, upper= 0.99)
    return(-mle$value)}
  sampled_maxes = sapply(sampled_points, max_at_b2)
  return(approxfun(x=sampled_points, y=sampled_maxes, method="linear", rule=1))
}

# loglike_prof_1(beta1=0.01, times=c(0, 0.2, 0.8, 3,  9, 1.2, 3.4, 5.1, 6, 7), T_1 = 5, T_f = 10, N = 60)