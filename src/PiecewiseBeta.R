
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
  else{
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

loglike_C(betas=c(0.5, 0.5), times=c(0, 0.2, 0.8, 3,  9, 1.2, 3.4, 5.1, 6, 7), T_1 = 5, T_f = 10, N = 60)
