library(deSolve)
# source("src/MarkhovChain.R")

end_time = 11
N = 60

# parameters <- list(beta = 0.7, gamma = 0)

initial_state_sir <- c(S = 59/60, I = 1/60, R = 0)
initial_state_sis <- c(S = 59/60, I = 1/60)

# Time vector for the simulation
times <- seq(0, end_time, by = 0.01)


SIS_model <- function(time, state, parms) { # Using backwards in time A for system
  with(as.list(c(state, parms)), {
    dS <- (-beta * S * I + gamma * I)
    dI <- (beta * S * I - gamma * I)
    return(list(c(dS, dI)))
  })
}

SIR_model <- function(time, state, parms) {
  with(as.list(c(state, parms)), {
    # print(parms)
    dS <- (-beta * S * I)
    dI <- (beta * S * I - gamma * I)
    dR <- (gamma * I)
    
    return(list(c(dS, dI, dR)))
  })
}
# 
# SIRoutput <- ode(y = initial_state_sir, times = times, func = SIR_model, parms = parameters)
# SIR_df <- as.data.frame((SIRoutput))
# 
# SIR_df$A <- SIR_df$I
# 
# print(SIR_df)
# 
# ggplot(SIR_df, aes(x = times)) +
#     geom_line(aes(y = S, color = "Susceptible")) +
#     geom_line(aes(y = I, color = "Infected")) +
#     geom_line(aes(y = R, color = "Recovered")) +
#     labs(y = "Number of Individuals / Ancestor Function", color = "Compartment") +
#     ggtitle("SIR Model and Ancestor Function")

# table <- markhov_virus(end_time=10,beta=0.7, gamma=0, S=59, I=1)
# tableTimes <- as.vector(table$time)
# print(tableTimes)

Abar <- function(data, beta, t){
  f_si <- beta * (data$S[data$time == t]) * (data$I[data$time == t])
  return(f_si)
}

Volz_SI <- function(beta, sample_times){
  parameters <- c(beta = beta, gamma = 0)
  SIRoutput <- ode(y = initial_state_sir, times = times, func = SIR_model, parms = parameters)
  SIR_df <- as.data.frame((SIRoutput))
  SIR_df$A <- SIR_df$I 
  
  n = length(SIR_df)
  AT = SIR_df$A[SIR_df$time == 10]
  At = SIR_df$A[SIR_df$time == 0]
  
  res <- -(n - 1)*log(AT - At)
  for (t in sample_times){
    t <- round(t, 0.01)
    res <- res + log(Abar(SIR_df, beta, t))
  }
  return(res)
}

Volz_SIS <- function(beta, gamma, sample_times){
  parameters <- c(beta = beta, gamma = gamma)
  SISoutput <- as.data.frame(ode(y = initial_state_sis, times = times, func = SIS_model, parms = parameters))
  # print(SISoutput$S)
  f_SI_points <- sapply(seq(1, end_time/0.01 + 1, by=1), function(i) beta * SISoutput$S[i] * SISoutput$I[i] )
  I_func <- approxfun(x=times, y=as.vector(SISoutput$I))
  f_SI <- approxfun(x=times, y=f_SI_points)
  # print("passed f_SI")
  A_initial <- c(A=I_func(10))
  A_model <- function(time, state, parms){ # parms should be nothing...
    with(as.list(c(state, parms)), {
      dA <- - I_func(10 - time) * ((A/I(10-time))^2)
      return(list(c(dA)))
    })
  }
  A_ode <- as.data.frame(ode(y=A_initial,times=seq(0, 10, by=0.01), func=A_model, parms = c(), method="ode45"))
  A_ode$A <- rev(A_ode$A)
  A_func <- approxfun(x = seq(0,10,by=0.01),y=as.vector(A_ode$A))
  # print("hi?")
  dA <- function(x) log(f_SI(x)*(A_func(x)/I_func(x))^2)
  # print(f_SI(0.5))
  # print(A_func(0.5))
  # print(I_func(0.5))
  # print(dA(0.5))
  # n <- length(sample_times)
  likelihood_times <- sapply(sample_times, dA)
  # print(sum(likelihood_times) - (n-1)*log(A_func(10)))
  # print(sum(likelihood_times) - (n-1)*log(A_func(10)-A_func(0.01)))
  return(sum(likelihood_times) + (n-1)*log(A_func(10)-A_func(0.01)))
} 

### TESTING CODE FOR VOLZ SIS
# 
# sampled_times <- markhov_virus(10, 0.7, 0.3, S=59, I=1, SIS=T)
# n <- nrow(sampled_times)
# coalescent_times <- sampled_times[2:n,]
# coalescent_times <- subset(coalescent_times, affected > 0)
# coalescent_times <- coalescent_times$time
# coalescent_times
# Volz_SIS(beta = 0.7, gamma=0.3, coalescent_times)

# SIS_like <- function(theta) -Volz_SIS(theta[1],theta[2],sample_times = coalescent_times)
# SIS_like <- function(theta) -Volz_SIS(theta,0.3,sample_times = coalescent_times)
# MLE <- optim(par = 0.6, fn= SIS_like, method="L-BFGS-B", lower = 0.01, upper = 2)
# MLE <- optim(par = c(0.5, 0.5), fn= SIS_like, method="L-BFGS-B", lower = c(0.01,0.01), upper = c(2,2))
# MLE$par
###

Volz_SIR <- function(beta, gamma, sample_times){
  parameters <- c(beta = beta, gamma = gamma)
  SIRoutput <- as.data.frame(ode(y = initial_state_sir, times = times, func = SIR_model, parms = parameters))
  # print(SISoutput$S)
  f_SI_points <- sapply(seq(1, end_time/0.01 + 1, by=1), function(i) beta * SIRoutput$S[i] * SIRoutput$I[i] )
  I_func <- approxfun(x=times, y=as.vector(SIRoutput$I))
  f_SI <- approxfun(x=times, y=f_SI_points)
  # print("passed f_SI")
  A_initial <- c(A=I_func(10))
  A_model <- function(time, state, parms){ # parms should be nothing...
    with(as.list(c(state, parms)), {
      dA <- - I_func(10 - time) * ((A/I(10-time))^2)
      return(list(c(dA)))
    })
  }
  A_ode <- as.data.frame(ode(y=A_initial,times=seq(0, 10, by=0.01), func=A_model, parms = c(), method="ode45"))
  # A_ode$A <- rev(A_ode$A)
  # I'm lazy so this is A backwards in time.
  A_func <- approxfun(x = seq(0,10,by=0.01),y=as.vector(A_ode$A))
  # print("hi?")
  dA <- function(x) log(f_SI(x)*(A_func(10-x)/I_func(x))^2)
  # print(f_SI(0.5))
  # print(A_func(0.5))
  # print(I_func(0.5))
  # print(dA(0.5))
  # n <- length(sample_times)
  likelihood_times <- sapply(sample_times, dA)
  # print(sum(likelihood_times) - (n-1)*log(A_func(10)))
  # print(sum(likelihood_times) - (n-1)*log(A_func(10)-A_func(0.01)))
  return(sum(likelihood_times) + (n-1)*log(A_func(0.01)-A_func(9.99)))
} 



SIR_like <- function(theta) -Volz_SIR(beta=theta[1],gamma=theta[2],sample_times = coalescent_times)
# SIR_like(c(0.7, 0.3))
# SIR_like(c(0.7, 0.1))
# SIR_like(c(1, 0.3))
# SIR_like(c(0.1,2))
# MLE <- optim(par = c(0.5, 0.5), fn= SIR_like, method="L-BFGS-B", lower = c(0.01,0.01), upper = c(1,1))
# MLE$par