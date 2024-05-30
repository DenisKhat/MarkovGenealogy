library(deSolve)
source("src/MarkhovChain.R")

end_time = 10
N = 60

# parameters <- list(beta = 0.7, gamma = 0)

initial_state_sir <- c(S = 59, I = 1, R = 0)

# Time vector for the simulation
times <- seq(0, end_time, by = 0.01)


SIR_model <- function(time, state, parms) {
  with(as.list(c(state, parms)), {
    # print(parms)
    dS <- (-beta * S * I)/N
    dI <- (beta * S * I - gamma * I)/N
    dR <- (gamma * I)/N
    
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

table <- markhov_virus(end_time=10,beta=0.7, gamma=0, S=59, I=1)
tableTimes <- as.vector(table$time)
print(tableTimes)

Abar <- function(data, beta, t){
  f_si <- beta * (data$S[data$time == t]) * (data$I[data$time == t])
  return(f_si)
}

Volz_likelihood <- function(beta){
  parameters <- c(beta = beta, gamma = 0)
  SIRoutput <- ode(y = initial_state_sir, times = times, func = SIR_model, parms = parameters)
  SIR_df <- as.data.frame((SIRoutput))
  
  SIR_df$A <- SIR_df$I
  
  n = length(SIR_df)
  AT = SIR_df$A[SIR_df$time == 10]
  At = SIR_df$A[SIR_df$time == 0]
  
  res <- -(n - 1)*log(AT - At)
  for (t in tableTimes){
    t <- round(t, 0.01)
    res <- res + log(Abar(SIR_df, beta, t))
  }
  return(res)
}

# plot(Volz_likelihood)

volz <- function(beta) -Volz_likelihood(beta)


mle <- optim(0.5, volz, method="L-BFGS-B", lower=0.01, upper= 1)

print(mle)


