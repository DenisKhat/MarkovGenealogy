library(yaml)

params = read_yaml("slurm_simulation/params.yaml")
N <- params$num

markov_virus <- function(beta, I=params$I, t=0, end_t=params$T_f){ #SI model
  I_list <- c(I)
  t_list <- c(t)
  while (t < end_t & I < N){
    S = N - I
    next_time <- rexp(1, beta * S *I / N)
    t <- t + next_time
    I <- I + 1
    I_list <- c(I_list, I)
    t_list <- c(t_list, t)
  }
  table <- data.frame(time=t_list, I=I_list)
  if (t > end_t){
    table <- head(table, -1)
  }
  return(table)
}