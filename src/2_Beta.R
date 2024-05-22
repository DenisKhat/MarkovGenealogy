source('src/MarkhovChain.R')
source('src/Experiment.R')
source('src/Trees_v2.R')
source('src/PiecewiseBeta.R')
height = 1


N = 60
b1 = 0.3
b2 = 0.6
gamma = 0
change_time = 5
end_time = 10

first <- markhov_virus(change_time, b1, gamma, 59, 1)
# print(first[, 1:6])
first_R <- as.numeric(first$R[nrow(first)])
first_S <- as.numeric(first$S[nrow(first)])
first_I <- as.numeric(first$I[nrow(first)])
# print(first$S_list)


last <- markhov_virus(end_time, b2, gamma, first_S, first_I, first_R, curr_time=change_time, S_list=first$S_list[nrow(first)], I_list = first$I_list[nrow(first)], R_list=first$R_list[nrow(first)])

first <- first[-1,]
last <- last[-1,] # to include/not include the 5
# print(last[, 1:6])

full_table <- rbind(first, last)
p_0 <- full_table$infector[2]
# phylog(full_table, p_0, end_time)

# print(full_table$time)
# print(full_table[, 1:6])
infection_times <- full_table$time
likelihood_C <- function(beta) -loglike_C(beta, times = infection_times, T_1 = change_time, T_f=end_time, N = N)
likelihood_P1 <- get_loglike_prof_1(times=infection_times, T_1=change_time, T_f=end_time, N=N, precision = 0.01)
likelihood_P2 <- get_loglike_prof_2(times=infection_times, T_1=change_time, T_f=end_time, N=N, precision = 0.01)
# likelihood_C(c(0.1,0.4))
# likelihood_C(c(0.5, 0.5))
# print(likelihood_C(c(0.5, 0.5)))
# print(likelihood_C(c(1, 1)))

# two variable optimization
betas_hat = optim(c(0.5,0.5), likelihood_C, method = "L-BFGS-B", lower=c(0.01, 0.01), upper=c(0.99, 0.99))
betas_hat$par

#fix beta1, find max beta 2 (profile likelihood for beta 1)
beta1_hat = optim(0.5, function(x) -likelihood_P1(x), method = "L-BFGS-B", lower=0.01, upper=0.99)
beta1_hat$par
plot(likelihood_P1)
#fix beta2, find max beta 1 (profile likelihood for beta 2)
beta2_hat = optim(0.5, function(x) -likelihood_P2(x), method = "L-BFGS-B", lower=0.01, upper=0.99)
beta2_hat$par
plot(likelihood_P2)
