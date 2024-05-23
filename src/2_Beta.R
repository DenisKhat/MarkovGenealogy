source('src/MarkhovChain.R')
source('src/Experiment.R')
# source('src/Trees_v2.R')
source('src/PiecewiseBeta.R')
library(gridExtra)
height = 1

N = 60
b1 = 0.6
b2 = 0.3
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

# first <- first[-1,]
last <- last[-1,] # to include/not include the 5
# print(last[, 1:6])

# PRINT GRAPH
full_table <- rbind(first, last)
p_0 <- full_table$I_list[1]
# phylog(full_table, p_0, end_time)

# print(full_table$time)
print(full_table[, 1:6])
infection_times <- full_table$time
likelihood_C <- function(beta) -loglike_C(beta, times = infection_times, T_1 = change_time, T_f=end_time, N = N)
likelihood_P1 <- get_loglike_prof_1(times=infection_times, T_1=change_time, T_f=end_time, N=N, precision = 0.001)
likelihood_P2 <- get_loglike_prof_2(times=infection_times, T_1=change_time, T_f=end_time, N=N, precision = 0.001)
# likelihood_C(c(0.1,0.4))
# likelihood_C(c(0.5, 0.5))
# print(likelihood_C(c(0.5, 0.5)))
# print(likelihood_C(c(1, 1)))

# two variable optimization
betas_hat = optim(c(0.5,0.5), likelihood_C, method = "L-BFGS-B", lower=c(0.01, 0.01), upper=c(0.99, 0.99))
betas_hat$par

#fix beta1, find max beta 2 (profile likelihood for beta 1)
beta1_hat = optim(0.5, function(x) -likelihood_P1(x), method = "L-BFGS-B", lower=0.005, upper=0.99)
beta1_hat$par
#fix beta2, find max beta 1 (profile likelihood for beta 2)
beta2_hat = optim(0.5, function(x) -likelihood_P2(x), method = "L-BFGS-B", lower=0.005, upper=0.99)
beta2_hat$par

chi <- qchisq(p = 0.95, df = 1)/2
shifted_L1 <- function(x) likelihood_P1(x) + beta1_hat$value + chi
shifted_L2 <- function(x) likelihood_P2(x) + beta2_hat$value + chi

# plot(likelihood_P1)
# plot(shifted_L1)

# confidence interval b1 (Wilks' estimate)
upper_b1 <- 1
lower_b1 <- 0
if (shifted_L1(0.01) < 0){
lower_b1 <- uniroot(shifted_L1, lower = 0.001, upper=beta1_hat$par)$root }
if (shifted_L1(1) <= 0){
  upper_b1 = uniroot(shifted_L1, lower = beta1_hat$par, upper = 0.99)$root }

#confidence interval b2 (Wilks' estimate)
upper_b2 <- 1
lower_b2 <- 0
if (shifted_L2(0.01) < 0){
  lower_b2 = uniroot(shifted_L2, lower = 0.001, upper=beta2_hat$par)$root}
if (shifted_L2(1) <= 0){
  upper_b2 = uniroot(shifted_L2, lower = beta2_hat$par, upper = 0.99)$root
}
xs <- seq(0.01, 1, by=0.01)
y1s <- sapply(xs, likelihood_P1)
y2s <- sapply(xs, likelihood_P2)
print(y2s)
dta <- data.frame(xs = xs, y1s = y1s, y2s = y2s)
b1_plot <- ggplot(dta, aes(x = xs, y = y1s)) +
    # geom_point() +
    geom_line() +
    annotate("segment", x = lower_b1, y = -Inf, xend = lower_b1, yend = -beta1_hat$value - chi, linetype = "dashed", color = "red") +
    annotate("segment", x = upper_b1, y = -Inf, xend = upper_b1, yend = -beta1_hat$value - chi, linetype = "dashed", color = "red") +
    annotate("segment", x = beta1_hat$par, y = -Inf, xend = beta1_hat$par, yend = -beta1_hat$value, color = "red") +
    annotate("segment", x = b1, y = -Inf, xend = b1, yend = likelihood_P1(b1), color = "blue") + 
    # coord_cartesian(xlim=c(0.01, 1), ylim=c(-100, y_max)) +
    # ggtitle("M") + 
    xlab('Beta1') +
    ylab('Profile Log-Likelihood')

#Plot b2
b2_plot <- ggplot(dta, aes(x = xs, y = y2s)) +
    # geom_point() +
    geom_line() +
    annotate("segment", x = lower_b2, y = -Inf, xend = lower_b2, yend = -beta2_hat$value - chi, linetype = "dashed", color = "red") +
    annotate("segment", x = upper_b2, y = -Inf, xend = upper_b2, yend = -beta2_hat$value - chi, linetype = "dashed", color = "red") +
    annotate("segment", x = beta2_hat$par, y = -Inf, xend = beta2_hat$par, yend = -beta2_hat$value, color = "red") +
    annotate("segment", x = b2, y = -Inf, xend = b2, yend = likelihood_P2(b2), color = "blue") + 
    # coord_cartesian(xlim=c(0.01, 1), ylim=c(-100, y_max)) +
    # ggtitle("M") + 
    xlab('Beta2') +
    ylab('Profile Log-Likelihood')


heatmap_C <- function(times, T1, TF, N){
  b1 <- seq(0.01, 1, by = 0.01)
  b2 <- seq(0.01, 1, by = 0.01)
  
  LF <- function(b1, b2){
    betas <- c(b1, b2)
    loglike_C(betas, times = times, T_1 = T1, T_f=TF, N = N)
  } 
  matrix_result <- outer(b1, b2, Vectorize(LF))
  print(matrix_result)
  
  # NORMALIZE
  # min_val <- min(matrix_result)
  # max_val <- max(matrix_result)
  # normalized_result <- (matrix_result - min_val) / (max_val - min_val)
  
  data <- expand.grid(b1 = b1, b2 = b2)
  # data$Z <- as.vector(normalized_result)
  data$Z <- as.vector(matrix_result)

  p <- ggplot(data, aes(b1, b2, fill= Z)) +
    geom_tile()+
    scale_fill_gradient(low="blue", high="red")
  print(p)
}

# b1_and_b2_plot <- heatmap(infection_times, change_time, end_time, N)

# METHOD C
# heatmap_C(infection_times, change_time, end_time, N)

grid.arrange( b1_plot, b2_plot, nrow = 2)
