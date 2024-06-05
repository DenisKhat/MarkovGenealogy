source("slurm_simulation/virusSimulation.R")  # Uncomment locally
# source("virusSimulation") # Uncomment on slurm
dir <- getwd()
job_id = 1  # Uncomment locally
# args = commandArgs(TRUE)  # Uncomment on slurm
# job_id = as.numeric(args[1])  # Uncomment on slurm
m <- ceiling(params$nsim/params$ncores)
sim_numbers <- ((job_id-1)*m + 1):(min(job_id*m,params$nsim)) 
  
get_rate <- function(beta) return(function(i) beta * (N-i) * i / N)


log_likelihood <- function(beta, data){
  times <- data$time
  M <- length(times)
  if (times[1] == 0.0){
    times <- times[2:M]
    M <- M - 1
  }
  
  rate <- get_rate(beta)
  out <- 0
  
  # Nothing happens before T_f
  if (M==0){
    out <- out - rate(1)*(params$T_f)
  }
  # SOMETHING happens before T_f
  else {
    out <- out + log(rate(1)) - rate(1)*(times[1])
    if(M > 1){
      for (i in seq(2,M)){
        out <- out + log(rate(i)) - rate(i)*(times[i]-times[i-1])
      }
    }
    # stopped sampling before we infected everyone
    if (M < N){
      out <- out - rate(M+1)*(params$T_f - times[M])
    }
    return(out)
  }
}


full_knowledge_estimation <- function(data, alpha=0.05){
  LC <- function(beta) -log_likelihood(beta, data)
  MLE <- optim(0.5, LC, method = "L-BFGS-B", lower = 0.01, upper= 0.99, hessian = FALSE)
  in_interval <- 1
  chi <- qchisq(p = 1-alpha, df = 1)
  wilks_cutoff <-  -MLE$val - chi/2
  x_points <- seq(0.01, 1.01, by=0.01)
  y_points <- sapply(x_points, function(x) -LC(x) - wilks_cutoff)
  curve <- approxfun(x_points, y_points)
  # the below if is when simulation has no infections
  if (MLE$par == 0.01){
    interval_width <- NA
    in_interval <- NA
    MLE$par <- 0
  }
  else{
    lower_wilks = uniroot(curve, lower = 0.01, upper=MLE$par)$root
    upper_wilks = uniroot(curve, lower = MLE$par, upper = 1)$root
    interval_width <- upper_wilks - lower_wilks
    if (params$beta < lower_wilks || params$beta > upper_wilks){
      in_interval <- 0
    }
  }
  
  return(list(mle=MLE$par, in_interval=in_interval, interval_width=interval_width))
}


estimate <- c()
in_interval <- c()
ids <- c()
int_widths <- c()

whole_table <- read.csv(file = "simulation1beta.csv")

# for (i in 1:5000){
#   if (length(which(whole_table$id==i)) == 0){
#     print(i)
#   }
# }
for (i in sim_numbers){
  # dta <- markov_virus(params$beta)
  dta <- whole_table[which(whole_table$id == i),]
  row_data <- full_knowledge_estimation(dta)
  ids <- c(ids, i)
  estimate <- c(estimate, row_data$mle)
  in_interval <- c(in_interval, row_data$in_interval)
  int_widths <- c(int_widths, row_data$interval_width)
}


path = file.path(dir, "Results")
if (!file.exists(path)){
  dir.create(path)
}


output_table <- data.frame(id=ids, estimate=estimate, in_interval=in_interval, interval_width=int_widths)
write.csv(output_table,file=paste0("Results/full_knowledge_estimate_",job_id,".csv"),row.names=F)
