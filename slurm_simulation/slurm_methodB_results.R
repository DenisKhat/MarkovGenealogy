methodB_df <- read.csv("slurm_simulation/combined.csv")

estimate <- mean(methodB_df$estimate)
CI <- mean(methodB_df$inCI, na.rm = TRUE)
width <- mean(methodB_df$width, na.rm = TRUE)
tot_length <- length(methodB_df[,1])
sims_w_events <- methodB_df[!is.na(methodB_df$inCI), ]
length_sims_w_events <- length(sims_w_events[,1])


paste("avg estimate:", estimate)
paste("CI:", CI)
paste("avg width:", width)
paste("num of total sims:", tot_length)
paste("num of sims w events:", length_sims_w_events)
