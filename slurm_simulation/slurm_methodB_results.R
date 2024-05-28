methodB_df <- read.csv("slurm_simulation/combined.csv")

estimate <- mean(methodB_df$estimate)
CI <- mean(methodB_df$inCI, na.rm = TRUE)
width <- mean(methodB_df$width, na.rm = TRUE)
length <- length(methodB_df[,1])

paste("avg estimate:", estimate)
paste("CI:", CI)
paste("avg width:", width)
paste("num of sims:", length)
