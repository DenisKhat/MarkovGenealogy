methodB_df <- read.csv("slurm_simulation/combined.csv")

twoBetaMethodB_df <- read.csv("slurm_simulation/2betaMethodBcombined.csv")


# 1 beta 
estimate <- mean(methodB_df$estimate)
CI <- mean(methodB_df$inCI, na.rm = TRUE)
width <- mean(methodB_df$width, na.rm = TRUE)
tot_length <- length(methodB_df[,1])
sims_wout_outliers <- methodB_df[!is.na(methodB_df$inCI), ]
length_sims_wout_outliers <- length(sims_wout_outliers[,1])


paste("avg estimate:", estimate)
paste("CI:", CI)
paste("avg width:", width)
paste("num of total sims:", tot_length)
paste("num of sims wout outliers:", length_sims_wout_outliers)


# 2beta

b1_estimate <- mean(twoBetaMethodB_df$b1_estimate)
b2_estimate <- mean(twoBetaMethodB_df$b2_estimate)
b1_CI <- mean(twoBetaMethodB_df$b1_inCI_col, na.rm = TRUE)
b2_CI <- mean(twoBetaMethodB_df$b2_inCI_col, na.rm = TRUE)
b1_width <- mean(twoBetaMethodB_df$b1_CI_width_col, na.rm = TRUE)
b2_width <- mean(twoBetaMethodB_df$b2_CI_width_col, na.rm = TRUE)
tot_length2 <- length(twoBetaMethodB_df[,1])
sims_wout_outliers2 <- twoBetaMethodB_df[!is.na(twoBetaMethodB_df$b1_CI) & !is.na(twoBetaMethodB_df$b2_CI), ]

length_sims_wout_outliers2 <- length(sims_wout_outliers2[,1])

paste("b1 estimate:", b1_estimate)
paste("b2 estimate:", b2_estimate)
paste("b1_in_CI:", b1_CI)
paste("b2_in_CI:", b2_CI)
paste("b1 width:", b1_width)
paste("b2 width:", b2_width)
paste("num of total sims:", tot_length2)
paste("num of sims wout outliers:", length_sims_wout_outliers2)
