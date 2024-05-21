source('src/MarkhovChain.R')
source('src/Experiment.R')
source('src/Trees_v2.R')
height = 1


N = 60
b1 = .99
b2 = 0.1
gamma = 0.3
change_time = 5
end_time = 10

first <- markhov_virus(change_time, b1, gamma, 59, 1)
# print(first[, 1:6])
first_R <- as.numeric(first$R[nrow(first)])
first_S <- as.numeric(first$S[nrow(first)])
first_I <- as.numeric(first$I[nrow(first)])
# print(first$S_list)


last <- markhov_virus(end_time, b2, gamma, first_S, first_I, first_R, curr_time=change_time, S_list=first$S_list[nrow(first)], I_list = first$I_list[nrow(first)], R_list=first$R_list[nrow(first)])


# last <- last[-1, ]
print(last[, 1:6])

full_table <- rbind(first, last)
p_0 <- full_table$infector[2]
phylog(full_table, p_0, end_time)


print(full_table[, 1:6])

