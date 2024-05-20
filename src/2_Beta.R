source('src/MarkhovChain.R')
source('src/Experiment.R')

N = 60
b1 = 0.7
b2 = 0.4
change_time = 5
end_time = 10

first <- markhov_virus(change_time, b1, 0.2, 59, 1)
print(first[, 1:6])
first_R <- as.numeric(first$R[nrow(first)])
first_S <- as.numeric(first$S[nrow(first)])
first_I <- as.numeric(first$I[nrow(first)])
# 
# print(first_R)
# print(first_S)
# print(first_I)
# 
# print(first$R_list)
# print(first$S_list)
# print(first$I_list)



last <- markhov_virus(end_time, b2, first_R, first_S, first_I, curr_time=change_time, S_list=first$S_list[nrow(first)], I_list = first$I_list[nrow(first)], R_list=first$R_list[nrow(first)])
print(last[, 1:6])
# print(last$R_list)
# print(last$S_list)
# print(last$I_list)
