# This is the function that compares accuracy of rank estimation.
# The methods we compare are profile likelihood (PL), edge distribution (ED) and permutation method used in JIVE
source("../../../MyFunction/Angle_Calculation.R")
source("../../../MyFunction/Profile_Likelihood_Rank_Selection.R")
source("../../../MyFunction/DoubleMatchedMatrixDecomposition.R")
source("../../../MyFunction/FindOptMatrix.R")
source("../../../MyFunction/Preliminary_Functions.R")
source("../../../MyFunction/Select_ED_Rank.R")

library(foreach)
library(doParallel)

# Get the generated data
load("../Data/Data.RData")

set.seed(37)

n = 88
p = 736
nrep = 100
total_rank1 = rep(8,nrep)
total_rank2 = rep(6,nrep)

joint_rank_col = rep(0,nrep)
joint_rank_row = rep(2,nrep)

cl = makeCluster(28)
registerDoParallel(cl)

# Main for loop

output <- foreach (i = 1:nrep) %dopar% {
  library(ajive)
  X1 = X1_list[[i]]
  X2 = X2_list[[i]]
  signal1 = signal1_list[[i]]
  signal2 = signal2_list[[i]]
  
  # PL method
  my_result_PL = DMMD_v2(X1,X2,r1 = total_rank1[i], total_rank1[i])
  # Estimation error of joint rank
  my_error_joint_rank_col_PL = my_result_PL$Rank$"Joint Column Rank" - joint_rank_col[i]
  my_error_joint_rank_row_PL = my_result_PL$Rank$"Joint Row Rank" - joint_rank_row[i]

  data_col = list(X1,X2)
  data_row = list(t(X1),t(X2))
  initial_signal_ranks <- c(total_rank1[i], total_rank2[i])
  # Run AJIVE on the row direction
  ajive_result_row <- ajive(data_row, initial_signal_ranks)
  # Estimation error of joint rank
  ajive_error_joint_rank_r = ajive_result_row$joint_rank - joint_rank_row[i]
  
  # Run AJIVE on the column direction
  ajive_result_col <- ajive(data_col, initial_signal_ranks)
  # Estimation error of joint rank
  ajive_error_joint_rank_c = ajive_result_col$joint_rank - joint_rank_col[i]

  list(my_error_joint_rank_col_PL = my_error_joint_rank_col_PL,
       my_error_joint_rank_row_PL = my_error_joint_rank_row_PL,
       ajive_error_joint_rank_c = ajive_error_joint_rank_c,
       ajive_error_joint_rank_r = ajive_error_joint_rank_r)

}

stopCluster(cl)

save(output,total_rank1,total_rank2,joint_rank_col,joint_rank_row, file = "AJIVEoutput.RData")




