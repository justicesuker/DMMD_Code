# This is the function that compares accuracy of signal estimation given true ranks between DMMD and JIVE. 
rm(list=ls())
source("../MyFunction/Preliminary_Functions.R")
source("../MyFunction/Angle_Calculation.R")
source("../MyFunction/Profile_Likelihood_Rank_Selection.R")
source("../MyFunction/DoubleMatchedMatrixDecomposition.R")
source("../MyFunction/FindOptMatrix.R")

library(foreach)
library(doParallel)

# Get the generated data
load("../Data_Setting4_FixRank/Data.RData")

set.seed(37)
n = 240
p = 200
nrep = 140
total_rank1 = rep(20,nrep)
total_rank2 = rep(18,nrep)

joint_rank_col = rep(4,nrep)
joint_rank_row = rep(3,nrep)

angles <- function(Z, Zhat){
  svdZ = svd(Z)
  svdZhat = svd(Zhat)
  rZ = sum(svdZ$d > 1e-6)
  rZhat = sum(svdZhat$d > 1e-6)
  k = min(rZ, rZhat)
  cosines = svd(crossprod(svdZ$u[, 1:rZ], svdZhat$u[, 1:rZhat]))$d[1:k]
  return(list(cosines = cosines, dist = sqrt(sum(1-cosines^2))/sqrt(k)))
}

cl = makeCluster(28)
registerDoParallel(cl)

# Main for loop
output_small <- foreach (i = 1:nrep, .errorhandling = 'pass') %dopar% {
  print(i)
  library(ajive)
  X1 = X1_list[[i]]
  X2 = X2_list[[i]]
  signal1 = signal1_list[[i]]
  signal2 = signal2_list[[i]]
  J1_r = J1r_list[[i]]
  J2_r = J2r_list[[i]]
  J1_c = J1c_list[[i]]
  J2_c = J2c_list[[i]]
  
  I1_r = signal1 - J1_r
  I2_r = signal2 - J2_r
  I1_c = signal1 - J1_c
  I2_c = signal2 - J2_c
  
  # AJIVE
  data_col = list(X1,X2)
  data_row = list(t(X1),t(X2))
  
  # Run AJIVE on the row direction - row-wise decomposition
  initial_signal_ranks <- c(total_rank1[i] - 1, total_rank2[i] - 1)
  
  ajive_result_row <- ajive(data_row, initial_signal_ranks, joint_rank = 2)
  
  # Run AJIVE on the column direction - column-wise decomposition
  ajive_result_col <- ajive(data_col, initial_signal_ranks, joint_rank = 3)
  
  ajive_J1r <- ajive_result_row$block_decomps[[1]]$joint$full
  ajive_J2r <- ajive_result_row$block_decomps[[2]]$joint$full
  ajive_I1r <- ajive_result_row$block_decomps[[1]]$individual$full
  ajive_I2r <- ajive_result_row$block_decomps[[2]]$individual$full
  
  ajive_J1c <- ajive_result_col$block_decomps[[1]]$joint$full
  ajive_J2c <- ajive_result_col$block_decomps[[2]]$joint$full
  ajive_I1c <- ajive_result_col$block_decomps[[1]]$individual$full
  ajive_I2c <- ajive_result_col$block_decomps[[2]]$individual$full
  
  # The estimation error of signals in row-wise AJIVE 
  ajive_row_error1 = t(ajive_J1r) + t(ajive_I1r) - signal1
  ajive_row_error2 = t(ajive_J2r) + t(ajive_I2r) - signal2
  
  # Relative error
  ajive_row_error1 = Fnorm(ajive_row_error1)^2/Fnorm(signal1)^2
  ajive_row_error2 = Fnorm(ajive_row_error2)^2/Fnorm(signal2)^2
  
  # The estimation error of signals in column-wise AJIVE 
  ajive_col_error1 = ajive_J1c + ajive_I1c - signal1
  ajive_col_error2 = ajive_J2c + ajive_I2c - signal2
  
  # Relative error
  ajive_col_error1 = Fnorm(ajive_col_error1)^2/Fnorm(signal1)^2
  ajive_col_error2 = Fnorm(ajive_col_error2)^2/Fnorm(signal2)^2
  
  # Record the errors of row decomposition in AJIVE
  ajive_row_joint_error1 = angles(ajive_J1r + ajive_I1r, t(J1_r))$dist
  ajive_row_joint_error2 = angles(ajive_J2r + ajive_I2r, t(J2_r))$dist
  
  # Record the errors of column decomposition in JIVE
  ajive_col_joint_error1 = angles(ajive_J1c + ajive_I1c, J1_c)$dist
  ajive_col_joint_error2 = angles(ajive_J2c + ajive_I2c, J2_c)$dist

  list(ajive_row_error1 = ajive_row_error1,
       ajive_row_error2 = ajive_row_error2,
       ajive_row_joint_error1 = ajive_row_joint_error1,
       ajive_row_joint_error2 = ajive_row_joint_error2,
       ajive_col_error1 = ajive_col_error1,
       ajive_col_error2 = ajive_col_error2,
       ajive_col_joint_error1 = ajive_col_joint_error1,
       ajive_col_joint_error2 = ajive_col_joint_error2)
}

output_large <- foreach (i = 1:nrep, .errorhandling = 'pass') %dopar% {
  print(i)
  library(ajive)
  X1 = X1_list[[i]]
  X2 = X2_list[[i]]
  signal1 = signal1_list[[i]]
  signal2 = signal2_list[[i]]
  J1_r = J1r_list[[i]]
  J2_r = J2r_list[[i]]
  J1_c = J1c_list[[i]]
  J2_c = J2c_list[[i]]
  
  I1_r = signal1 - J1_r
  I2_r = signal2 - J2_r
  I1_c = signal1 - J1_c
  I2_c = signal2 - J2_c
  
  # AJIVE
  data_col = list(X1,X2)
  data_row = list(t(X1),t(X2))
  
  # Run AJIVE on the row direction - row-wise decomposition
  initial_signal_ranks <- c(total_rank1[i] + 1, total_rank2[i] + 1)
  
  ajive_result_row <- ajive(data_row, initial_signal_ranks, joint_rank = 3 + 1)
  
  # Run AJIVE on the column direction - column-wise decomposition
  ajive_result_col <- ajive(data_col, initial_signal_ranks, joint_rank = 4 + 1)
  
  ajive_J1r <- ajive_result_row$block_decomps[[1]]$joint$full
  ajive_J2r <- ajive_result_row$block_decomps[[2]]$joint$full
  ajive_I1r <- ajive_result_row$block_decomps[[1]]$individual$full
  ajive_I2r <- ajive_result_row$block_decomps[[2]]$individual$full
  
  ajive_J1c <- ajive_result_col$block_decomps[[1]]$joint$full
  ajive_J2c <- ajive_result_col$block_decomps[[2]]$joint$full
  ajive_I1c <- ajive_result_col$block_decomps[[1]]$individual$full
  ajive_I2c <- ajive_result_col$block_decomps[[2]]$individual$full
  
  # The estimation error of signals in row-wise AJIVE 
  ajive_row_error1 = t(ajive_J1r) + t(ajive_I1r) - signal1
  ajive_row_error2 = t(ajive_J2r) + t(ajive_I2r) - signal2
  
  # Relative error
  ajive_row_error1 = Fnorm(ajive_row_error1)^2/Fnorm(signal1)^2
  ajive_row_error2 = Fnorm(ajive_row_error2)^2/Fnorm(signal2)^2
  
  # The estimation error of signals in column-wise AJIVE 
  ajive_col_error1 = ajive_J1c + ajive_I1c - signal1
  ajive_col_error2 = ajive_J2c + ajive_I2c - signal2
  
  # Relative error
  ajive_col_error1 = Fnorm(ajive_col_error1)^2/Fnorm(signal1)^2
  ajive_col_error2 = Fnorm(ajive_col_error2)^2/Fnorm(signal2)^2
  
  # Record the errors of row decomposition in AJIVE
  ajive_row_joint_error1 = angles(ajive_J1r + ajive_I1r, t(J1_r))$dist
  ajive_row_joint_error2 = angles(ajive_J2r + ajive_I2r, t(J2_r))$dist
  
  # Record the errors of column decomposition in JIVE
  ajive_col_joint_error1 = angles(ajive_J1c + ajive_I1c, J1_c)$dist
  ajive_col_joint_error2 = angles(ajive_J2c + ajive_I2c, J2_c)$dist
  
  list(ajive_row_error1 = ajive_row_error1,
       ajive_row_error2 = ajive_row_error2,
       ajive_row_joint_error1 = ajive_row_joint_error1,
       ajive_row_joint_error2 = ajive_row_joint_error2,
       ajive_col_error1 = ajive_col_error1,
       ajive_col_error2 = ajive_col_error2,
       ajive_col_joint_error1 = ajive_col_joint_error1,
       ajive_col_joint_error2 = ajive_col_joint_error2)
}

output_jointsmall <- foreach (i = 1:nrep, .errorhandling = 'pass') %dopar% {
  print(i)
  library(ajive)
  X1 = X1_list[[i]]
  X2 = X2_list[[i]]
  signal1 = signal1_list[[i]]
  signal2 = signal2_list[[i]]
  J1_r = J1r_list[[i]]
  J2_r = J2r_list[[i]]
  J1_c = J1c_list[[i]]
  J2_c = J2c_list[[i]]
  
  I1_r = signal1 - J1_r
  I2_r = signal2 - J2_r
  I1_c = signal1 - J1_c
  I2_c = signal2 - J2_c
  
  # AJIVE
  data_col = list(X1,X2)
  data_row = list(t(X1),t(X2))
  
  # Run AJIVE on the row direction - row-wise decomposition
  initial_signal_ranks <- c(total_rank1[i], total_rank2[i])
  
  ajive_result_row <- ajive(data_row, initial_signal_ranks, joint_rank = 3 - 1)
  
  # Run AJIVE on the column direction - column-wise decomposition
  ajive_result_col <- ajive(data_col, initial_signal_ranks, joint_rank = 4 - 1)
  
  ajive_J1r <- ajive_result_row$block_decomps[[1]]$joint$full
  ajive_J2r <- ajive_result_row$block_decomps[[2]]$joint$full
  ajive_I1r <- ajive_result_row$block_decomps[[1]]$individual$full
  ajive_I2r <- ajive_result_row$block_decomps[[2]]$individual$full
  
  ajive_J1c <- ajive_result_col$block_decomps[[1]]$joint$full
  ajive_J2c <- ajive_result_col$block_decomps[[2]]$joint$full
  ajive_I1c <- ajive_result_col$block_decomps[[1]]$individual$full
  ajive_I2c <- ajive_result_col$block_decomps[[2]]$individual$full
  
  # The estimation error of signals in row-wise AJIVE 
  ajive_row_error1 = t(ajive_J1r) + t(ajive_I1r) - signal1
  ajive_row_error2 = t(ajive_J2r) + t(ajive_I2r) - signal2
  
  # Relative error
  ajive_row_error1 = Fnorm(ajive_row_error1)^2/Fnorm(signal1)^2
  ajive_row_error2 = Fnorm(ajive_row_error2)^2/Fnorm(signal2)^2
  
  # The estimation error of signals in column-wise AJIVE 
  ajive_col_error1 = ajive_J1c + ajive_I1c - signal1
  ajive_col_error2 = ajive_J2c + ajive_I2c - signal2
  
  # Relative error
  ajive_col_error1 = Fnorm(ajive_col_error1)^2/Fnorm(signal1)^2
  ajive_col_error2 = Fnorm(ajive_col_error2)^2/Fnorm(signal2)^2
  
  # Record the errors of row decomposition in AJIVE
  ajive_row_joint_error1 = angles(ajive_J1r + ajive_I1r, t(J1_r))$dist
  ajive_row_joint_error2 = angles(ajive_J2r + ajive_I2r, t(J2_r))$dist
  
  # Record the errors of column decomposition in JIVE
  ajive_col_joint_error1 = angles(ajive_J1c + ajive_I1c, J1_c)$dist
  ajive_col_joint_error2 = angles(ajive_J2c + ajive_I2c, J2_c)$dist
  
  list(ajive_row_error1 = ajive_row_error1,
       ajive_row_error2 = ajive_row_error2,
       ajive_row_joint_error1 = ajive_row_joint_error1,
       ajive_row_joint_error2 = ajive_row_joint_error2,
       ajive_col_error1 = ajive_col_error1,
       ajive_col_error2 = ajive_col_error2,
       ajive_col_joint_error1 = ajive_col_joint_error1,
       ajive_col_joint_error2 = ajive_col_joint_error2)
}

stopCluster(cl)
save(output_small, output_large, output_jointsmall, file = "AJIVE_output.RData")
