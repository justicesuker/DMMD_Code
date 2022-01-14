# This is the function that compares accuracy of signal estimation given true ranks between DMMD and JIVE. 
rm(list=ls())

source("../../../MyFunction/Angle_Calculation.R")
source("../../../MyFunction/Profile_Likelihood_Rank_Selection.R")
source("../../../MyFunction/DoubleMatchedMatrixDecomposition.R")
source("../../../MyFunction/FindOptMatrix.R")
source("../../../MyFunction/Preliminary_Functions.R")


library(foreach)
library(doParallel)

# Get the generated data
load("../Data/Data1.RData")
load("../Data/Data2.RData")

set.seed(37)
n = 88
p = 736

nrep = 100
total_rank1 = rep(8,nrep)
total_rank2 = rep(6,nrep)

joint_rank_col = rep(0,nrep)
joint_rank_row = rep(2,nrep)

cl = makeCluster(24)
registerDoParallel(cl)

# Main for loop
output <- foreach (i = 1:nrep, .errorhandling = 'pass') %dopar% {
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
  
  # Run JIVE on the row direction - row-wise decomposition
  initial_signal_ranks <- c(total_rank1[i], total_rank2[i])
  
  ajive_result_row <- ajive(data_row, initial_signal_ranks, joint_rank = joint_rank_row[i])
  
  # Run JIVE on the column direction - column-wise decomposition
  ajive_result_col <- ajive(data_col, initial_signal_ranks, joint_rank = joint_rank_col[i])
  
  ajive_J1r <- ajive_result_row$block_decomps[[1]]$joint$full
  ajive_J2r <- ajive_result_row$block_decomps[[2]]$joint$full
  ajive_I1r <- ajive_result_row$block_decomps[[1]]$individual$full
  ajive_I2r <- ajive_result_row$block_decomps[[2]]$individual$full
  
  ajive_J1c <- ajive_result_col$block_decomps[[1]]$joint$full
  ajive_J2c <- ajive_result_col$block_decomps[[2]]$joint$full
  ajive_I1c <- ajive_result_col$block_decomps[[1]]$individual$full
  ajive_I2c <- ajive_result_col$block_decomps[[2]]$individual$full
  
  # The estimation error of signals in row-wise AJIVE 
  if (is.na(ajive_J1r)){
    ajive_row_joint_error1 =  Fnorm(J1_r)^2
    ajive_row_ind_error1 = Fnorm(t(ajive_I1r) - I1_r)^2
    ajive_row_error1 = t(ajive_I1r) - signal1
    ajive_row_error1 = Fnorm(ajive_row_error1)^2/Fnorm(signal1)^2
  }
  if (!is.na(ajive_J1r)){
    ajive_row_joint_error1 = Fnorm(t(ajive_J1r) - J1_r)^2
    if (is.na(ajive_I1r)){
      ajive_row_ind_error1 = Fnorm(I1_r)^2
      ajive_row_error1 = t(ajive_J1r) - signal1
      ajive_row_error1 = Fnorm(ajive_row_error1)^2/Fnorm(signal1)^2
    }
    if (!is.na(ajive_I1r)){
      ajive_row_ind_error1 = Fnorm(t(ajive_I1r) - I1_r)^2
      ajive_row_error1 = t(ajive_J1r) + t(ajive_I1r) - signal1
      ajive_row_error1 = Fnorm(ajive_row_error1)^2/Fnorm(signal1)^2
    }
  }
  
  if (is.na(ajive_J2r)){
    ajive_row_joint_error2 =  Fnorm(J2_r)^2
    ajive_row_ind_error2 = Fnorm(t(ajive_I2r) - I2_r)^2
    ajive_row_error2 = t(ajive_I2r) - signal2
    ajive_row_error2 = Fnorm(ajive_row_error2)^2/Fnorm(signal2)^2
  }
  if (!is.na(ajive_J2r)){
    ajive_row_joint_error2 = Fnorm(t(ajive_J2r) - J2_r)^2
    # In this case, ajive_I2r should be 0 because of the design.
    if (is.na(ajive_I2r)){
      ajive_row_ind_error2 = Fnorm(I2_r)^2
      ajive_row_error2 = t(ajive_J2r) - signal2
      ajive_row_error2 = Fnorm(ajive_row_error2)^2/Fnorm(signal2)^2
    }
    if (!is.na(ajive_I2r)){
      ajive_row_ind_error2 = Fnorm(t(ajive_I2r) - I2_r)^2
      ajive_row_error2 = t(ajive_J2r) + t(ajive_I2r) - signal2
      ajive_row_error2 = Fnorm(ajive_row_error2)^2/Fnorm(signal2)^2
    }
  }
  
  # The estimation error of signals in column-wise AJIVE 
  if (is.na(ajive_J1c)){
    ajive_col_joint_error1 =  Fnorm(J1_c)^2
    ajive_col_ind_error1 = Fnorm(ajive_I1c - I1_c)^2
    ajive_col_error1 = ajive_I1c - signal1
    ajive_col_error1 = Fnorm(ajive_col_error1)^2/Fnorm(signal1)^2
  }
  if (!is.na(ajive_J1c)){
    ajive_col_joint_error1 = Fnorm(ajive_J1c - J1_c)^2
    if (is.na(ajive_I1c)){
      ajive_col_ind_error1 = Fnorm(I1_c)^2
      ajive_col_error1 = ajive_J1c - signal1
      ajive_col_error1 = Fnorm(ajive_col_error1)^2/Fnorm(signal1)^2
    }
    if (!is.na(ajive_I1c)){
      ajive_col_ind_error1 = Fnorm(ajive_I1c - I1_c)^2
      ajive_col_error1 = ajive_J1c + ajive_I1c - signal1
      ajive_col_error1 = Fnorm(ajive_col_error1)^2/Fnorm(signal1)^2
    }
  }
  
  if (is.na(ajive_J2c)){
    ajive_col_joint_error2 =  Fnorm(J2_c)^2
    ajive_col_ind_error2 = Fnorm(ajive_I2c - I2_c)^2
    ajive_col_error2 = ajive_I2c - signal2
    ajive_col_error2 = Fnorm(ajive_col_error2)^2/Fnorm(signal2)^2
  }
  if (!is.na(ajive_J2c)){
    ajive_col_joint_error2 = Fnorm(ajive_J2c - J2_c)^2
    if (is.na(ajive_I2c)){
      ajive_col_ind_error2 = Fnorm(I2_c)^2
      ajive_col_error2 = ajive_J2c - signal2
      ajive_col_error2 = Fnorm(ajive_col_error2)^2/Fnorm(signal2)^2
    }
    if (!is.na(ajive_I2c)){
      ajive_col_ind_error2 = Fnorm(ajive_I2c - I2_c)^2
      ajive_col_error2 = ajive_J2c + ajive_I2c - signal2
      ajive_col_error2 = Fnorm(ajive_col_error2)^2/Fnorm(signal2)^2
    }
  }
  
  list(ajive_row_error1 = ajive_row_error1,
       ajive_row_error2 = ajive_row_error2,
       ajive_row_joint_error1 = ajive_row_joint_error1,
       ajive_row_joint_error2 = ajive_row_joint_error2,
       ajive_row_ind_error1 = ajive_row_ind_error1,
       ajive_row_ind_error2 = ajive_row_ind_error2,
       ajive_col_error1 = ajive_col_error1,
       ajive_col_error2 = ajive_col_error2,
       ajive_col_joint_error1 = ajive_col_joint_error1,
       ajive_col_joint_error2 = ajive_col_joint_error2,
       ajive_col_ind_error1 = ajive_col_ind_error1,
       ajive_col_ind_error2 = ajive_col_ind_error2)
}

stopCluster(cl)
save(output, file = "AJIVE_output.RData")
