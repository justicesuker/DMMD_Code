# This is the function that compares accuracy of signal estimation given true ranks between DMMD and AJIVE. 
rm(list=ls())
function_path = "DMMDFunctions/"
source(paste(function_path,"Angle_Calculation.R",sep=''))
source(paste(function_path,"Profile_Likelihood_Rank_Selection.R",sep=''))
source(paste(function_path,"DoubleMatchedMatrixDecomposition.R",sep=''))
source(paste(function_path,"FindOptMatrix.R",sep=''))
source(paste(function_path,"Preliminary_Functions.R",sep=''))
library(foreach)
library(doParallel)

# Get the generated data
load("Simulations/SimulationData_Setting6/Data1.RData")
load("Simulations/SimulationData_Setting6/Data2.RData")

set.seed(37)
n = 240
p = 200
nrep = 140
total_rank1 = rep(20,nrep)
total_rank2 = rep(18,nrep)

min_total = rep(18,nrep)

joint_rank_col = c(rep(0, nrep/2), min_total[(nrep/2+1):nrep])
joint_rank_row = c(rep(0, nrep/4), min_total[(nrep/4+1):(nrep/2)], rep(0, nrep/4), min_total[(3*nrep/4+1):nrep])

cl = makeCluster(4)
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
save(output, file = "Simulations/Signal_Identification_Setting6_Special/AJIVE_output.RData")
