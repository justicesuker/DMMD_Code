# This is the function that compares accuracy of signal estimation given true ranks between DMMD and JIVE. 
function_path = "DMMDFunctions/"
source(paste(function_path,"Angle_Calculation.R",sep=''))
source(paste(function_path,"Profile_Likelihood_Rank_Selection.R",sep=''))
source(paste(function_path,"DoubleMatchedMatrixDecomposition.R",sep=''))
source(paste(function_path,"FindOptMatrix.R",sep=''))
source(paste(function_path,"Preliminary_Functions.R",sep=''))

library(foreach)
library(doParallel)

# Get the generated data
load("Simulations/Data_Setting5_FixRank_SNR0.5/Data.RData")

set.seed(37)
n = 240
p = 200

# Total rank 20, 18
# Joint column rank 4 
# Joint row rank 3

nrep = 140
total_rank1 = rep(20,nrep)
total_rank2 = rep(18,nrep)

joint_rank_col = rep(4,nrep)
joint_rank_row = rep(3,nrep)

cl = makeCluster(28)
registerDoParallel(cl)

# Main for loop
output <- foreach (i = 1:nrep) %dopar% {
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
  
  # DMMD
  my_result = DMMD_v2(X1, X2, r1 = total_rank1[i], r2 = total_rank2[i], joint_rank_c = joint_rank_col[i], joint_rank_r = joint_rank_row[i])
  
  # The estimation error of signals in DMMD
  temp_error1 = X1 - my_result$Error$Error1 - signal1
  temp_error2 = X2 - my_result$Error$Error2 - signal2
  
  # Record the relative signal errors
  my_signal_error1 = Fnorm(temp_error1)^2/Fnorm(signal1)^2
  my_signal_error2 = Fnorm(temp_error2)^2/Fnorm(signal2)^2
  
  # Record the errors of row decomposition in DMMD
  my_joint_row_error1 = Fnorm(my_result$`Row Decomposition`$`Joint Row 1` - J1_r)^2/Fnorm(J1_r)^2
  my_joint_row_error2 = Fnorm(my_result$`Row Decomposition`$`Joint Row 2` - J2_r)^2/Fnorm(J2_r)^2
  my_ind_row_error1 = Fnorm(my_result$`Row Decomposition`$`Individual Row 1` - I1_r)^2/Fnorm(I1_r)^2
  my_ind_row_error2 = Fnorm(my_result$`Row Decomposition`$`Individual Row 2` - I2_r)^2/Fnorm(I2_r)^2
  
  # Record the errors of column decomposition in DMMD
  my_joint_col_error1 = Fnorm(my_result$`Column Decomposition`$`Joint Column 1` - J1_c)^2/Fnorm(J1_c)^2
  my_joint_col_error2 = Fnorm(my_result$`Column Decomposition`$`Joint Column 2` - J2_c)^2/Fnorm(J2_c)^2
  my_ind_col_error1 = Fnorm(my_result$`Column Decomposition`$`Individual Column 1` - I1_c)^2/Fnorm(I1_c)^2
  my_ind_col_error2 = Fnorm(my_result$`Column Decomposition`$`Individual Column 2` - I2_c)^2/Fnorm(I2_c)^2
  
  # JIVE
  library(r.jive)
  data_row = list(X1,X2)
  data_col = list(t(X1),t(X2))
  
  # Run JIVE on the row direction - row-wise decomposition
  temprankA_row = c(total_rank1[i], total_rank2[i]) - joint_rank_row[i]
  jive_result_row = jive(data_row, rankJ = joint_rank_row[i], rankA = temprankA_row, method = 'given', center = FALSE, scale = FALSE)
  
  # Run JIVE on the column direction - column-wise decomposition
  temprankA_col = c(total_rank1[i], total_rank2[i]) - joint_rank_col[i]
  jive_result_col = jive(data_col, rankJ = joint_rank_col[i], rankA = temprankA_col, method = 'given', center = FALSE, scale = FALSE)
  
  # The estimation error of signals in row-wise JIVE 
  jive_row_error1 = jive_result_row$joint[[1]] + jive_result_row$individual[[1]] - signal1
  jive_row_error2 = jive_result_row$joint[[2]] + jive_result_row$individual[[2]] - signal2
  
  # Relative error
  jive_row_error1 = Fnorm(jive_row_error1)^2/Fnorm(signal1)^2
  jive_row_error2 = Fnorm(jive_row_error2)^2/Fnorm(signal2)^2
  
  # The estimation error of signals in column-wise JIVE 
  jive_col_error1 = jive_result_col$joint[[1]] + jive_result_col$individual[[1]] - t(signal1)
  jive_col_error2 = jive_result_col$joint[[2]] + jive_result_col$individual[[2]] - t(signal2)
  
  # Relative error
  jive_col_error1 = Fnorm(jive_col_error1)^2/Fnorm(signal1)^2
  jive_col_error2 = Fnorm(jive_col_error2)^2/Fnorm(signal2)^2
  
  # Record the errors of row decomposition in JIVE
  jive_row_joint_error1 = Fnorm(jive_result_row$joint[[1]] - J1_r)^2/Fnorm(J1_r)^2
  jive_row_joint_error2 = Fnorm(jive_result_row$joint[[2]] - J2_r)^2/Fnorm(J2_r)^2
  jive_row_ind_error1 = Fnorm(jive_result_row$individual[[1]] - I1_r)^2/Fnorm(I1_r)^2
  jive_row_ind_error2 = Fnorm(jive_result_row$individual[[2]] - I2_r)^2/Fnorm(I2_r)^2
  
  # Record the errors of column decomposition in JIVE
  jive_col_joint_error1 = Fnorm(t(jive_result_col$joint[[1]]) - J1_c)^2/Fnorm(J1_c)^2
  jive_col_joint_error2 = Fnorm(t(jive_result_col$joint[[2]]) - J2_c)^2/Fnorm(J2_c)^2
  jive_col_ind_error1 = Fnorm(t(jive_result_col$individual[[1]]) - I1_c)^2/Fnorm(I1_c)^2
  jive_col_ind_error2 = Fnorm(t(jive_result_col$individual[[2]]) - I2_c)^2/Fnorm(I2_c)^2
  
  list(my_signal_error1 = my_signal_error1,
       my_signal_error2 = my_signal_error2,
       my_joint_row_error1 = my_joint_row_error1,
       my_joint_row_error2 = my_joint_row_error2,
       my_ind_row_error1 = my_ind_row_error1,
       my_ind_row_error2 = my_ind_row_error2,
       my_joint_col_error1 = my_joint_col_error1,
       my_joint_col_error2 = my_joint_col_error2,
       my_ind_col_error1 = my_ind_col_error1,
       my_ind_col_error2 = my_ind_col_error2,
      
       jive_row_error1 = jive_row_error1,
       jive_row_error2 = jive_row_error2,
       jive_row_joint_error1 = jive_row_joint_error1,
       jive_row_joint_error2 = jive_row_joint_error2,
       jive_row_ind_error1 = jive_row_ind_error1,
       jive_row_ind_error2 = jive_row_ind_error2,
       jive_col_error1 = jive_col_error1,
       jive_col_error2 = jive_col_error2,
       jive_col_joint_error1 = jive_col_joint_error1,
       jive_col_joint_error2 = jive_col_joint_error2,
       jive_col_ind_error1 = jive_col_ind_error1,
       jive_col_ind_error2 = jive_col_ind_error2)
}

stopCluster(cl)
save(output, file = "Simulations/Signal_Identification_Setting5/output.RData")
