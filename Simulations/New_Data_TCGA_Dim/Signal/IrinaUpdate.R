rm(list = ls())
# Check the difference between my original algorithm and DMMD-i.
function_path = "DMMDFunctions/"
source(paste(function_path,"Angle_Calculation.R",sep=''))
source(paste(function_path,"Profile_Likelihood_Rank_Selection.R",sep=''))
source(paste(function_path,"DoubleMatchedMatrixDecomposition.R",sep=''))
source(paste(function_path,"FindOptMatrix.R",sep=''))
source(paste(function_path,"Preliminary_Functions.R",sep=''))
source(paste(function_path,"DMMD_iterative.R",sep=''))

# Get the generated data
load("Simulations/New_Data_TCGA_Dim/Data/Data1.RData")
load("Simulations/New_Data_TCGA_Dim/Data/Data2.RData")

set.seed(37)
n = 88
p = 736

nrep = 100
total_rank1 = rep(8,nrep)
total_rank2 = rep(6,nrep)

joint_rank_col = rep(0,nrep)
joint_rank_row = rep(2,nrep)

Irina_signal_error1 = rep(NA,nrep)
Irina_signal_error2 = rep(NA,nrep)
Irina_joint_row_error1 = rep(NA,nrep)
Irina_joint_row_error2 = rep(NA,nrep)
Irina_ind_row_error1 = rep(NA,nrep)
Irina_ind_row_error2 = rep(NA,nrep)
Irina_joint_col_error1 = rep(NA,nrep)
Irina_joint_col_error2 = rep(NA,nrep)
Irina_ind_col_error1 = rep(NA,nrep)
Irina_ind_col_error2 = rep(NA,nrep)

# Main for loop
for (i in 1:nrep){
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
  Irina_result = DMMD_i(X1, X2, r1 = total_rank1[i], r2 = total_rank2[i], rc = joint_rank_col[i], rr = joint_rank_row[i])
  # The estimation error of signals in DMMD
  Irina_Jc1 = projection(Irina_result$M,ortho = TRUE) %*% Irina_result$A1
  Irina_Jc2 = projection(Irina_result$M,ortho = TRUE) %*% Irina_result$A2
  Irina_Ic1 = Irina_result$A1 - Irina_Jc1
  Irina_Ic2 = Irina_result$A2 - Irina_Jc2
  
  Irina_Jr1 = Irina_result$A1 %*% projection(Irina_result$N,ortho = TRUE)
  Irina_Jr2 = Irina_result$A2 %*% projection(Irina_result$N,ortho = TRUE)
  Irina_Ir1 = Irina_result$A1 - Irina_Jr1
  Irina_Ir2 = Irina_result$A2 - Irina_Jr2

  # Record the error of Irina's update
  Irina_signal_error1[i] = Fnorm(Irina_result$A1 - signal1)^2/Fnorm(signal1)^2
  Irina_signal_error2[i] = Fnorm(Irina_result$A2 - signal2)^2/Fnorm(signal2)^2
  
  Irina_joint_row_error1[i] = Fnorm(Irina_Jr1 - J1_r)^2
  Irina_joint_row_error2[i] = Fnorm(Irina_Jr2 - J2_r)^2
  Irina_ind_row_error1[i] = Fnorm(Irina_Ir1 - I1_r)^2
  Irina_ind_row_error2[i] = Fnorm(Irina_Ir2 - I2_r)^2
  
  Irina_joint_col_error1[i] = Fnorm(Irina_Jc1 - J1_c)^2
  Irina_joint_col_error2[i] = Fnorm(Irina_Jc2 - J2_c)^2
  Irina_ind_col_error1[i] = Fnorm(Irina_Ic1 - I1_c)^2
  Irina_ind_col_error2[i] = Fnorm(Irina_Ic2 - I2_c)^2
}

save(Irina_signal_error1, Irina_signal_error2,
     Irina_joint_row_error1, Irina_joint_row_error2,
     Irina_ind_row_error1, Irina_ind_row_error2,
     Irina_joint_col_error1, Irina_joint_col_error2,
     Irina_ind_col_error1, Irina_ind_col_error2, 
     file = "Simulations/New_Data_TCGA_Dim/Signal/Irina_output.RData")