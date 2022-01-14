# Check the difference between my original algorithm and Irina's update.
# This is the function that compares accuracy of signal estimation given true ranks between DMMD and JIVE. 
source("../../../MyFunction/Angle_Calculation.R")
source("../../../MyFunction/Profile_Likelihood_Rank_Selection.R")
source("../../../MyFunction/DoubleMatchedMatrixDecomposition.R")
source("../../../MyFunction/FindOptMatrix.R")
source("../../../MyFunction/Preliminary_Functions.R")
source("../../../MyFunction/Select_ED_Rank.R")
source("../../../IrinaFunction/DMMD_Irina.R")

library(foreach)
library(doParallel)

# Get the generated data
load("../../../SimulationData_Setting6/Data1.RData")
load("../../../SimulationData_Setting6/Data2.RData")

set.seed(37)
n = 240
p = 200

# Total rank 20, 18
# Joint column rank 4 
# Joint row rank 3

nrep = 140
total_rank1 = rep(20,nrep)
total_rank2 = rep(18,nrep)

min_total = rep(18,nrep)

joint_rank_col = c(rep(0, nrep/2), min_total[(nrep/2+1):nrep])
joint_rank_row = c(rep(0, nrep/4), min_total[(nrep/4+1):(nrep/2)], rep(0, nrep/4), min_total[(3*nrep/4+1):nrep])

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
  my_joint_row_error1 = Fnorm(my_result$`Row Decomposition`$`Joint Row 1` - J1_r)^2
  my_joint_row_error2 = Fnorm(my_result$`Row Decomposition`$`Joint Row 2` - J2_r)^2
  my_ind_row_error1 = Fnorm(my_result$`Row Decomposition`$`Individual Row 1` - I1_r)^2
  my_ind_row_error2 = Fnorm(my_result$`Row Decomposition`$`Individual Row 2` - I2_r)^2
  
  # Record the errors of column decomposition in DMMD
  my_joint_col_error1 = Fnorm(my_result$`Column Decomposition`$`Joint Column 1` - J1_c)^2
  my_joint_col_error2 = Fnorm(my_result$`Column Decomposition`$`Joint Column 2` - J2_c)^2
  my_ind_col_error1 = Fnorm(my_result$`Column Decomposition`$`Individual Column 1` - I1_c)^2
  my_ind_col_error2 = Fnorm(my_result$`Column Decomposition`$`Individual Column 2` - I2_c)^2
  
  # Irina update
  Irina_result = DMMD_Irina(X1, X2, r1 = total_rank1[i], r2 = total_rank2[i], rc = joint_rank_col[i], rr = joint_rank_row[i])
  Alliter_result = DMMD_All(X1, X2, r1 = total_rank1[i], r2 = total_rank2[i], rc = joint_rank_col[i], rr = joint_rank_row[i])
  
  # Get all decomposition sections from Irina update
  Irina_Jc1 = projection(Irina_result$M,ortho = TRUE) %*% Irina_result$A1
  Irina_Jc2 = projection(Irina_result$M,ortho = TRUE) %*% Irina_result$A2
  Irina_Ic1 = Irina_result$A1 - Irina_Jc1
  Irina_Ic2 = Irina_result$A2 - Irina_Jc2
  
  Irina_Jr1 = Irina_result$A1 %*% projection(Irina_result$N,ortho = TRUE)
  Irina_Jr2 = Irina_result$A2 %*% projection(Irina_result$N,ortho = TRUE)
  Irina_Ir1 = Irina_result$A1 - Irina_Jr1
  Irina_Ir2 = Irina_result$A2 - Irina_Jr2
  
  Alliter_Jc1 = projection(Alliter_result$M,ortho = TRUE) %*% Alliter_result$A1
  Alliter_Jc2 = projection(Alliter_result$M,ortho = TRUE) %*% Alliter_result$A2
  Alliter_Ic1 = Alliter_result$A1 - Alliter_Jc1
  Alliter_Ic2 = Alliter_result$A2 - Alliter_Jc2
  
  Alliter_Jr1 = Alliter_result$A1 %*% projection(Alliter_result$N,ortho = TRUE)
  Alliter_Jr2 = Alliter_result$A2 %*% projection(Alliter_result$N,ortho = TRUE)
  Alliter_Ir1 = Alliter_result$A1 - Alliter_Jr1
  Alliter_Ir2 = Alliter_result$A2 - Alliter_Jr2
  
  # Record the error of Irina's update
  Irina_signal_error1 = Fnorm(Irina_result$A1 - signal1)^2/Fnorm(signal1)^2
  Irina_signal_error2 = Fnorm(Irina_result$A2 - signal2)^2/Fnorm(signal2)^2
  
  Irina_joint_row_error1 = Fnorm(Irina_Jr1 - J1_r)^2
  Irina_joint_row_error2 = Fnorm(Irina_Jr2 - J2_r)^2
  Irina_ind_row_error1 = Fnorm(Irina_Ir1 - I1_r)^2
  Irina_ind_row_error2 = Fnorm(Irina_Ir2 - I2_r)^2
  
  Irina_joint_col_error1 = Fnorm(Irina_Jc1 - J1_c)^2
  Irina_joint_col_error2 = Fnorm(Irina_Jc2 - J2_c)^2
  Irina_ind_col_error1 = Fnorm(Irina_Ic1 - I1_c)^2
  Irina_ind_col_error2 = Fnorm(Irina_Ic2 - I2_c)^2
  
  # Record the error of updating all the components.
  Alliter_signal_error1 = Fnorm(Alliter_result$A1 - signal1)^2/Fnorm(signal1)^2
  Alliter_signal_error2 = Fnorm(Alliter_result$A2 - signal2)^2/Fnorm(signal2)^2
  
  Alliter_joint_row_error1 = Fnorm(Alliter_Jr1 - J1_r)^2
  Alliter_joint_row_error2 = Fnorm(Alliter_Jr2 - J2_r)^2
  Alliter_ind_row_error1 = Fnorm(Alliter_Ir1 - I1_r)^2
  Alliter_ind_row_error2 = Fnorm(Alliter_Ir2 - I2_r)^2
  
  Alliter_joint_col_error1 = Fnorm(Alliter_Jc1 - J1_c)^2
  Alliter_joint_col_error2 = Fnorm(Alliter_Jc2 - J2_c)^2
  Alliter_ind_col_error1 = Fnorm(Alliter_Ic1 - I1_c)^2
  Alliter_ind_col_error2 = Fnorm(Alliter_Ic2 - I2_c)^2
  
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
       
       Irina_signal_error1 = Irina_signal_error1,
       Irina_signal_error2 = Irina_signal_error2,
       Irina_joint_row_error1 = Irina_joint_row_error1,
       Irina_joint_row_error2 = Irina_joint_row_error2,
       Irina_ind_row_error1 = Irina_ind_row_error1,
       Irina_ind_row_error2 = Irina_ind_row_error2,
       Irina_joint_col_error1 = Irina_joint_col_error1,
       Irina_joint_col_error2 = Irina_joint_col_error2,
       Irina_ind_col_error1 = Irina_ind_col_error1,
       Irina_ind_col_error2 = Irina_ind_col_error2,
       
       Alliter_signal_error1 = Alliter_signal_error1,
       Alliter_signal_error2 = Alliter_signal_error2,
       Alliter_joint_row_error1 = Alliter_joint_row_error1,
       Alliter_joint_row_error2 = Alliter_joint_row_error2,
       Alliter_ind_row_error1 = Alliter_ind_row_error1,
       Alliter_ind_row_error2 = Alliter_ind_row_error2,
       Alliter_joint_col_error1 = Alliter_joint_col_error1,
       Alliter_joint_col_error2 = Alliter_joint_col_error2,
       Alliter_ind_col_error1 = Alliter_ind_col_error1,
       Alliter_ind_col_error2 = Alliter_ind_col_error2)
}

stopCluster(cl)
save(output, file = "output.RData")