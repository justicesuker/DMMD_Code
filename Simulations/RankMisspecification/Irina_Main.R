# This is the function that compares accuracy of signal estimation given true ranks between DMMD and JIVE. 
source("../MyFunction/Angle_Calculation.R")
source("../MyFunction/Profile_Likelihood_Rank_Selection.R")
source("../MyFunction/DoubleMatchedMatrixDecomposition.R")
source("../MyFunction/FindOptMatrix.R")
source("../MyFunction/Preliminary_Functions.R")
source("../MyFunction/Select_ED_Rank.R")
source("../IrinaFunction/DMMD_Irina.R")

library(foreach)
library(doParallel)

# Get the generated data
load("../Data_Setting4_FixRank/Data.RData")

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
output_small <- foreach (i = 1:nrep) %dopar% {
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
  
  # DMMD-i
  Irina_result = DMMD_Irina(X1, X2, r1 = total_rank1[i]-1, r2 = total_rank2[i]-1, rc = joint_rank_col[i]-1, rr = joint_rank_row[i]-1)
  
  Irina_Jc1 = projection(Irina_result$M,ortho = TRUE) %*% Irina_result$A1
  Irina_Jc2 = projection(Irina_result$M,ortho = TRUE) %*% Irina_result$A2
  Irina_Ic1 = Irina_result$A1 - Irina_Jc1
  Irina_Ic2 = Irina_result$A2 - Irina_Jc2
  
  Irina_Jr1 = Irina_result$A1 %*% projection(Irina_result$N,ortho = TRUE)
  Irina_Jr2 = Irina_result$A2 %*% projection(Irina_result$N,ortho = TRUE)
  Irina_Ir1 = Irina_result$A1 - Irina_Jr1
  Irina_Ir2 = Irina_result$A2 - Irina_Jr2
  
  # Record the error of Irina's update
  Irina_signal_error1 = Fnorm(Irina_result$A1 - signal1)^2/Fnorm(signal1)^2
  Irina_signal_error2 = Fnorm(Irina_result$A2 - signal2)^2/Fnorm(signal2)^2
  
  Irina_joint_row_error1 = angles(t(Irina_result$A1), t(J1_r))$dist
  Irina_joint_row_error2 = angles(t(Irina_result$A2), t(J2_r))$dist

  Irina_joint_col_error1 = angles(Irina_result$A1, J1_c)$dist
  Irina_joint_col_error2 = angles(Irina_result$A2, J2_c)$dist
  
  list(Irina_signal_error1 = Irina_signal_error1,
       Irina_signal_error2 = Irina_signal_error2,
       Irina_joint_row_error1 = Irina_joint_row_error1,
       Irina_joint_row_error2 = Irina_joint_row_error2,
       Irina_joint_col_error1 = Irina_joint_col_error1,
       Irina_joint_col_error2 = Irina_joint_col_error2)
}

output_large <- foreach (i = 1:nrep) %dopar% {
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
  
  # DMMD-i
  Irina_result = DMMD_Irina(X1, X2, r1 = total_rank1[i]+1, r2 = total_rank2[i]+1, rc = joint_rank_col[i]+1, rr = joint_rank_row[i]+1)
  
  Irina_Jc1 = projection(Irina_result$M,ortho = TRUE) %*% Irina_result$A1
  Irina_Jc2 = projection(Irina_result$M,ortho = TRUE) %*% Irina_result$A2
  Irina_Ic1 = Irina_result$A1 - Irina_Jc1
  Irina_Ic2 = Irina_result$A2 - Irina_Jc2
  
  Irina_Jr1 = Irina_result$A1 %*% projection(Irina_result$N,ortho = TRUE)
  Irina_Jr2 = Irina_result$A2 %*% projection(Irina_result$N,ortho = TRUE)
  Irina_Ir1 = Irina_result$A1 - Irina_Jr1
  Irina_Ir2 = Irina_result$A2 - Irina_Jr2
  
  # Record the error of Irina's update
  Irina_signal_error1 = Fnorm(Irina_result$A1 - signal1)^2/Fnorm(signal1)^2
  Irina_signal_error2 = Fnorm(Irina_result$A2 - signal2)^2/Fnorm(signal2)^2
  
  Irina_joint_row_error1 = angles(t(Irina_result$A1), t(J1_r))$dist
  Irina_joint_row_error2 = angles(t(Irina_result$A2), t(J2_r))$dist
  
  Irina_joint_col_error1 = angles(Irina_result$A1, J1_c)$dist
  Irina_joint_col_error2 = angles(Irina_result$A2, J2_c)$dist
  
  list(Irina_signal_error1 = Irina_signal_error1,
       Irina_signal_error2 = Irina_signal_error2,
       Irina_joint_row_error1 = Irina_joint_row_error1,
       Irina_joint_row_error2 = Irina_joint_row_error2,
       Irina_joint_col_error1 = Irina_joint_col_error1,
       Irina_joint_col_error2 = Irina_joint_col_error2)
}

output_jointsmall <- foreach (i = 1:nrep) %dopar% {
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
  
  # DMMD-i
  Irina_result = DMMD_Irina(X1, X2, r1 = total_rank1[i], r2 = total_rank2[i], rc = joint_rank_col[i]-1, rr = joint_rank_row[i]-1)
  
  Irina_Jc1 = projection(Irina_result$M,ortho = TRUE) %*% Irina_result$A1
  Irina_Jc2 = projection(Irina_result$M,ortho = TRUE) %*% Irina_result$A2
  Irina_Ic1 = Irina_result$A1 - Irina_Jc1
  Irina_Ic2 = Irina_result$A2 - Irina_Jc2
  
  Irina_Jr1 = Irina_result$A1 %*% projection(Irina_result$N,ortho = TRUE)
  Irina_Jr2 = Irina_result$A2 %*% projection(Irina_result$N,ortho = TRUE)
  Irina_Ir1 = Irina_result$A1 - Irina_Jr1
  Irina_Ir2 = Irina_result$A2 - Irina_Jr2
  
  # Record the error of Irina's update
  Irina_signal_error1 = Fnorm(Irina_result$A1 - signal1)^2/Fnorm(signal1)^2
  Irina_signal_error2 = Fnorm(Irina_result$A2 - signal2)^2/Fnorm(signal2)^2
  
  Irina_joint_row_error1 = angles(t(Irina_result$A1), t(J1_r))$dist
  Irina_joint_row_error2 = angles(t(Irina_result$A2), t(J2_r))$dist
  
  Irina_joint_col_error1 = angles(Irina_result$A1, J1_c)$dist
  Irina_joint_col_error2 = angles(Irina_result$A2, J2_c)$dist
  
  list(Irina_signal_error1 = Irina_signal_error1,
       Irina_signal_error2 = Irina_signal_error2,
       Irina_joint_row_error1 = Irina_joint_row_error1,
       Irina_joint_row_error2 = Irina_joint_row_error2,
       Irina_joint_col_error1 = Irina_joint_col_error1,
       Irina_joint_col_error2 = Irina_joint_col_error2)
}

stopCluster(cl)
save(output_small, output_large, output_jointsmall, file = "Irina_output.RData")