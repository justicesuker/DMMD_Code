# This is the function that compares accuracy of signal estimation given true ranks between DMMD and JIVE. 
source("Simulations/MyFunction/Angle_Calculation.R")
source("Simulations/MyFunction/Profile_Likelihood_Rank_Selection.R")
source("Simulations/MyFunction/DoubleMatchedMatrixDecomposition.R")
source("Simulations/MyFunction/FindOptMatrix.R")
source("Simulations/MyFunction/Preliminary_Functions.R")
source("Simulations/MyFunction/Select_ED_Rank.R")

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

cl = makeCluster(28)
registerDoParallel(cl)

angles <- function(Z, Zhat){
  svdZ = svd(Z)
  svdZhat = svd(Zhat)
  rZ = sum(svdZ$d > 1e-6)
  rZhat = sum(svdZhat$d > 1e-6)
  k = min(rZ, rZhat)
  cosines = svd(crossprod(svdZ$u[, 1:rZ], svdZhat$u[, 1:rZhat]))$d[1:k]
  return(list(cosines = cosines, dist = sqrt(sum(1-cosines^2))/sqrt(k)))
}

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
  
  # DMMD
  my_result = DMMD_v2(X1, X2, r1 = total_rank1[i]-1, r2 = total_rank2[i]-1, joint_rank_c = joint_rank_col[i]-1, joint_rank_r = joint_rank_row[i]-1)
  
  # The estimation error of signals in DMMD
  temp_error1 = X1 - my_result$Error$Error1 - signal1
  temp_error2 = X2 - my_result$Error$Error2 - signal2
  
  # Record the relative signal errors
  my_signal_error1 = Fnorm(temp_error1)^2/Fnorm(signal1)^2
  my_signal_error2 = Fnorm(temp_error2)^2/Fnorm(signal2)^2
  
  # Record the chordal distance between total estimated signal row space with true joint row space.
  my_joint_row_error1 = angles(t(X1 - my_result$Error$Error1), t(J1_r))$dist
  my_joint_row_error2 = angles(t(X2 - my_result$Error$Error2), t(J2_r))$dist
  
  # Record the chordal distance between total estimated signal column space with true joint column space.
  my_joint_col_error1 = angles(X1 - my_result$Error$Error1, J1_c)$dist
  my_joint_col_error2 = angles(X2 - my_result$Error$Error2, J2_c)$dist
  
  # JIVE
  library(r.jive)
  data_row = list(X1,X2)
  data_col = list(t(X1),t(X2))
  
  # Run JIVE on the row direction - row-wise decomposition
  temprankA_row = c(total_rank1[i], total_rank2[i]) - joint_rank_row[i]
  jive_result_row = jive(data_row, rankJ = joint_rank_row[i]-1, rankA = temprankA_row, method = 'given', center = FALSE, scale = FALSE)
  
  # Run JIVE on the column direction - column-wise decomposition
  temprankA_col = c(total_rank1[i], total_rank2[i]) - joint_rank_col[i]
  jive_result_col = jive(data_col, rankJ = joint_rank_col[i]-1, rankA = temprankA_col, method = 'given', center = FALSE, scale = FALSE)
  
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
  
  # Record the chordal distance between total estimated signal row space with true joint row space.
  jive_row_joint_error1 = angles(t(jive_result_row$joint[[1]] + jive_result_row$individual[[1]]), t(J1_r))$dist
  jive_row_joint_error2 = angles(t(jive_result_row$joint[[2]] + jive_result_row$individual[[2]]), t(J2_r))$dist
  
  # Record the chordal distance between total estimated signal column space with true joint column space.
  jive_col_joint_error1 = angles(jive_result_row$joint[[1]] + jive_result_row$individual[[1]], J1_c)$dist
  jive_col_joint_error2 = angles(jive_result_row$joint[[2]] + jive_result_row$individual[[2]], J2_c)$dist
  
  list(my_signal_error1 = my_signal_error1,
       my_signal_error2 = my_signal_error2,
       my_joint_row_error1 = my_joint_row_error1,
       my_joint_row_error2 = my_joint_row_error2,
       my_joint_col_error1 = my_joint_col_error1,
       my_joint_col_error2 = my_joint_col_error2,
      
       jive_row_error1 = jive_row_error1,
       jive_row_error2 = jive_row_error2,
       jive_row_joint_error1 = jive_row_joint_error1,
       jive_row_joint_error2 = jive_row_joint_error2,
       jive_col_error1 = jive_col_error1,
       jive_col_error2 = jive_col_error2,
       jive_col_joint_error1 = jive_col_joint_error1,
       jive_col_joint_error2 = jive_col_joint_error2)
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
  
  # DMMD
  my_result = DMMD_v2(X1, X2, r1 = total_rank1[i]+1, r2 = total_rank2[i]+1, joint_rank_c = joint_rank_col[i]+1, joint_rank_r = joint_rank_row[i]+1)
  
  # The estimation error of signals in DMMD
  temp_error1 = X1 - my_result$Error$Error1 - signal1
  temp_error2 = X2 - my_result$Error$Error2 - signal2
  
  # Record the relative signal errors
  my_signal_error1 = Fnorm(temp_error1)^2/Fnorm(signal1)^2
  my_signal_error2 = Fnorm(temp_error2)^2/Fnorm(signal2)^2
  
  # Record the chordal distance between total estimated signal row space with true joint row space.
  my_joint_row_error1 = angles(t(X1 - my_result$Error$Error1), t(J1_r))$dist
  my_joint_row_error2 = angles(t(X2 - my_result$Error$Error2), t(J2_r))$dist
  
  # Record the chordal distance between total estimated signal column space with true joint column space.
  my_joint_col_error1 = angles(X1 - my_result$Error$Error1, J1_c)$dist
  my_joint_col_error2 = angles(X2 - my_result$Error$Error2, J2_c)$dist
  
  # JIVE
  library(r.jive)
  data_row = list(X1,X2)
  data_col = list(t(X1),t(X2))
  
  # Run JIVE on the row direction - row-wise decomposition
  temprankA_row = c(total_rank1[i], total_rank2[i]) - joint_rank_row[i]
  jive_result_row = jive(data_row, rankJ = joint_rank_row[i]+1, rankA = temprankA_row, method = 'given', center = FALSE, scale = FALSE)
  
  # Run JIVE on the column direction - column-wise decomposition
  temprankA_col = c(total_rank1[i], total_rank2[i]) - joint_rank_col[i]
  jive_result_col = jive(data_col, rankJ = joint_rank_col[i]+1, rankA = temprankA_col, method = 'given', center = FALSE, scale = FALSE)
  
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
  
  # Record the chordal distance between total estimated signal row space with true joint row space.
  jive_row_joint_error1 = angles(t(jive_result_row$joint[[1]] + jive_result_row$individual[[1]]), t(J1_r))$dist
  jive_row_joint_error2 = angles(t(jive_result_row$joint[[2]] + jive_result_row$individual[[2]]), t(J2_r))$dist
  
  # Record the chordal distance between total estimated signal column space with true joint column space.
  jive_col_joint_error1 = angles(jive_result_row$joint[[1]] + jive_result_row$individual[[1]], J1_c)$dist
  jive_col_joint_error2 = angles(jive_result_row$joint[[2]] + jive_result_row$individual[[2]], J2_c)$dist
  
  list(my_signal_error1 = my_signal_error1,
       my_signal_error2 = my_signal_error2,
       my_joint_row_error1 = my_joint_row_error1,
       my_joint_row_error2 = my_joint_row_error2,
       my_joint_col_error1 = my_joint_col_error1,
       my_joint_col_error2 = my_joint_col_error2,
       
       jive_row_error1 = jive_row_error1,
       jive_row_error2 = jive_row_error2,
       jive_row_joint_error1 = jive_row_joint_error1,
       jive_row_joint_error2 = jive_row_joint_error2,
       jive_col_error1 = jive_col_error1,
       jive_col_error2 = jive_col_error2,
       jive_col_joint_error1 = jive_col_joint_error1,
       jive_col_joint_error2 = jive_col_joint_error2)
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
  
  # DMMD
  my_result = DMMD_v2(X1, X2, r1 = total_rank1[i], r2 = total_rank2[i], joint_rank_c = joint_rank_col[i]-1, joint_rank_r = joint_rank_row[i]-1)
  
  # The estimation error of signals in DMMD
  temp_error1 = X1 - my_result$Error$Error1 - signal1
  temp_error2 = X2 - my_result$Error$Error2 - signal2
  
  # Record the relative signal errors
  my_signal_error1 = Fnorm(temp_error1)^2/Fnorm(signal1)^2
  my_signal_error2 = Fnorm(temp_error2)^2/Fnorm(signal2)^2
  
  # Record the chordal distance between total estimated signal row space with true joint row space.
  my_joint_row_error1 = angles(t(X1 - my_result$Error$Error1), t(J1_r))$dist
  my_joint_row_error2 = angles(t(X2 - my_result$Error$Error2), t(J2_r))$dist
  
  # Record the chordal distance between total estimated signal column space with true joint column space.
  my_joint_col_error1 = angles(X1 - my_result$Error$Error1, J1_c)$dist
  my_joint_col_error2 = angles(X2 - my_result$Error$Error2, J2_c)$dist
  
  # JIVE
  library(r.jive)
  data_row = list(X1,X2)
  data_col = list(t(X1),t(X2))
  
  # Run JIVE on the row direction - row-wise decomposition
  temprankA_row = c(total_rank1[i], total_rank2[i]) - (joint_rank_row[i] - 1)
  jive_result_row = jive(data_row, rankJ = joint_rank_row[i]-1, rankA = temprankA_row, method = 'given', center = FALSE, scale = FALSE)
  
  # Run JIVE on the column direction - column-wise decomposition
  temprankA_col = c(total_rank1[i], total_rank2[i]) - (joint_rank_col[i] - 1)
  jive_result_col = jive(data_col, rankJ = joint_rank_col[i]-1, rankA = temprankA_col, method = 'given', center = FALSE, scale = FALSE)
  
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
  
  # Record the chordal distance between total estimated signal row space with true joint row space.
  jive_row_joint_error1 = angles(t(jive_result_row$joint[[1]] + jive_result_row$individual[[1]]), t(J1_r))$dist
  jive_row_joint_error2 = angles(t(jive_result_row$joint[[2]] + jive_result_row$individual[[2]]), t(J2_r))$dist
  
  # Record the chordal distance between total estimated signal column space with true joint column space.
  jive_col_joint_error1 = angles(jive_result_row$joint[[1]] + jive_result_row$individual[[1]], J1_c)$dist
  jive_col_joint_error2 = angles(jive_result_row$joint[[2]] + jive_result_row$individual[[2]], J2_c)$dist
  
  list(my_signal_error1 = my_signal_error1,
       my_signal_error2 = my_signal_error2,
       my_joint_row_error1 = my_joint_row_error1,
       my_joint_row_error2 = my_joint_row_error2,
       my_joint_col_error1 = my_joint_col_error1,
       my_joint_col_error2 = my_joint_col_error2,
       
       jive_row_error1 = jive_row_error1,
       jive_row_error2 = jive_row_error2,
       jive_row_joint_error1 = jive_row_joint_error1,
       jive_row_joint_error2 = jive_row_joint_error2,
       jive_col_error1 = jive_col_error1,
       jive_col_error2 = jive_col_error2,
       jive_col_joint_error1 = jive_col_joint_error1,
       jive_col_joint_error2 = jive_col_joint_error2)
}

stopCluster(cl)
save(output_small, output_large, output_jointsmall, file = "output.RData")
