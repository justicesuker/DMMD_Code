# This is the function that compares accuracy of rank estimation.
# The methods we compare are profile likelihood (PL), edge distribution (ED) and permutation method used in JIVE
function_path1 = "DMMDFunctions/"
function_path2 = "OtherFunctions/"
source(paste(function_path1,"Angle_Calculation.R",sep=''))
source(paste(function_path1,"Profile_Likelihood_Rank_Selection.R",sep=''))
source(paste(function_path1,"DoubleMatchedMatrixDecomposition.R",sep=''))
source(paste(function_path1,"FindOptMatrix.R",sep=''))
source(paste(function_path1,"Preliminary_Functions.R",sep=''))
source(paste(function_path2,"Select_ED_Rank.R",sep=''))

library(foreach)
library(doParallel)

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

cl = makeCluster(28)
registerDoParallel(cl)

# Main for loop

output <- foreach (i = 1:nrep) %dopar% {

  X1 = X1_list[[i]]
  X2 = X2_list[[i]]
  signal1 = signal1_list[[i]]
  signal2 = signal2_list[[i]]
  
  # PL method
  my_result_PL = DMMD_v2(X1,X2)
  # Estimation error of joint rank
  my_error_joint_rank_col_PL = my_result_PL$Rank$"Joint Column Rank" - joint_rank_col[i]
  my_error_joint_rank_row_PL = my_result_PL$Rank$"Joint Row Rank" - joint_rank_row[i]
  # Estimation error of total rank
  my_error_total_rank1_PL = my_result_PL$Rank$"Rank 1" - total_rank1[i]
  my_error_total_rank2_PL = my_result_PL$Rank$"Rank 2" - total_rank2[i]

  # ED method
  my_result_ED = DMMD_v2(X1,X2, method = "ED")
  # Estimation error of joint rank
  my_error_joint_rank_col_ED = my_result_ED$Rank$"Joint Column Rank" - joint_rank_col[i]
  my_error_joint_rank_row_ED = my_result_ED$Rank$"Joint Row Rank" - joint_rank_row[i]
  # Estimation error of total rank
  my_error_total_rank1_ED = my_result_ED$Rank$"Rank 1" - total_rank1[i]
  my_error_total_rank2_ED = my_result_ED$Rank$"Rank 2" - total_rank2[i]

  library(r.jive)
  # Run JIVE on the column direction
  data_temp = list(t(X1),t(X2))
  jive_result_c = jive(data_temp,scale = FALSE, center = FALSE)
  # Estimation error of joint rank
  jive_error_joint_rank_c = jive_result_c$rankJ - joint_rank_col[i]
  # Estimation error of total rank
  jive_error_total_rank1_c = jive_result_c$rankJ + jive_result_c$rankA[1] - total_rank1[i]
  jive_error_total_rank2_c = jive_result_c$rankJ + jive_result_c$rankA[2] - total_rank2[i]

  # Run JIVE on the row direction
  data_temp2 = list(X1,X2)
  jive_result_r = jive(data_temp2,scale = FALSE, center = FALSE)
  # Estimation error of joint rank
  jive_error_joint_rank_r = jive_result_r$rankJ - joint_rank_row[i]
  # Estimation error of total rank
  jive_error_total_rank1_r = jive_result_r$rankJ + jive_result_r$rankA[1] - total_rank1[i]
  jive_error_total_rank2_r = jive_result_r$rankJ + jive_result_r$rankA[2] - total_rank2[i]

  list(my_error_joint_rank_col_PL = my_error_joint_rank_col_PL,
       my_error_joint_rank_row_PL = my_error_joint_rank_row_PL,
       my_error_total_rank1_PL = my_error_total_rank1_PL,
       my_error_total_rank2_PL = my_error_total_rank2_PL,

       my_error_joint_rank_col_ED = my_error_joint_rank_col_ED,
       my_error_joint_rank_row_ED = my_error_joint_rank_row_ED,
       my_error_total_rank1_ED = my_error_total_rank1_ED,
       my_error_total_rank2_ED = my_error_total_rank2_ED,

       jive_error_joint_rank_c = jive_error_joint_rank_c,
       jive_error_total_rank1_c = jive_error_total_rank1_c,
       jive_error_total_rank2_c = jive_error_total_rank2_c,
       jive_error_joint_rank_r = jive_error_joint_rank_r,
       jive_error_total_rank1_r = jive_error_total_rank1_r,
       jive_error_total_rank2_r = jive_error_total_rank2_r)

}

stopCluster(cl)

save(output,total_rank1,total_rank2,joint_rank_col,joint_rank_row, file = "Simulations/New_Data_TCGA_Dim/Rank/output.RData")