# This is the function that gets signal identification results from SLIDE.

library(foreach)
library(doParallel)
# Get the generated data
load("Simulations/Data_Setting4_FixRank/Data.RData")
source("Simulations/MyFunction/slide_prelim.R")

set.seed(37)

n = 240
p = 200

nrep = 140
total_rank1 = rep(20,nrep)
total_rank2 = rep(18,nrep)

joint_rank_col = rep(4,nrep)
joint_rank_row = rep(3,nrep)

cl = makeCluster(4)
registerDoParallel(cl)

# Main for loop

output_small <- foreach (i = 1:nrep, .errorhandling = 'pass') %dopar% {
  library(SLIDE)
  r_c = joint_rank_col[i] - 1
  r_r = joint_rank_row[i] - 1
  r1 = total_rank1[i] - 1
  r2 = total_rank2[i] - 1
  pvec = c(200,200)
  nvec = c(240,240)
  X1 = X1_list[[i]]
  X2 = X2_list[[i]]
  data_temp_c = cbind(X1,X2)
  data_temp_r = cbind(t(X1),t(X2))
  S_c = getS(r1, r2, r_c, ndata = 2)
  S_r = getS(r1, r2, r_r, ndata = 2)
  slide_result_c = slide_givenS(data_temp_c, pvec, S_c, standardized = T)
  slide_result_r = slide_givenS(data_temp_r, nvec, S_r, standardized = T)
  list(slide_result_c = slide_result_c, 
       slide_result_r = slide_result_r)
}

output_large <- foreach (i = 1:nrep, .errorhandling = 'pass') %dopar% {
  library(SLIDE)
  r_c = joint_rank_col[i] + 1
  r_r = joint_rank_row[i] + 1
  r1 = total_rank1[i] + 1
  r2 = total_rank2[i] + 1
  pvec = c(200,200)
  nvec = c(240,240)
  X1 = X1_list[[i]]
  X2 = X2_list[[i]]
  data_temp_c = cbind(X1,X2)
  data_temp_r = cbind(t(X1),t(X2))
  S_c = getS(r1, r2, r_c, ndata = 2)
  S_r = getS(r1, r2, r_r, ndata = 2)
  slide_result_c = slide_givenS(data_temp_c, pvec, S_c, standardized = T)
  slide_result_r = slide_givenS(data_temp_r, nvec, S_r, standardized = T)
  list(slide_result_c = slide_result_c, 
       slide_result_r = slide_result_r)
}

output_jointsmall <- foreach (i = 1:nrep, .errorhandling = 'pass') %dopar% {
  library(SLIDE)
  r_c = joint_rank_col[i] - 1
  r_r = joint_rank_row[i] - 1
  r1 = total_rank1[i] 
  r2 = total_rank2[i] 
  pvec = c(200,200)
  nvec = c(240,240)
  X1 = X1_list[[i]]
  X2 = X2_list[[i]]
  data_temp_c = cbind(X1,X2)
  data_temp_r = cbind(t(X1),t(X2))
  S_c = getS(r1, r2, r_c, ndata = 2)
  S_r = getS(r1, r2, r_r, ndata = 2)
  slide_result_c = slide_givenS(data_temp_c, pvec, S_c, standardized = T)
  slide_result_r = slide_givenS(data_temp_r, nvec, S_r, standardized = T)
  list(slide_result_c = slide_result_c, 
       slide_result_r = slide_result_r)
}

stopCluster(cl)

save(output_small, output_large, output_jointsmall, total_rank1,total_rank2,joint_rank_col,joint_rank_row, file = "SLIDEoutput.RData")