# This is the function that gets signal identification results from SLIDE.

library(foreach)
library(doParallel)

# Get the generated data
load("Simulations/New_Data_TCGA_Dim/Data/Data1.RData")
load("Simulations/New_Data_TCGA_Dim/Data/Data2.RData")
source("OtherFunctions/slide_prelim.R")

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
  library(SLIDE)
  r_c = joint_rank_col[i]
  r_r = joint_rank_row[i]
  r1 = total_rank1[i]
  r2 = total_rank2[i]
  pvec = c(p,p)
  nvec = c(n,n)
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

save(output,total_rank1,total_rank2,joint_rank_col,joint_rank_row, file = "Simulations/New_Data_TCGA_Dim/Signal/SLIDEoutput.RData")