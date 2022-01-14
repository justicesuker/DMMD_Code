# This is the function that gets rank estimation result from SLIDE.

library(foreach)
library(doParallel)

# Get the generated data
load("Simulations/SimulationData_Setting3/Data.RData")

set.seed(37)

n = 240
p = 200

nrep = 140

cl = makeCluster(28)
registerDoParallel(cl)

# Main for loop

output <- foreach (i = 1:nrep, .errorhandling = 'pass') %dopar% {
  library(SLIDE)
  pvec = c(200,200)
  nvec = c(240,240)
  X1 = X1_list[[i]]
  X2 = X2_list[[i]]
  data_temp_c = cbind(X1,X2)
  data_temp_r = cbind(t(X1),t(X2))
  slide_result_c = slide(data_temp_c, pvec, center = F)$S
  slide_result_r = slide(data_temp_r, nvec, center = F)$S
  list(slide_result_c = slide_result_c,
       slide_result_r = slide_result_r)
}

stopCluster(cl)
save(output,total_rank1,total_rank2,joint_rank_col,joint_rank_row, file = "Simulations/RankEstimation_Setting3_Special_Final/SLIDEoutput.RData")