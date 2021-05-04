# This is the function that gets signal identification results from SLIDE.

library(foreach)
library(doParallel)

# Get the generated data
load("Simulations/Data_Setting5_FixRank_SNR0.5/Data.RData")

set.seed(37)

n = 240
p = 200

nrep = 140
total_rank1 = rep(20,nrep)
total_rank2 = rep(18,nrep)

joint_rank_col = rep(4,nrep)
joint_rank_row = rep(3,nrep)

cl = makeCluster(28)
registerDoParallel(cl)

# Main for loop

output <- foreach (i = 1:nrep, .errorhandling = 'pass') %dopar% {
  library(SLIDE)
  r_c = joint_rank_col[i]
  r_r = joint_rank_row[i]
  r1 = total_rank1[i]
  r2 = total_rank2[i]
  pvec = c(200,200)
  nvec = c(240,240)
  X1 = X1_list[[i]]
  X2 = X2_list[[i]]
  data_temp_c = cbind(X1,X2)
  data_temp_r = cbind(t(X1),t(X2))
  S_c = matrix(rep(0,2*(r1+r2-r_c)),nrow = 2)
  if (r_c != 0){
    for (i in 1:r_c){
      S_c[,i] = c(1,1)
    }
    if (r1-r_c>0){
      for (i in 1:(r1-r_c)){
        S_c[,i+r_c] = c(1,0)
      }
      if (r2-r_c>0){
        for (i in 1:(r2-r_c)){
          S_c[,i+r1] = c(0,1)
        }
      }
    }
    if(r1-r_c==0){
      if (r2-r_c>0){
        for (i in 1:(r2-r_c)){
          S_c[,i+r1] = c(0,1)
        }
      }
    }
  }
  if (r_c == 0){
    if (r1-r_c>0){
      for (i in 1:(r1-r_c)){
        S_c[,i+r_c] = c(1,0)
      }
    }
    if (r1-r_c==0){
      if (r2-r_c>0){
        for (i in 1:(r2-r_c)){
          S_c[,i+r1] = c(0,1)
        }
      }
    }
  }
  S_r = matrix(rep(0,2*(r1+r2-r_r)),nrow = 2)
  if (r_r != 0){
    for (i in 1:r_r){
      S_r[,i] = c(1,1)
    }
    if (r1-r_r>0){
      for (i in 1:(r1-r_r)){
        S_r[,i+r_r] = c(1,0)
      }
      if (r2-r_r>0){
        for (i in 1:(r2-r_r)){
          S_r[,i+r1] = c(0,1)
        }
      }
    }
    if(r1-r_r==0){
      if (r2-r_r>0){
        for (i in 1:(r2-r_r)){
          S_r[,i+r1] = c(0,1)
        }
      }
    }
  }
  if (r_r == 0){
    if (r1-r_r>0){
      for (i in 1:(r1-r_r)){
        S_r[,i+r_r] = c(1,0)
      }
    }
    if (r1-r_r==0){
      if (r2-r_r>0){
        for (i in 1:(r2-r_r)){
          S_r[,i+r1] = c(0,1)
        }
      }
    }
  }
  slide_result_c = slide_givenS(data_temp_c, pvec, S_c, standardized = T)
  slide_result_r = slide_givenS(data_temp_r, nvec, S_r, standardized = T)
  list(slide_result_c = slide_result_c, 
       slide_result_r = slide_result_r)
}

stopCluster(cl)

save(output,total_rank1,total_rank2,joint_rank_col,joint_rank_row, file = "Simulations/Signal_Identification_Setting5/SLIDEoutput.RData")





