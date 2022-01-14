# Load the result from SLIDE and DMMD
rm(list = ls())
# The SLIDE output is too large, which is stored outside our project folder.
load("Simulations/AdditionalSimulations/SLIDEoutput.RData")
output_slide_small = output_small
output_slide_large = output_large
output_slide_jointsmall = output_jointsmall
rm(output_small, output_large, output_jointsmall)

load("Simulations/AdditionalSimulations/AJIVE_output.RData")
output_ajive_small = output_small
output_ajive_large = output_large
output_ajive_jointsmall = output_jointsmall
rm(output_small, output_large, output_jointsmall)

load("Simulations/AdditionalSimulations/Irina_output.RData")
output_Irina_small = output_small
output_Irina_large = output_large
output_Irina_jointsmall = output_jointsmall
rm(output_small, output_large, output_jointsmall)

load("Simulations/AdditionalSimulations/output.RData")
output_small = output_small
output_large = output_large
output_jointsmall = output_jointsmall

# Load my functions
source("Simulations/MyFunction/Preliminary_Functions.R")

# Size of each data matrix.
n = 240
p = 200

# Load the original generated data to get the signal information.
load("Simulations/Data_Setting4_FixRank/Data.RData")

# Here is the true rank information.
nrep = 140
r1 = 20
r2 = 18
rc = 4
rr = 3

# Some comments about structure S that I used in SLIDE:
# The first r_c columns are (1,1) - joint
# Then the r1 - r_c columns are (1,0) - individual for 1st matrix
# Then the r2 - r_c columns are (0,1) - individual for 2nd matrix

library(ggplot2)
my_signal_error1_small = rep(NA,nrep)
my_signal_error2_small = rep(NA,nrep)
my_joint_row_error1_small = rep(NA,nrep)
my_joint_row_error2_small = rep(NA,nrep)
my_joint_col_error1_small = rep(NA,nrep)
my_joint_col_error2_small = rep(NA,nrep)

Irina_signal_error1_small = rep(NA,nrep)
Irina_signal_error2_small = rep(NA,nrep)
Irina_joint_row_error1_small = rep(NA,nrep)
Irina_joint_row_error2_small = rep(NA,nrep)
Irina_joint_col_error1_small = rep(NA,nrep)
Irina_joint_col_error2_small = rep(NA,nrep)

jive_row_error1_small = rep(NA,nrep)
jive_row_error2_small = rep(NA,nrep)
jive_row_joint_error1_small = rep(NA,nrep)
jive_row_joint_error2_small = rep(NA,nrep)
jive_col_error1_small = rep(NA,nrep)
jive_col_error2_small = rep(NA,nrep)
jive_col_joint_error1_small = rep(NA,nrep)
jive_col_joint_error2_small = rep(NA,nrep)

slide_row_error1_small = rep(NA,nrep)
slide_row_error2_small = rep(NA,nrep)
slide_row_joint_error1_small = rep(NA,nrep)
slide_row_joint_error2_small = rep(NA,nrep)
slide_col_error1_small = rep(NA,nrep)
slide_col_error2_small = rep(NA,nrep)
slide_col_joint_error1_small = rep(NA,nrep)
slide_col_joint_error2_small = rep(NA,nrep)

ajive_row_error1_small = rep(NA,nrep)
ajive_row_error2_small = rep(NA,nrep)
ajive_row_joint_error1_small = rep(NA,nrep)
ajive_row_joint_error2_small = rep(NA,nrep)
ajive_col_error1_small = rep(NA,nrep)
ajive_col_error2_small = rep(NA,nrep)
ajive_col_joint_error1_small = rep(NA,nrep)
ajive_col_joint_error2_small = rep(NA,nrep)

my_signal_error1_jointsmall = rep(NA,nrep)
my_signal_error2_jointsmall = rep(NA,nrep)
my_joint_row_error1_jointsmall = rep(NA,nrep)
my_joint_row_error2_jointsmall = rep(NA,nrep)
my_joint_col_error1_jointsmall = rep(NA,nrep)
my_joint_col_error2_jointsmall = rep(NA,nrep)

Irina_signal_error1_jointsmall = rep(NA,nrep)
Irina_signal_error2_jointsmall = rep(NA,nrep)
Irina_joint_row_error1_jointsmall = rep(NA,nrep)
Irina_joint_row_error2_jointsmall = rep(NA,nrep)
Irina_joint_col_error1_jointsmall = rep(NA,nrep)
Irina_joint_col_error2_jointsmall = rep(NA,nrep)

jive_row_error1_jointsmall = rep(NA,nrep)
jive_row_error2_jointsmall = rep(NA,nrep)
jive_row_joint_error1_jointsmall = rep(NA,nrep)
jive_row_joint_error2_jointsmall = rep(NA,nrep)
jive_col_error1_jointsmall = rep(NA,nrep)
jive_col_error2_jointsmall = rep(NA,nrep)
jive_col_joint_error1_jointsmall = rep(NA,nrep)
jive_col_joint_error2_jointsmall = rep(NA,nrep)

slide_row_error1_jointsmall = rep(NA,nrep)
slide_row_error2_jointsmall = rep(NA,nrep)
slide_row_joint_error1_jointsmall = rep(NA,nrep)
slide_row_joint_error2_jointsmall = rep(NA,nrep)
slide_col_error1_jointsmall = rep(NA,nrep)
slide_col_error2_jointsmall = rep(NA,nrep)
slide_col_joint_error1_jointsmall = rep(NA,nrep)
slide_col_joint_error2_jointsmall = rep(NA,nrep)

ajive_row_error1_jointsmall = rep(NA,nrep)
ajive_row_error2_jointsmall = rep(NA,nrep)
ajive_row_joint_error1_jointsmall = rep(NA,nrep)
ajive_row_joint_error2_jointsmall = rep(NA,nrep)
ajive_col_error1_jointsmall = rep(NA,nrep)
ajive_col_error2_jointsmall = rep(NA,nrep)
ajive_col_joint_error1_jointsmall = rep(NA,nrep)
ajive_col_joint_error2_jointsmall = rep(NA,nrep)

my_signal_error1_large = rep(NA,nrep)
my_signal_error2_large = rep(NA,nrep)
my_joint_row_error1_large = rep(NA,nrep)
my_joint_row_error2_large = rep(NA,nrep)
my_joint_col_error1_large = rep(NA,nrep)
my_joint_col_error2_large = rep(NA,nrep)

Irina_signal_error1_large = rep(NA,nrep)
Irina_signal_error2_large = rep(NA,nrep)
Irina_joint_row_error1_large = rep(NA,nrep)
Irina_joint_row_error2_large = rep(NA,nrep)
Irina_joint_col_error1_large = rep(NA,nrep)
Irina_joint_col_error2_large = rep(NA,nrep)

jive_row_error1_large = rep(NA,nrep)
jive_row_error2_large = rep(NA,nrep)
jive_row_joint_error1_large = rep(NA,nrep)
jive_row_joint_error2_large = rep(NA,nrep)
jive_col_error1_large = rep(NA,nrep)
jive_col_error2_large = rep(NA,nrep)
jive_col_joint_error1_large = rep(NA,nrep)
jive_col_joint_error2_large = rep(NA,nrep)

slide_row_error1_large = rep(NA,nrep)
slide_row_error2_large = rep(NA,nrep)
slide_row_joint_error1_large = rep(NA,nrep)
slide_row_joint_error2_large = rep(NA,nrep)
slide_col_error1_large = rep(NA,nrep)
slide_col_error2_large = rep(NA,nrep)
slide_col_joint_error1_large = rep(NA,nrep)
slide_col_joint_error2_large = rep(NA,nrep)

ajive_row_error1_large = rep(NA,nrep)
ajive_row_error2_large = rep(NA,nrep)
ajive_row_joint_error1_large = rep(NA,nrep)
ajive_row_joint_error2_large = rep(NA,nrep)
ajive_col_error1_large = rep(NA,nrep)
ajive_col_error2_large = rep(NA,nrep)
ajive_col_joint_error1_large = rep(NA,nrep)
ajive_col_joint_error2_large = rep(NA,nrep)

# Load the simulation result got from high performance server.
for (i in 1:nrep){
  my_signal_error1_small[i] = output_small[[i]]$my_signal_error1
  my_signal_error2_small[i] = output_small[[i]]$my_signal_error2
  my_joint_row_error1_small[i] = output_small[[i]]$my_joint_row_error1
  my_joint_row_error2_small[i] = output_small[[i]]$my_joint_row_error2
  my_joint_col_error1_small[i] = output_small[[i]]$my_joint_col_error1
  my_joint_col_error2_small[i] = output_small[[i]]$my_joint_col_error2

  jive_row_error1_small[i] = output_small[[i]]$jive_row_error1
  jive_row_error2_small[i] = output_small[[i]]$jive_row_error2
  jive_row_joint_error1_small[i] = output_small[[i]]$jive_row_joint_error1
  jive_row_joint_error2_small[i] = output_small[[i]]$jive_row_joint_error2
  jive_col_error1_small[i] = output_small[[i]]$jive_col_error1
  jive_col_error2_small[i] = output_small[[i]]$jive_col_error2
  jive_col_joint_error1_small[i] = output_small[[i]]$jive_col_joint_error1
  jive_col_joint_error2_small[i] = output_small[[i]]$jive_col_joint_error2
  
  my_signal_error1_jointsmall[i] = output_jointsmall[[i]]$my_signal_error1
  my_signal_error2_jointsmall[i] = output_jointsmall[[i]]$my_signal_error2
  my_joint_row_error1_jointsmall[i] = output_jointsmall[[i]]$my_joint_row_error1
  my_joint_row_error2_jointsmall[i] = output_jointsmall[[i]]$my_joint_row_error2
  my_joint_col_error1_jointsmall[i] = output_jointsmall[[i]]$my_joint_col_error1
  my_joint_col_error2_jointsmall[i] = output_jointsmall[[i]]$my_joint_col_error2
  
  jive_row_error1_jointsmall[i] = output_jointsmall[[i]]$jive_row_error1
  jive_row_error2_jointsmall[i] = output_jointsmall[[i]]$jive_row_error2
  jive_row_joint_error1_jointsmall[i] = output_jointsmall[[i]]$jive_row_joint_error1
  jive_row_joint_error2_jointsmall[i] = output_jointsmall[[i]]$jive_row_joint_error2
  jive_col_error1_jointsmall[i] = output_jointsmall[[i]]$jive_col_error1
  jive_col_error2_jointsmall[i] = output_jointsmall[[i]]$jive_col_error2
  jive_col_joint_error1_jointsmall[i] = output_jointsmall[[i]]$jive_col_joint_error1
  jive_col_joint_error2_jointsmall[i] = output_jointsmall[[i]]$jive_col_joint_error2
  
  my_signal_error1_large[i] = output_large[[i]]$my_signal_error1
  my_signal_error2_large[i] = output_large[[i]]$my_signal_error2
  my_joint_row_error1_large[i] = output_large[[i]]$my_joint_row_error1
  my_joint_row_error2_large[i] = output_large[[i]]$my_joint_row_error2
  my_joint_col_error1_large[i] = output_large[[i]]$my_joint_col_error1
  my_joint_col_error2_large[i] = output_large[[i]]$my_joint_col_error2
  
  jive_row_error1_large[i] = output_large[[i]]$jive_row_error1
  jive_row_error2_large[i] = output_large[[i]]$jive_row_error2
  jive_row_joint_error1_large[i] = output_large[[i]]$jive_row_joint_error1
  jive_row_joint_error2_large[i] = output_large[[i]]$jive_row_joint_error2
  jive_col_error1_large[i] = output_large[[i]]$jive_col_error1
  jive_col_error2_large[i] = output_large[[i]]$jive_col_error2
  jive_col_joint_error1_large[i] = output_large[[i]]$jive_col_joint_error1
  jive_col_joint_error2_large[i] = output_large[[i]]$jive_col_joint_error2
}

# Get SLIDE result. This takes some time. 3 mins.
angles <- function(Z, Zhat){
  svdZ = svd(Z)
  svdZhat = svd(Zhat)
  rZ = sum(svdZ$d > 1e-6)
  rZhat = sum(svdZhat$d > 1e-6)
  k = min(rZ, rZhat)
  cosines = svd(crossprod(svdZ$u[, 1:rZ], svdZhat$u[, 1:rZhat]))$d[1:k]
  return(list(cosines = cosines, dist = sqrt(sum(1-cosines^2))/sqrt(k)))
}

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
  
  if (!is.null(output_slide_small[[i]])){
    # True signals
    r1 = 20 - 1
    r2 = 18 - 1
    rc = 4 - 1
    rr = 3 - 1
    
    # Result from SLIDE
    Uc = output_slide_small[[i]]$slide_result_c$U
    Vc = output_slide_small[[i]]$slide_result_c$V
    temp_c = tcrossprod(Uc,Vc) # The estimated concatenation of [A1,A2] from SLIDE
    Ur = output_slide_small[[i]]$slide_result_r$U
    Vr = output_slide_small[[i]]$slide_result_r$V
    temp_r = tcrossprod(Ur,Vr) # The estimated concatenation of [t(A1),t(A2)] from SLIDE
    
    slide_col_error1_small[i] = Fnorm(temp_c[,1:p] - signal1)^2/Fnorm(signal1)^2
    slide_col_error2_small[i] = Fnorm(temp_c[,(p+1):(2*p)] - signal2)^2/Fnorm(signal2)^2
    
    slide_row_error1_small[i] = Fnorm(temp_r[,1:n] - t(signal1))^2/Fnorm(signal1)^2
    slide_row_error2_small[i] = Fnorm(temp_r[,(n+1):(2*n)] - t(signal2))^2/Fnorm(signal2)^2
    
    slide_Jc = tcrossprod(Uc[,1:rc],Vc[,1:rc])
    slide_Jr = tcrossprod(Ur[,1:rr],Vr[,1:rr])
    slide_Jc1 = slide_Jc[,1:p]
    slide_Jc2 = slide_Jc[,(p+1):(2*p)]
    slide_Jr1 = t(slide_Jr[,1:n])
    slide_Jr2 = t(slide_Jr[,(n+1):(2*n)])

    slide_col_joint_error1_small[i] = angles(temp_c[,1:p], J1_c)$dist
    slide_col_joint_error2_small[i] = angles(temp_c[,(p+1):(2*p)], J2_c)$dist
    
    slide_row_joint_error1_small[i] = angles(temp_r[,1:n], t(J1_r))$dist
    slide_row_joint_error2_small[i] = angles(temp_r[,(n+1):(2*n)], t(J2_r))$dist
  }
  
  if (!is.null(output_slide_jointsmall[[i]])){
    # True signals
    r1 = 20
    r2 = 18
    rc = 4 - 1
    rr = 3 - 1
    
    # Result from SLIDE
    Uc = output_slide_jointsmall[[i]]$slide_result_c$U
    Vc = output_slide_jointsmall[[i]]$slide_result_c$V
    temp_c = tcrossprod(Uc,Vc) # The estimated concatenation of [A1,A2] from SLIDE
    Ur = output_slide_jointsmall[[i]]$slide_result_r$U
    Vr = output_slide_jointsmall[[i]]$slide_result_r$V
    temp_r = tcrossprod(Ur,Vr) # The estimated concatenation of [t(A1),t(A2)] from SLIDE
    
    slide_col_error1_jointsmall[i] = Fnorm(temp_c[,1:p] - signal1)^2/Fnorm(signal1)^2
    slide_col_error2_jointsmall[i] = Fnorm(temp_c[,(p+1):(2*p)] - signal2)^2/Fnorm(signal2)^2
    
    slide_row_error1_jointsmall[i] = Fnorm(temp_r[,1:n] - t(signal1))^2/Fnorm(signal1)^2
    slide_row_error2_jointsmall[i] = Fnorm(temp_r[,(n+1):(2*n)] - t(signal2))^2/Fnorm(signal2)^2
    
    slide_Jc = tcrossprod(Uc[,1:rc],Vc[,1:rc])
    slide_Jr = tcrossprod(Ur[,1:rr],Vr[,1:rr])
    slide_Jc1 = slide_Jc[,1:p]
    slide_Jc2 = slide_Jc[,(p+1):(2*p)]
    slide_Jr1 = t(slide_Jr[,1:n])
    slide_Jr2 = t(slide_Jr[,(n+1):(2*n)])

    slide_col_joint_error1_jointsmall[i] = angles(temp_c[,1:p], J1_c)$dist
    slide_col_joint_error2_jointsmall[i] = angles(temp_c[,(p+1):(2*p)], J2_c)$dist
    
    slide_row_joint_error1_jointsmall[i] = angles(temp_r[,1:n], t(J1_r))$dist
    slide_row_joint_error2_jointsmall[i] = angles(temp_r[,(n+1):(2*n)], t(J2_r))$dist
  }
  
  if (!is.null(output_slide_large[[i]])){
    # True signals
    r1 = 20 + 1
    r2 = 18 + 1
    rc = 4 + 1
    rr = 3 + 1
    
    # Result from SLIDE
    Uc = output_slide_large[[i]]$slide_result_c$U
    Vc = output_slide_large[[i]]$slide_result_c$V
    temp_c = tcrossprod(Uc,Vc) # The estimated concatenation of [A1,A2] from SLIDE
    Ur = output_slide_large[[i]]$slide_result_r$U
    Vr = output_slide_large[[i]]$slide_result_r$V
    temp_r = tcrossprod(Ur,Vr) # The estimated concatenation of [t(A1),t(A2)] from SLIDE
    
    slide_col_error1_large[i] = Fnorm(temp_c[,1:p] - signal1)^2/Fnorm(signal1)^2
    slide_col_error2_large[i] = Fnorm(temp_c[,(p+1):(2*p)] - signal2)^2/Fnorm(signal2)^2
    
    slide_row_error1_large[i] = Fnorm(temp_r[,1:n] - t(signal1))^2/Fnorm(signal1)^2
    slide_row_error2_large[i] = Fnorm(temp_r[,(n+1):(2*n)] - t(signal2))^2/Fnorm(signal2)^2
    
    slide_Jc = tcrossprod(Uc[,1:rc],Vc[,1:rc])
    slide_Jr = tcrossprod(Ur[,1:rr],Vr[,1:rr])
    slide_Jc1 = slide_Jc[,1:p]
    slide_Jc2 = slide_Jc[,(p+1):(2*p)]
    slide_Jr1 = t(slide_Jr[,1:n])
    slide_Jr2 = t(slide_Jr[,(n+1):(2*n)])

    slide_col_joint_error1_large[i] = angles(temp_c[,1:p], J1_c)$dist
    slide_col_joint_error2_large[i] = angles(temp_c[,(p+1):(2*p)], J2_c)$dist
    
    slide_row_joint_error1_large[i] = angles(temp_r[,1:n], t(J1_r))$dist
    slide_row_joint_error2_large[i] = angles(temp_r[,(n+1):(2*n)], t(J2_r))$dist
  }
}

# Get AJIVE result.
for (i in 1:nrep){
  if (!is.null(output_ajive_small[[i]]$ajive_row_error1)){
    ajive_row_error1_small[i] = output_ajive_small[[i]]$ajive_row_error1
    ajive_row_error2_small[i] = output_ajive_small[[i]]$ajive_row_error2
    ajive_row_joint_error1_small[i] = output_ajive_small[[i]]$ajive_row_joint_error1
    ajive_row_joint_error2_small[i] = output_ajive_small[[i]]$ajive_row_joint_error2
    ajive_col_error1_small[i] = output_ajive_small[[i]]$ajive_col_error1
    ajive_col_error2_small[i] = output_ajive_small[[i]]$ajive_col_error2
    ajive_col_joint_error1_small[i] = output_ajive_small[[i]]$ajive_col_joint_error1
    ajive_col_joint_error2_small[i] = output_ajive_small[[i]]$ajive_col_joint_error2
  }
  if (!is.null(output_ajive_large[[i]]$ajive_row_error1)){
    ajive_row_error1_large[i] = output_ajive_large[[i]]$ajive_row_error1
    ajive_row_error2_large[i] = output_ajive_large[[i]]$ajive_row_error2
    ajive_row_joint_error1_large[i] = output_ajive_large[[i]]$ajive_row_joint_error1
    ajive_row_joint_error2_large[i] = output_ajive_large[[i]]$ajive_row_joint_error2
    ajive_col_error1_large[i] = output_ajive_large[[i]]$ajive_col_error1
    ajive_col_error2_large[i] = output_ajive_large[[i]]$ajive_col_error2
    ajive_col_joint_error1_large[i] = output_ajive_large[[i]]$ajive_col_joint_error1
    ajive_col_joint_error2_large[i] = output_ajive_large[[i]]$ajive_col_joint_error2
  }
  if (!is.null(output_ajive_jointsmall[[i]]$ajive_row_error1)){
    ajive_row_error1_jointsmall[i] = output_ajive_jointsmall[[i]]$ajive_row_error1
    ajive_row_error2_jointsmall[i] = output_ajive_jointsmall[[i]]$ajive_row_error2
    ajive_row_joint_error1_jointsmall[i] = output_ajive_jointsmall[[i]]$ajive_row_joint_error1
    ajive_row_joint_error2_jointsmall[i] = output_ajive_jointsmall[[i]]$ajive_row_joint_error2
    ajive_col_error1_jointsmall[i] = output_ajive_jointsmall[[i]]$ajive_col_error1
    ajive_col_error2_jointsmall[i] = output_ajive_jointsmall[[i]]$ajive_col_error2
    ajive_col_joint_error1_jointsmall[i] = output_ajive_jointsmall[[i]]$ajive_col_joint_error1
    ajive_col_joint_error2_jointsmall[i] = output_ajive_jointsmall[[i]]$ajive_col_joint_error2
  }
}

for (i in 1:nrep){
  Irina_signal_error1_small[i] = output_Irina_small[[i]]$Irina_signal_error1
  Irina_signal_error2_small[i] = output_Irina_small[[i]]$Irina_signal_error2
  Irina_joint_row_error1_small[i] = output_Irina_small[[i]]$Irina_joint_row_error1
  Irina_joint_row_error2_small[i] = output_Irina_small[[i]]$Irina_joint_row_error2
  Irina_joint_col_error1_small[i] = output_Irina_small[[i]]$Irina_joint_col_error1
  Irina_joint_col_error2_small[i] = output_Irina_small[[i]]$Irina_joint_col_error2

  Irina_signal_error1_large[i] = output_Irina_large[[i]]$Irina_signal_error1
  Irina_signal_error2_large[i] = output_Irina_large[[i]]$Irina_signal_error2
  Irina_joint_row_error1_large[i] = output_Irina_large[[i]]$Irina_joint_row_error1
  Irina_joint_row_error2_large[i] = output_Irina_large[[i]]$Irina_joint_row_error2
  Irina_joint_col_error1_large[i] = output_Irina_large[[i]]$Irina_joint_col_error1
  Irina_joint_col_error2_large[i] = output_Irina_large[[i]]$Irina_joint_col_error2
  
  Irina_signal_error1_jointsmall[i] = output_Irina_jointsmall[[i]]$Irina_signal_error1
  Irina_signal_error2_jointsmall[i] = output_Irina_jointsmall[[i]]$Irina_signal_error2
  Irina_joint_row_error1_jointsmall[i] = output_Irina_jointsmall[[i]]$Irina_joint_row_error1
  Irina_joint_row_error2_jointsmall[i] = output_Irina_jointsmall[[i]]$Irina_joint_row_error2
  Irina_joint_col_error1_jointsmall[i] = output_Irina_jointsmall[[i]]$Irina_joint_col_error1
  Irina_joint_col_error2_jointsmall[i] = output_Irina_jointsmall[[i]]$Irina_joint_col_error2
}

# Get some immediate sense on the comparison.
boxplot(my_signal_error1_small, jive_row_error1_small, slide_col_error1_small, slide_row_error1_small, ajive_col_error1_small, Irina_signal_error1_small)

# Plot the comparison result based on overall signal error.
y_row = c(my_signal_error1_small, my_signal_error2_small, Irina_signal_error1_small, Irina_signal_error2_small, slide_row_error1_small, 
          slide_row_error2_small, jive_row_error1_small, jive_row_error2_small, ajive_row_error1_small, ajive_row_error2_small,
          my_signal_error1_large, my_signal_error2_large, Irina_signal_error1_large, Irina_signal_error2_large, slide_row_error1_large, 
          slide_row_error2_large, jive_row_error1_large, jive_row_error2_large, ajive_row_error1_large, ajive_row_error2_large,
          my_signal_error1_jointsmall, my_signal_error2_jointsmall, Irina_signal_error1_jointsmall, Irina_signal_error2_jointsmall, slide_row_error1_jointsmall, 
          slide_row_error2_jointsmall, jive_row_error1_jointsmall, jive_row_error2_jointsmall, ajive_row_error1_jointsmall, ajive_row_error2_jointsmall)

y_col = c(my_signal_error1_small, my_signal_error2_small, Irina_signal_error1_small, Irina_signal_error2_small, slide_col_error1_small, 
          slide_col_error2_small, jive_col_error1_small, jive_col_error2_small, ajive_col_error1_small, ajive_col_error2_small,
          my_signal_error1_large, my_signal_error2_large, Irina_signal_error1_large, Irina_signal_error2_large, slide_col_error1_large, 
          slide_col_error2_large, jive_col_error1_large, jive_col_error2_large, ajive_col_error1_large, ajive_col_error2_large,
          my_signal_error1_jointsmall, my_signal_error2_jointsmall, Irina_signal_error1_jointsmall, Irina_signal_error2_jointsmall, slide_col_error1_jointsmall, 
          slide_col_error2_jointsmall, jive_col_error1_jointsmall, jive_col_error2_jointsmall, ajive_col_error1_jointsmall, ajive_col_error2_jointsmall)

method = as.factor(rep((c(rep("DMMD",2*nrep), rep("DMMD-i",2*nrep), rep("SLIDE",2*nrep), rep("JIVE",2*nrep), rep("AJIVE",2*nrep))),3))

section_setting = as.factor(c(rep(c(rep("1st signal setting (a)", nrep), rep("2nd signal setting (a)", nrep)),5), 
                              rep(c(rep("1st signal setting (b)", nrep), rep("2nd signal setting (b)", nrep)),5),
                              rep(c(rep("1st signal setting (c)", nrep), rep("2nd signal setting (c)", nrep)),5)))

gg_rowmat = data.frame(Error = y_row, Method = method, Section = section_setting)
gg_colmat = data.frame(Error = y_col, Method = method, Section = section_setting)

gg_row <- ggplot(gg_rowmat, aes(x=Method, y=Error))
gg_row <- gg_row + geom_boxplot(aes(color=Method))
gg_row <- gg_row + facet_wrap(~Section)
gg_row <- gg_row + theme_bw()
gg_row <- gg_row + ylab("Relative Error")
gg_row <- gg_row + theme(strip.background=element_rect(fill="black"))
gg_row <- gg_row + theme(strip.text=element_text(color="white", face="bold", size = 25))
gg_row <- gg_row + theme(text=element_text(size = 25))
#gg_row <- gg_row + theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())
gg_row <- gg_row + scale_colour_manual(values = c("blue", "black", "maroon", "red","green"))
gg_row = gg_row + theme(axis.text.x = element_text(angle = 45, size = 25, vjust = 0.5)) + guides(col="none")

gg_col <- ggplot(gg_colmat, aes(x=Method, y=Error))
gg_col <- gg_col + geom_boxplot(aes(color=Method))
gg_col <- gg_col + facet_wrap(~Section)
gg_col <- gg_col + theme_bw()
gg_col <- gg_col + ylab("Relative Error")
gg_col <- gg_col + theme(strip.background=element_rect(fill="black"))
gg_col <- gg_col + theme(strip.text=element_text(color="white", face="bold", size = 25))
gg_col <- gg_col + theme(text=element_text(size = 25))
#gg_col <- gg_col + theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())
gg_col <- gg_col + scale_colour_manual(values = c("blue", "black", "maroon", "red","green"))

gg_col = gg_col + theme(axis.text.x = element_text(angle = 45, size = 25, vjust = 0.5)) + guides(col="none")

pdf(file = paste(fig.path,"Signal Identification Rank Misspecified Column Decomposition Setting4.pdf",sep=""), width = 16, height = 9)
print(gg_col)
dev.off()

# Save plots
fig.path = "Simulations/AdditionalSimulations/Figures/RankMis/"
pdf(file = paste(fig.path,"Signal Identification Rank Misspecified Row Decomposition Setting4.pdf",sep=""), width = 13, height = 8)
print(gg_row)
dev.off()

pdf(file = paste(fig.path,"Signal Identification Rank Misspecified Column Decomposition Setting4.pdf",sep=""), width = 13, height = 8)
print(gg_col)
dev.off()

# Summary statistics of joint and individual estimation error (chordal distance invariant to subspace dimension)
joint_row = data.frame(my_joint_row_error1_small, my_joint_row_error2_small, Irina_joint_row_error1_small, Irina_joint_row_error2_small, slide_row_joint_error1_small, 
                       slide_row_joint_error2_small, jive_row_joint_error1_small, jive_row_joint_error2_small, ajive_row_joint_error1_small, ajive_row_joint_error2_small,
                       my_joint_row_error1_large, my_joint_row_error2_large, Irina_joint_row_error1_large, Irina_joint_row_error2_large, slide_row_joint_error1_large, 
                       slide_row_joint_error2_large, jive_row_joint_error1_large, jive_row_joint_error2_large, ajive_row_joint_error1_large, ajive_row_joint_error2_large,
                       my_joint_row_error1_jointsmall, my_joint_row_error2_jointsmall, Irina_joint_row_error1_jointsmall, Irina_joint_row_error2_jointsmall, slide_row_joint_error1_jointsmall, 
                       slide_row_joint_error2_jointsmall, jive_row_joint_error1_jointsmall, jive_row_joint_error2_jointsmall, ajive_row_joint_error1_jointsmall, ajive_row_joint_error2_jointsmall)
colMeans(joint_row)

joint_col = data.frame(my_joint_col_error1_small, my_joint_col_error2_small, Irina_joint_col_error1_small, Irina_joint_col_error2_small, slide_col_joint_error1_small, 
                       slide_col_joint_error2_small, jive_col_joint_error1_small, jive_col_joint_error2_small, ajive_col_joint_error1_small, ajive_col_joint_error2_small,
                       my_joint_col_error1_large, my_joint_col_error2_large, Irina_joint_col_error1_large, Irina_joint_col_error2_large, slide_col_joint_error1_large, 
                       slide_col_joint_error2_large, jive_col_joint_error1_large, jive_col_joint_error2_large, ajive_col_joint_error1_large, ajive_col_joint_error2_large,
                       my_joint_col_error1_jointsmall, my_joint_col_error2_jointsmall, Irina_joint_col_error1_jointsmall, Irina_joint_col_error2_jointsmall, slide_col_joint_error1_jointsmall, 
                       slide_col_joint_error2_jointsmall, jive_col_joint_error1_jointsmall, jive_col_joint_error2_jointsmall, ajive_col_joint_error1_jointsmall, ajive_col_joint_error2_jointsmall)
colMeans(joint_col)

##### Automatic creation of Table R7 with rank misspecification and joint results
# Use joint_row and joint_col created above with proper organization
############################################################################

library(tidyr)
library(tidyverse)

# Make longer
results_row = pivot_longer(joint_row, cols = 1:30, names_to = "type", values_to = "dist")
results_col = pivot_longer(joint_col, cols = 1:30, names_to = "type", values_to = "dist")
results = rbind(results_row, results_col)
# Assign method
results$method = NA
results$method[grep("my_joint", results$type)] = "DMMD"
results$method[grep("slide", results$type)] = "SLIDE"
results$method[grep("Irina", results$type)] = "DMMD-i"
results$method[grep("jive", results$type)] = "JIVE"
results$method[grep("ajive", results$type)] = "AJIVE"
# The levels of factor ensure final table has order as desired
results$method = factor(results$method, levels = c("DMMD", "DMMD-i", "SLIDE", "JIVE", "AJIVE"))
# Assign setting
results$setting = NA
results$setting[grep("_small", results$type)] = "(a)"
results$setting[grep("_large", results$type)] = "(b)"
results$setting[grep("_jointsmall", results$type)] = "(c)"
# Assign section
results$section = NA
results$section[intersect(grep("_row_", results$type), grep("_error1_", results$type))] = "1st row subspace"
results$section[intersect(grep("_row_", results$type), grep("_error2_", results$type))] = "2nd row subspace"
results$section[intersect(grep("_col_", results$type), grep("_error1_", results$type))] = "1st column subspace"
results$section[intersect(grep("_col_", results$type), grep("_error2_", results$type))] = "2nd column subspace"

# Create a table
table = results %>%
  dplyr::group_by(section, setting, method) %>%
  dplyr::summarize(meand = round(mean(dist), 2), se = round(sd(dist)/sqrt(140), 3)) %>%
  dplyr::mutate(output = paste0(meand, " (", se, ")"))

# Rearrange a table
table2 = table %>%
  dplyr::select(section, setting, method, output) %>%
  pivot_wider(names_from = "method", values_from = "output")

table2

# Make it in latex
library(xtable)
t = xtable(table2)
print(t, include.rownames = FALSE)
