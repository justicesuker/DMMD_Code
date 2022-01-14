# Let's run a smaller example.
rm(list=ls())
function_path1 = "DMMDFunctions/"
function_path2 = "OtherFunctions/"

source(paste(function_path1,"Angle_Calculation.R",sep=''))
source(paste(function_path1,"Profile_Likelihood_Rank_Selection.R",sep=''))
source(paste(function_path1,"DoubleMatchedMatrixDecomposition.R",sep=''))
source(paste(function_path1,"DMMD_iterative.R",sep=''))
source(paste(function_path1,"FindOptMatrix.R",sep=''))
source(paste(function_path1,"Preliminary_Functions.R",sep=''))
source(paste(function_path2,"DoubleMatchedDataGen.R",sep=''))
source(paste(function_path2,"slide_prelim.R",sep=''))

library(microbenchmark)
library(r.jive)
library(ajive)
library(SLIDE)

set.seed(37)

n = 100
# Setting when p = 800
p_l = 800

total_rank1_l = 25
total_rank2_l = 20
r_c_l = 10
r_r_l = 5

std_est1_l = sqrt(total_rank1_l/(n*p_l))
std_est2_l = sqrt(total_rank2_l/(n*p_l))

# Generate data
data_l = DoubleDataGen3(n = n, p = p_l, rank = c(total_rank1_l, total_rank2_l), 
                      joint_rank_col = r_c_l, joint_rank_row = r_r_l, nrep = 1, 
                      std1 = std_est1_l, std2 = std_est2_l)
X1_l = data_l$X1_list[[1]]
X2_l = data_l$X2_list[[1]]

# Set up the parameters
################# JIVE ##################
# Matched columns - row decomposition
data_row_l = list(X1_l,X2_l)

# Matched rows - column decomposition
data_col_l = list(t(X1_l),t(X2_l))

initial_signal_ranks_l <- c(total_rank1_l, total_rank2_l)
temprankA_row_l = initial_signal_ranks_l - r_r_l
temprankA_col_l = initial_signal_ranks_l - r_c_l

################# SLIDE ###################
r1_l = total_rank1_l
r2_l = total_rank2_l
pvec_l = c(p_l,p_l)
nvec = c(n,n)
# Column concatenation - matched rows - column decomposition
data_temp_c_l = cbind(X1_l,X2_l)
# row concatenation - matched columns - row decomposition
data_temp_r_l = cbind(t(X1_l),t(X2_l))
S_c_l = getS(r1_l, r2_l, r_c_l, ndata = 2)
S_r_l = getS(r1_l, r2_l, r_r_l, ndata = 2)

# Setting when p = 80
p_s = 80
total_rank1_s = 10
total_rank2_s = 8
r_c_s = 4
r_r_s = 3

std_est1_s = sqrt(total_rank1_s/(n*p_s))
std_est2_s = sqrt(total_rank2_s/(n*p_s))

# Generate data
data_s = DoubleDataGen3(n = n, p = p_s, rank = c(total_rank1_s, total_rank2_s), 
                        joint_rank_col = r_c_s, joint_rank_row = r_r_s, nrep = 1, 
                        std1 = std_est1_s, std2 = std_est2_s)
X1_s = data_s$X1_list[[1]]
X2_s = data_s$X2_list[[1]]

# Set up the parameters
################# JIVE ##################
# Matched columns - row decomposition
data_row_s = list(X1_s,X2_s)

# Matched rows - column decomposition
data_col_s = list(t(X1_s),t(X2_s))

initial_signal_ranks_s <- c(total_rank1_s, total_rank2_s)
temprankA_row_s = initial_signal_ranks_s - r_r_s
temprankA_col_s = initial_signal_ranks_s - r_c_s

################# SLIDE ###################
r1_s = total_rank1_s
r2_s = total_rank2_s
pvec_s = c(p_s,p_s)
nvec = c(n,n)
# Column concatenation - matched rows - column decomposition
data_temp_c_s = cbind(X1_s,X2_s)
# row concatenation - matched columns - row decomposition
data_temp_r_s = cbind(t(X1_s),t(X2_s))
S_c_s = getS(r1_s, r2_s, r_c_s, ndata = 2)
S_r_s = getS(r1_s, r2_s, r_r_s, ndata = 2)

# Compare DMMD with JIVE, SLIDE, AJIVE on given ranks
# Will take ~1min to run.
microbenchmark(DMMD_l_given = DMMD_v2(X1_l, X2_l, r1 = r1_l, r2 = r2_l, joint_rank_c = r_c_l, joint_rank_r = r_r_l),
               Irina_result_l_given = DMMD_Irina(X1_l, X2_l, r1 = r1_l, r2 = r2_l, rc = r_c_l, rr = r_r_l),
               jive_result_r_l_given = jive(data_row_l, rankJ = r_r_l, rankA = temprankA_row_l, method = 'given', center = FALSE, scale = FALSE),
               jive_result_c_l_given = jive(data_col_l, rankJ = r_c_l, rankA = temprankA_col_l, method = 'given', center = FALSE, scale = FALSE),
               slide_result_c_l_given = slide_givenS(data_temp_c_l, pvec_l, S_c_l, standardized = T),
               slide_result_r_l_given = slide_givenS(data_temp_r_l, nvec, S_r_l, standardized = T),
               ajive_result_r_l_given = ajive(data_col_l, initial_signal_ranks_l, joint_rank = r_r_l),
               ajive_result_c_l_given = ajive(data_row_l, initial_signal_ranks_l, joint_rank = r_c_l),
               
               DMMD_l = DMMD_v2(X1_l, X2_l),
               Irina_result_l = DMMD_Irina_unknown(X1_l,X2_l),
               jive_result_r_l = jive(data_row_l, center = FALSE, scale = FALSE),
               jive_result_c_l = jive(data_col_l, center = FALSE, scale = FALSE),
               slide_result_c_l = slide(data_temp_c_l, pvec_l),
               slide_result_r_l = slide(data_temp_r_l, nvec),
               ajive_result_r_l = ajive(data_col_l, initial_signal_ranks_l),
               ajive_result_c_l = ajive(data_row_l, initial_signal_ranks_l),
               
               # p = 80
               DMMD_s_given = DMMD_v2(X1_s, X2_s, r1 = r1_s, r2 = r2_s, joint_rank_c = r_c_s, joint_rank_r = r_r_s),
               Irina_result_s_given = DMMD_Irina(X1_s, X2_s, r1 = r1_s, r2 = r2_s, rc = r_c_s, rr = r_r_s),
               jive_result_r_s_given = jive(data_row_s, rankJ = r_r_s, rankA = temprankA_row_s, method = 'given', center = FALSE, scale = FALSE),
               jive_result_c_s_given = jive(data_col_s, rankJ = r_c_s, rankA = temprankA_col_s, method = 'given', center = FALSE, scale = FALSE),
               slide_result_c_s_given = slide_givenS(data_temp_c_s, pvec_s, S_c_s, standardized = T),
               slide_result_r_s_given = slide_givenS(data_temp_r_s, nvec, S_r_s, standardized = T),
               ajive_result_r_s_given = ajive(data_col_s, initial_signal_ranks_s, joint_rank = r_r_s),
               ajive_result_c_s_given = ajive(data_row_s, initial_signal_ranks_s, joint_rank = r_c_s),
               
               DMMD_s = DMMD_v2(X1_s, X2_s),
               Irina_result_s = DMMD_Irina_unknown(X1_s,X2_s),
               jive_result_r_s = jive(data_row_s, center = FALSE, scale = FALSE),
               jive_result_c_s = jive(data_col_s, center = FALSE, scale = FALSE),
               slide_result_c_s = slide(data_temp_c_s, pvec_s),
               slide_result_r_s = slide(data_temp_r_s, nvec),
               ajive_result_r_s = ajive(data_col_s, initial_signal_ranks_s),
               ajive_result_c_s = ajive(data_row_s, initial_signal_ranks_s),
               
               times = 1
)

# Unit: milliseconds
# expr                            min           lq         mean       median           uq          max neval
# DMMD_l_given             11818.3256   11818.3256   11818.3256   11818.3256   11818.3256   11818.3256     1
# Irina_result_l_given     48394.5877   48394.5877   48394.5877   48394.5877   48394.5877   48394.5877     1
# jive_result_r_l_given    21076.8957   21076.8957   21076.8957   21076.8957   21076.8957   21076.8957     1
# jive_result_c_l_given     7838.7361    7838.7361    7838.7361    7838.7361    7838.7361    7838.7361     1
# slide_result_c_l_given    1101.5663    1101.5663    1101.5663    1101.5663    1101.5663    1101.5663     1
# slide_result_r_l_given    1036.4655    1036.4655    1036.4655    1036.4655    1036.4655    1036.4655     1
# ajive_result_r_l_given     443.8832     443.8832     443.8832     443.8832     443.8832     443.8832     1
# ajive_result_c_l_given     250.1484     250.1484     250.1484     250.1484     250.1484     250.1484     1

# DMMD_l                   24104.7921   24104.7921   24104.7921   24104.7921   24104.7921   24104.7921     1
# Irina_result_l          618263.4256  618263.4256  618263.4256  618263.4256  618263.4256  618263.4256     1
# jive_result_r_l         219729.0462  219729.0462  219729.0462  219729.0462  219729.0462  219729.0462     1
# jive_result_c_l          12630.0666   12630.0666   12630.0666   12630.0666   12630.0666   12630.0666     1
# slide_result_c_l        317270.6047  317270.6047  317270.6047  317270.6047  317270.6047  317270.6047     1
# slide_result_r_l       1226774.9264 1226774.9264 1226774.9264 1226774.9264 1226774.9264 1226774.9264     1
# ajive_result_r_l         66645.4573   66645.4573   66645.4573   66645.4573   66645.4573   66645.4573     1
# ajive_result_c_l         42452.0821   42452.0821   42452.0821   42452.0821   42452.0821   42452.0821     1

# DMMD_s_given               362.9350     362.9350     362.9350     362.9350     362.9350     362.9350     1
# Irina_result_s_given      1412.6003    1412.6003    1412.6003    1412.6003    1412.6003    1412.6003     1
# jive_result_r_s_given     8253.0106    8253.0106    8253.0106    8253.0106    8253.0106    8253.0106     1
# jive_result_c_s_given     6115.1304    6115.1304    6115.1304    6115.1304    6115.1304    6115.1304     1
# slide_result_c_s_given      74.3331      74.3331      74.3331      74.3331      74.3331      74.3331     1
# slide_result_r_s_given      53.4265      53.4265      53.4265      53.4265      53.4265      53.4265     1
# ajive_result_r_s_given      32.1215      32.1215      32.1215      32.1215      32.1215      32.1215     1
# ajive_result_c_s_given      22.8535      22.8535      22.8535      22.8535      22.8535      22.8535     1

# DMMD_s                     794.8818     794.8818     794.8818     794.8818     794.8818     794.8818     1
# Irina_result_s           10918.1485   10918.1485   10918.1485   10918.1485   10918.1485   10918.1485     1
# jive_result_r_s           2065.8257    2065.8257    2065.8257    2065.8257    2065.8257    2065.8257     1
# jive_result_c_s           2227.6161    2227.6161    2227.6161    2227.6161    2227.6161    2227.6161     1
# slide_result_c_s         63740.2384   63740.2384   63740.2384   63740.2384   63740.2384   63740.2384     1
# slide_result_r_s         36202.2834   36202.2834   36202.2834   36202.2834   36202.2834   36202.2834     1
# ajive_result_r_s          6321.5968    6321.5968    6321.5968    6321.5968    6321.5968    6321.5968     1
# ajive_result_c_s          6539.8837    6539.8837    6539.8837    6539.8837    6539.8837    6539.8837     1