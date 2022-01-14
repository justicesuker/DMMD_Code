rm(list = ls())
win_data = read.csv("Application/Soccer/Winning Team Soccer.csv")
lose_data = read.csv("Application/Soccer/Losing Team Soccer.csv")
function_path1 = "DMMDFunctions/"
function_path2 = "OtherFunctions/"

source(paste(function_path1,"Angle_Calculation.R",sep=''))
source(paste(function_path1,"Profile_Likelihood_Rank_Selection.R",sep=''))
source(paste(function_path1,"DoubleMatchedMatrixDecomposition.R",sep=''))
source(paste(function_path1,"DMMD_iterative.R",sep=''))
source(paste(function_path1,"FindOptMatrix.R",sep=''))
source(paste(function_path1,"Preliminary_Functions.R",sep=''))
source(paste(function_path2,"Select_ED_Rank.R",sep=''))

dim(win_data)
win_data[1:5,1:5]
colnames(win_data)
names_c = c("Full.Time.Goals","Half.Time.Goals","Shots","Shots.on.Target",
  "Hit.Woodwork", "Corners", "Fouls.Committed",
  "Offsides","Yellow.Cards", "Red.Cards")

p = dim(lose_data)[2]
w_data = win_data[,2:p]
l_data = lose_data[,2:p]
colnames(w_data) = names_c
colnames(l_data) = names_c
head(w_data)
result_equal = DMMD_v2(w_data,l_data)
result_equal$Rank

result_ed = DMMD_v2(w_data,l_data,method = "ED")
result_ed$Rank

result_unequal = DMMD_v2(w_data,l_data, variance1 = "unequal")
result_unequal$Rank

plot(svd(l_data)$d, main = "Scree plot for losing team")
plot(svd(w_data)$d, main = "Scree plot for winning team")
svd_w = svd(w_data)
svd_l = svd(l_data)

ProfileLikCluster(svd(w_data)$d,"unequal")
ProfileLikCluster(svd(w_data)$d,"equal")
ProfileLikCluster(svd(l_data)$d,"unequal")
ProfileLikCluster(svd(l_data)$d,"equal")

I1 = result_unequal$`Row Decomposition`$`Individual Row 1`
J1 = result_unequal$`Row Decomposition`$`Joint Row 1`
I2 = result_unequal$`Row Decomposition`$`Individual Row 2`
J2 = result_unequal$`Row Decomposition`$`Joint Row 2`
E1 = result_unequal$Error$Error1
E2 = result_unequal$Error$Error2

print(c(Fnorm(J1),Fnorm(I1),Fnorm(E1),Fnorm(w_data)))
print(c(Fnorm(J2),Fnorm(I2),Fnorm(E2),Fnorm(l_data)))

I1[1,]/I1[1,1]
J1[1,]/J1[1,1]
E2[1,]

# Create a function that adds % to each number in a vector
# Input: vec, the numeric singular value vector.
# Output: a string vector that has the percentage information.
AddPercentage <- function(vec,digit = 1){
  new_vec = round(vec^2/sum(vec^2)*100,digit)
  new_vec = format(new_vec,nsmall = digit)
  new_vec = as.character(new_vec)
  for (i in 1:length(vec)){
    new_vec[i] = paste(new_vec[i],"%",sep = "")
  }
  return(new_vec)
}

svd_result = svd(w_data)$d
variation_vec = AddPercentage(svd_result)
print(variation_vec)

svd_result2 = svd(l_data)$d
plot(svd_result2)
variation_vec2 = AddPercentage(svd_result2)
print(variation_vec2)

#If we want to have more than 90% variance explained, we need r1 = 2, 
# which is in agreement of PL using unequal variance assumption. 
I1_u = result_unequal$`Row Decomposition`$`Individual Row 1`
J1_u = result_unequal$`Row Decomposition`$`Joint Row 1`
print(J1_u[1,]/J1_u[1,1])
print(I1_u[1,]/I1_u[1,1])

# Now let's try to use Iterative DMMD:
r1 = 2
r2 = 1
rc = 1
rr = 1
DMMDi_result = DMMD_i(as.matrix(w_data),as.matrix(l_data), r1, r2, rc, rr)

J1_i = DMMDi_result$A1 %*% projection(DMMDi_result$N,ortho = TRUE)
I1_i = DMMDi_result$A1 - J1_i

J2_i = DMMDi_result$A2 %*% projection(DMMDi_result$N,ortho = TRUE)
I2_i = DMMDi_result$A2 - J2_i

I1_i[1,]/I1_i[1,1]
J1_i[1,]/J1_i[1,1]

rbind("DMMD win signal" = I1[1,]/I1[1,1], "DMMD-i win signal" = I1_i[1,]/I1_i[1,1])
rbind("DMMD joint signal" = J1[1,]/J1[1,1], "DMMD-i joint signal" = J1_i[1,]/J1_i[1,1])
