rm(list = ls())
load("Simulations/RankEstimation_Setting1/SLIDEoutput.RData")
output_slide = output
rm(output)
load("Simulations/RankEstimation_Setting1/AJIVEoutput.RData")
output_ajive = output
rm(output)
load("Simulations/RankEstimation_Setting1/output.RData")

# Load my Fnorm function
function_path = "DMMDFunctions/"
source(paste(function_path,"Preliminary_Functions.R",sep=''))

# n = 240
# p = 200
# nrep = 140

# Define a function that could extract the joint rank of a structure S output from SLIDE
# S is the output from SLIDE function
getRank <- function(S){
  p = dim(S)[2]
  r1 = 0
  r2 = 0
  rjoint = 0
  for (i in 1:p){
    if (S[1,i] == 1 & S[2,i] == 1){
      r1 = r1 + 1
      r2 = r2 + 1
      rjoint = rjoint + 1
    }
    if (S[1,i] == 1 & S[2,i] == 0){
      r1 = r1 + 1
    }
    if (S[1,i] == 0 & S[2,i] == 1){
      r2 = r2 + 1
    }
  }
  return(list(r1 = r1, r2 = r2, rjoint = rjoint))
}

library(ggplot2)
n = length(output)

my_error_joint_rank_col_PL = rep(NA,n)
my_error_joint_rank_row_PL = rep(NA,n)
my_error_total_rank1_PL = rep(NA,n)
my_error_total_rank2_PL = rep(NA,n)

jive_error_joint_rank_c = rep(NA,n)
jive_error_total_rank1_c = rep(NA,n)
jive_error_total_rank2_c = rep(NA,n)
jive_error_joint_rank_r = rep(NA,n)
jive_error_total_rank1_r = rep(NA,n)
jive_error_total_rank2_r = rep(NA,n)

slide_error_joint_rank_c = rep(NA,n)
slide_error_total_rank1_c = rep(NA,n)
slide_error_total_rank2_c = rep(NA,n)
slide_error_joint_rank_r = rep(NA,n)
slide_error_total_rank1_r = rep(NA,n)
slide_error_total_rank2_r = rep(NA,n)

ajive_error_joint_rank_c = rep(NA,n)
ajive_error_joint_rank_r = rep(NA,n)

my_error_total_rank1_ED = rep(NA,n)
my_error_total_rank2_ED = rep(NA,n)

for (i in 1:n){
  my_error_joint_rank_col_PL[i] = output[[i]]$my_error_joint_rank_col_PL
  my_error_joint_rank_row_PL[i] = output[[i]]$my_error_joint_rank_row_PL
  my_error_total_rank1_PL[i] = output[[i]]$my_error_total_rank1_PL
  my_error_total_rank2_PL[i] = output[[i]]$my_error_total_rank2_PL
  my_error_total_rank1_ED[i] = output[[i]]$my_error_total_rank1_ED
  my_error_total_rank2_ED[i] = output[[i]]$my_error_total_rank2_ED
  
  jive_error_joint_rank_c[i] = output[[i]]$jive_error_joint_rank_c
  jive_error_total_rank1_c[i] = output[[i]]$jive_error_total_rank1_c
  jive_error_total_rank2_c[i] = output[[i]]$jive_error_total_rank2_c
  jive_error_joint_rank_r[i] = output[[i]]$jive_error_joint_rank_r
  jive_error_total_rank1_r[i] = output[[i]]$jive_error_total_rank1_r
  jive_error_total_rank2_r[i] = output[[i]]$jive_error_total_rank2_r
  
  if (!is.null(output_slide[[i]]$slide_result_c)){
    slide_error_joint_rank_c[i] = getRank(output_slide[[i]]$slide_result_c)$rjoint - joint_rank_col[i]
    slide_error_total_rank1_c[i] = getRank(output_slide[[i]]$slide_result_c)$r1 - total_rank1[i]
    slide_error_total_rank2_c[i] = getRank(output_slide[[i]]$slide_result_c)$r2 - total_rank2[i]
  }
  if (!is.null(output_slide[[i]]$slide_result_r)){
    slide_error_joint_rank_r[i] = getRank(output_slide[[i]]$slide_result_r)$rjoint - joint_rank_row[i]
    slide_error_total_rank1_r[i] = getRank(output_slide[[i]]$slide_result_r)$r1 - total_rank1[i]
    slide_error_total_rank2_r[i] = getRank(output_slide[[i]]$slide_result_r)$r2 - total_rank2[i]
  }
  
  ajive_error_joint_rank_c[i] = output_ajive[[i]]$ajive_error_joint_rank_c
  ajive_error_joint_rank_r[i] = output_ajive[[i]]$ajive_error_joint_rank_r
}

y_joint = c(my_error_joint_rank_col_PL, jive_error_joint_rank_c, 
            ajive_error_joint_rank_c, slide_error_joint_rank_c,
            my_error_joint_rank_row_PL, jive_error_joint_rank_r,
            ajive_error_joint_rank_r, slide_error_joint_rank_r)

method_joint = as.factor(rep(c(rep("PL",n),rep("JIVE",n),
                               rep("AJIVE",n), rep("SLIDE",n)),2))

RowCol_joint = c(rep("Joint Column Rank",4*n),rep("Joint Row Rank",4*n))
joint_mat = data.frame(Error = y_joint,
                       Method = method_joint,
                       RowCol = RowCol_joint)

gg_joint <- ggplot(joint_mat, aes(x=Method, y=Error)) + 
  geom_boxplot(aes(fill = Method)) + facet_wrap(~RowCol) + theme_bw() +
  theme(strip.background=element_rect(fill="black"), text=element_text(size = 25))
gg_joint <- gg_joint + theme(strip.text=element_text(color="white", face="bold", size = 25))
gg_joint <- gg_joint + scale_fill_manual(values = c("#ffffcc", "#a1dab4", "#41b6c4", "#225ea8"))

gg_joint <- gg_joint + theme(axis.text.x = element_text(angle = 45, size = 20, vjust = 1, hjust = 1)) + guides(fill="none") + ylab("Estimated Rank - True Rank")
gg_joint

y_total = c(my_error_total_rank1_PL, my_error_total_rank1_ED, 
            jive_error_total_rank1_c, jive_error_total_rank1_r,
            slide_error_total_rank1_c,slide_error_total_rank1_r,
            my_error_total_rank2_PL, my_error_total_rank2_ED,
            jive_error_total_rank2_c, jive_error_total_rank2_r,
            slide_error_total_rank2_c,slide_error_total_rank1_r)

method_total = rep(c(rep("PL",n),
                     rep("ED",n),
                     rep("JIVE(R)",n),
                     rep("JIVE(C)",n),
                     rep("SLIDE(R)",n),
                     rep("SLIDE(C)",n)),2)
MatrixName = c(rep("First Matrix",6*n),rep("Second Matrix",6*n))

total_mat = data.frame(Error = y_total,
                       Method = method_total,
                       Matrix = MatrixName)

gg_total <- ggplot(total_mat, aes(x=Method, y=Error)) + geom_boxplot(aes(fill=Method)) + facet_wrap(~Matrix) + theme_bw()
gg_total <- gg_total + theme(strip.background=element_rect(fill="black"), text=element_text(size = 25)) + theme(strip.text=element_text(color="white", face="bold", size = 25))
gg_total <- gg_total + scale_fill_manual(values = c("#ffffcc","#c7e9b4","#7fcdbb","#41b6c4","#2c7fb8","#253494"))
gg_total <- gg_total + theme(axis.text.x = element_text(angle = 45, size = 20, vjust = 1, hjust = 1)) + guides(fill="none") + ylab("Estimated Rank - True Rank")
gg_total

# Draw boxplot
fig.path = "Simulations/RankEstimation_Setting1/Figures/Draft/"
pdf(file = paste(fig.path,"Comparison on Joint Rank Estimation_Setting1.pdf",sep=""), width = 10, height = 7)
print(gg_joint)
dev.off()

pdf(file = paste(fig.path,"Comparison on Total Rank Estimation_Setting1.pdf",sep=""), width = 10, height = 7)
print(gg_total)
dev.off()
