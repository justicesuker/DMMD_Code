rm(list = ls())
# Load my Fnorm function
function_path = "DMMDFunctions/"
source(paste(function_path,"Preliminary_Functions.R",sep=''))
# Load the original generated data to get the signal information.
load("Simulations/SimulationData_Setting6/Data1.RData")
load("Simulations/SimulationData_Setting6/Data2.RData")

load("Simulations/Signal_Identification_Setting6_Special/DMMD_i/output.RData")
# Size of each data matrix.
n = 240
p = 200

library(ggplot2)
# Here is the rank information used to generate the data.
nrep = 140
r1 = 20
r2 = 18

total_rank1 = rep(r1,nrep)
total_rank2 = rep(r2,nrep)
min_total = rep(18,nrep)

joint_rank_col = c(rep(0, nrep/2), min_total[(nrep/2+1):nrep])
joint_rank_row = c(rep(0, nrep/4), min_total[(nrep/4+1):(nrep/2)], rep(0, nrep/4), min_total[(3*nrep/4+1):nrep])

my_signal_error1 = rep(NA,nrep)
my_signal_error2 = rep(NA,nrep)
my_joint_row_error1 = rep(NA,nrep)
my_joint_row_error2 = rep(NA,nrep)
my_ind_row_error1 = rep(NA,nrep)
my_ind_row_error2 = rep(NA,nrep)
my_joint_col_error1 = rep(NA,nrep)
my_joint_col_error2 = rep(NA,nrep)
my_ind_col_error1 = rep(NA,nrep)
my_ind_col_error2 = rep(NA,nrep)

Irina_signal_error1 = rep(NA,nrep)
Irina_signal_error2 = rep(NA,nrep)
Irina_joint_row_error1 = rep(NA,nrep)
Irina_joint_row_error2 = rep(NA,nrep)
Irina_ind_row_error1 = rep(NA,nrep)
Irina_ind_row_error2 = rep(NA,nrep)
Irina_joint_col_error1 = rep(NA,nrep)
Irina_joint_col_error2 = rep(NA,nrep)
Irina_ind_col_error1 = rep(NA,nrep)
Irina_ind_col_error2 = rep(NA,nrep)

Alliter_signal_error1 = rep(NA,nrep)
Alliter_signal_error2 = rep(NA,nrep)
Alliter_joint_row_error1 = rep(NA,nrep)
Alliter_joint_row_error2 = rep(NA,nrep)
Alliter_ind_row_error1 = rep(NA,nrep)
Alliter_ind_row_error2 = rep(NA,nrep)
Alliter_joint_col_error1 = rep(NA,nrep)
Alliter_joint_col_error2 = rep(NA,nrep)
Alliter_ind_col_error1 = rep(NA,nrep)
Alliter_ind_col_error2 = rep(NA,nrep)

# Load the simulation result got from high performance server.
for (i in 1:nrep){
  my_signal_error1[i] = output[[i]]$my_signal_error1
  my_signal_error2[i] = output[[i]]$my_signal_error2
  my_joint_row_error1[i] = output[[i]]$my_joint_row_error1
  my_joint_row_error2[i] = output[[i]]$my_joint_row_error2
  my_ind_row_error1[i] = output[[i]]$my_ind_row_error1
  my_ind_row_error2[i] = output[[i]]$my_ind_row_error2
  my_joint_col_error1[i] = output[[i]]$my_joint_col_error1
  my_joint_col_error2[i] = output[[i]]$my_joint_col_error2
  my_ind_col_error1[i] = output[[i]]$my_ind_col_error1
  my_ind_col_error2[i] = output[[i]]$my_ind_col_error2
  
  Irina_signal_error1[i] = output[[i]]$Irina_signal_error1
  Irina_signal_error2[i] = output[[i]]$Irina_signal_error2
  Irina_joint_row_error1[i] = output[[i]]$Irina_joint_row_error1
  Irina_joint_row_error2[i] = output[[i]]$Irina_joint_row_error2
  Irina_ind_row_error1[i] = output[[i]]$Irina_ind_row_error1
  Irina_ind_row_error2[i] = output[[i]]$Irina_ind_row_error2
  Irina_joint_col_error1[i] = output[[i]]$Irina_joint_col_error1
  Irina_joint_col_error2[i] = output[[i]]$Irina_joint_col_error2
  Irina_ind_col_error1[i] = output[[i]]$Irina_ind_col_error1
  Irina_ind_col_error2[i] = output[[i]]$Irina_ind_col_error2
  
  Alliter_signal_error1[i] = output[[i]]$Alliter_signal_error1
  Alliter_signal_error2[i] = output[[i]]$Alliter_signal_error2
  Alliter_joint_row_error1[i] = output[[i]]$Alliter_joint_row_error1
  Alliter_joint_row_error2[i] = output[[i]]$Alliter_joint_row_error2
  Alliter_ind_row_error1[i] = output[[i]]$Alliter_ind_row_error1
  Alliter_ind_row_error2[i] = output[[i]]$Alliter_ind_row_error2
  Alliter_joint_col_error1[i] = output[[i]]$Alliter_joint_col_error1
  Alliter_joint_col_error2[i] = output[[i]]$Alliter_joint_col_error2
  Alliter_ind_col_error1[i] = output[[i]]$Alliter_ind_col_error1
  Alliter_ind_col_error2[i] = output[[i]]$Alliter_ind_col_error2
}

rm(output)
load("Simulations/Signal_Identification_Setting6_Special/SLIDEoutput.RData")
output_slide = output
rm(output)
load("Simulations/Signal_Identification_Setting6_Special/AJIVE_output.RData")
output_ajive = output
rm(output)
load("Simulations/Signal_Identification_Setting6_Special/output.RData")

jive_row_error1 = rep(NA,nrep)
jive_row_error2 = rep(NA,nrep)
jive_row_joint_error1 = rep(NA,nrep)
jive_row_joint_error2 = rep(NA,nrep)
jive_row_ind_error1 = rep(NA,nrep)
jive_row_ind_error2 = rep(NA,nrep)
jive_col_error1 = rep(NA,nrep)
jive_col_error2 = rep(NA,nrep)
jive_col_joint_error1 = rep(NA,nrep)
jive_col_joint_error2 = rep(NA,nrep)
jive_col_ind_error1 = rep(NA,nrep)
jive_col_ind_error2 = rep(NA,nrep)

slide_row_error1 = rep(NA,nrep)
slide_row_error2 = rep(NA,nrep)
slide_row_joint_error1 = rep(NA,nrep)
slide_row_joint_error2 = rep(NA,nrep)
slide_row_ind_error1 = rep(NA,nrep)
slide_row_ind_error2 = rep(NA,nrep)
slide_col_error1 = rep(NA,nrep)
slide_col_error2 = rep(NA,nrep)
slide_col_joint_error1 = rep(NA,nrep)
slide_col_joint_error2 = rep(NA,nrep)
slide_col_ind_error1 = rep(NA,nrep)
slide_col_ind_error2 = rep(NA,nrep)

ajive_row_error1 = rep(NA,nrep)
ajive_row_error2 = rep(NA,nrep)
ajive_row_joint_error1 = rep(NA,nrep)
ajive_row_joint_error2 = rep(NA,nrep)
ajive_row_ind_error1 = rep(NA,nrep)
ajive_row_ind_error2 = rep(NA,nrep)
ajive_col_error1 = rep(NA,nrep)
ajive_col_error2 = rep(NA,nrep)
ajive_col_joint_error1 = rep(NA,nrep)
ajive_col_joint_error2 = rep(NA,nrep)
ajive_col_ind_error1 = rep(NA,nrep)
ajive_col_ind_error2 = rep(NA,nrep)

# Load the simulation result got from high performance server.
for (i in 1:nrep){
  jive_row_error1[i] = output[[i]]$jive_row_error1
  jive_row_error2[i] = output[[i]]$jive_row_error2
  jive_row_joint_error1[i] = output[[i]]$jive_row_joint_error1
  jive_row_joint_error2[i] = output[[i]]$jive_row_joint_error2
  jive_row_ind_error1[i] = output[[i]]$jive_row_ind_error1
  jive_row_ind_error2[i] = output[[i]]$jive_row_ind_error2
  jive_col_error1[i] = output[[i]]$jive_col_error1
  jive_col_error2[i] = output[[i]]$jive_col_error2
  jive_col_joint_error1[i] = output[[i]]$jive_col_joint_error1
  jive_col_joint_error2[i] = output[[i]]$jive_col_joint_error2
  jive_col_ind_error1[i] = output[[i]]$jive_col_ind_error1
  jive_col_ind_error2[i] = output[[i]]$jive_col_ind_error2
}

# Get SLIDE result.
for (i in 1:nrep){
  if (!is.null(output_slide[[i]])){
    # True signals
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
    
    # Result from SLIDE
    Uc = output_slide[[i]]$slide_result_c$U
    Vc = output_slide[[i]]$slide_result_c$V
    temp_c = tcrossprod(Uc,Vc) # The estimated concatenation of [A1,A2] from SLIDE
    Ur = output_slide[[i]]$slide_result_r$U
    Vr = output_slide[[i]]$slide_result_r$V
    temp_r = tcrossprod(Ur,Vr) # The estimated concatenation of [t(A1),t(A2)] from SLIDE
    rc = joint_rank_col[i]
    rr = joint_rank_row[i]
    
    slide_col_error1[i] = Fnorm(temp_c[,1:p] - signal1)^2/Fnorm(signal1)^2
    slide_col_error2[i] = Fnorm(temp_c[,(p+1):(2*p)] - signal2)^2/Fnorm(signal2)^2
    
    slide_row_error1[i] = Fnorm(temp_r[,1:n] - t(signal1))^2/Fnorm(signal1)^2
    slide_row_error2[i] = Fnorm(temp_r[,(n+1):(2*n)] - t(signal2))^2/Fnorm(signal2)^2
    
    if (rc != 0){
      slide_Jc = tcrossprod(Uc[,1:rc],Vc[,1:rc])
      slide_Jc1 = slide_Jc[,1:p]
      slide_Jc2 = slide_Jc[,(p+1):(2*p)]
    }
    if (rc == 0){
      slide_Jc1 = matrix(rep(0,n*p),nrow = n)
      slide_Jc2 = matrix(rep(0,n*p),nrow = n)
    }
    if (rr != 0){
      slide_Jr = tcrossprod(Ur[,1:rr],Vr[,1:rr])
      slide_Jr1 = t(slide_Jr[,1:n])
      slide_Jr2 = t(slide_Jr[,(n+1):(2*n)])
    }
    if (rr == 0){
      slide_Jr1 = matrix(rep(0,n*p),nrow = n)
      slide_Jr2 = matrix(rep(0,n*p),nrow = n)
    }
    slide_Ic1 = tcrossprod(Uc[,(rc+1):r1],Vc[,(rc+1):r1])[,1:p]
    slide_Ir1 = t(tcrossprod(Ur[,(rr+1):r1],Vr[,(rr+1):r1])[,1:n])
    if (rc == r2){
      slide_Ic2 = matrix(rep(0,n*p),nrow = n)
    }
    if (rc != r2){
      slide_Ic2 = tcrossprod(Uc[,(r1+1):(r1+r2-rc)],Vc[,(r1+1):(r1+r2-rc)])[,(p+1):(2*p)]
    }
    if (rr == r2){
      slide_Ir2 = matrix(rep(0,n*p),nrow = n)
    }
    if (rr != r2){
      slide_Ir2 = t(tcrossprod(Ur[,(r1+1):(r1+r2-rr)],Vr[,(r1+1):(r1+r2-rr)])[,(n+1):(2*n)])
    }
    slide_col_joint_error1[i] = Fnorm(slide_Jc1 - J1_c)^2
    slide_col_joint_error2[i] = Fnorm(slide_Jc2 - J2_c)^2
    slide_col_ind_error1[i] = Fnorm(slide_Ic1 - I1_c)^2
    slide_col_ind_error2[i] = Fnorm(slide_Ic2 - I2_c)^2
    
    slide_row_joint_error1[i] = Fnorm(slide_Jr1 - J1_r)^2
    slide_row_joint_error2[i] = Fnorm(slide_Jr2 - J2_r)^2
    slide_row_ind_error1[i] = Fnorm(slide_Ir1 - I1_r)^2
    slide_row_ind_error2[i] = Fnorm(slide_Ir2 - I2_r)^2
  }
}

# Get AJIVE result.
for (i in 1:nrep){
  if (!is.null(output_ajive[[i]]$ajive_row_error1)){
    ajive_row_error1[i] = output_ajive[[i]]$ajive_row_error1
    ajive_row_error2[i] = output_ajive[[i]]$ajive_row_error2
    ajive_row_joint_error1[i] = output_ajive[[i]]$ajive_row_joint_error1
    ajive_row_joint_error2[i] = output_ajive[[i]]$ajive_row_joint_error2
    ajive_row_ind_error1[i] = output_ajive[[i]]$ajive_row_ind_error1
    ajive_row_ind_error2[i] = output_ajive[[i]]$ajive_row_ind_error2
    ajive_col_error1[i] = output_ajive[[i]]$ajive_col_error1
    ajive_col_error2[i] = output_ajive[[i]]$ajive_col_error2
    ajive_col_joint_error1[i] = output_ajive[[i]]$ajive_col_joint_error1
    ajive_col_joint_error2[i] = output_ajive[[i]]$ajive_col_joint_error2
    ajive_col_ind_error1[i] = output_ajive[[i]]$ajive_col_ind_error1
    ajive_col_ind_error2[i] = output_ajive[[i]]$ajive_col_ind_error2
  }
}

# Plot the comparison result.
Signal_Error_Row = c(my_joint_row_error1, Irina_joint_row_error1, jive_row_joint_error1, slide_row_joint_error1, ajive_row_joint_error1,
                     my_ind_row_error1, Irina_ind_row_error1, jive_row_ind_error1, slide_row_ind_error1, ajive_row_ind_error1,
                     my_joint_row_error2, Irina_joint_row_error2, jive_row_joint_error2, slide_row_joint_error2, ajive_row_joint_error2,
                     my_ind_row_error2, Irina_ind_row_error2, jive_row_ind_error2, slide_row_ind_error2, ajive_row_ind_error2)

Method = rep(c(rep('DMMD',nrep),rep('DMMD-i',nrep),rep('JIVE',nrep),rep('SLIDE',nrep),rep('AJIVE',nrep)),4)

Section = factor(c(rep('1st Joint',5*nrep),rep('1st Ind',5*nrep),
                   rep('2nd Joint',5*nrep),rep('2nd Ind',5*nrep)))

Note = rep(c(rep('1-35',35),rep('36-70',35),rep('71-105',35),rep('106-140',35)),20)

ggplot_row_mat = data.frame(Signal_Error = Signal_Error_Row, Method, Section, Note)

levels(ggplot_row_mat$Method)<- c('AJIVE', 'DMMD-i', 'DMMD', 'JIVE', 'SLIDE')

levels(ggplot_row_mat$Section)[levels(ggplot_row_mat$Section)=="1st Joint"] <- "'1st Joint'"
levels(ggplot_row_mat$Section)[levels(ggplot_row_mat$Section)=="1st Ind"] <- "'1st Ind'"
levels(ggplot_row_mat$Section)[levels(ggplot_row_mat$Section)=="2nd Joint"] <- "'2nd Joint'"
levels(ggplot_row_mat$Section)[levels(ggplot_row_mat$Section)=="2nd Ind"] <- "'2nd Ind'"

ggplot_row_mat$Note <- factor(ggplot_row_mat$Note,
                              levels=c("1-35","36-70","71-105", "106-140"),
                              labels=c(expression(paste(r[c] == 0, ",", r[r] == 0)),
                                       expression(paste(r[c] == 0, ",", r[r] == r[2])),
                                       expression(paste(r[c] == r[2], ",", r[r] == 0)),
                                       expression(paste(r[c], " = ", r[r] == r[2]))))

gg_row <- ggplot(ggplot_row_mat, aes(x=Method, y=Signal_Error)) + geom_boxplot(aes(fill=Method)) + facet_grid(rows = vars(Note), cols = vars(Section), labeller = label_parsed, scales = "free") + theme_bw() + ylab("Absolute Estimation Error")+ theme(strip.background=element_rect(fill="black")) + theme(strip.text=element_text(color="white", face="bold",size = 25), text=element_text(size = 25)) + theme(axis.title.x=element_blank(), axis.ticks.x=element_blank(), axis.text.x=element_blank(),  legend.position = "bottom")
gg_row <- gg_row + scale_fill_manual(values = c("#ffffcc", "#a1dab4", "#41b6c4", "#2c7fb8","#253494"))
gg_row

Signal_Error_Col = c(my_joint_col_error1, Irina_joint_col_error1, jive_col_joint_error1, slide_col_joint_error1, ajive_col_joint_error1,
                     my_ind_col_error1, Irina_ind_col_error1, jive_col_ind_error1, slide_col_ind_error1, ajive_col_ind_error1,
                     my_joint_col_error2, Irina_joint_col_error2, jive_col_joint_error2, slide_col_joint_error2, ajive_col_joint_error2,
                     my_ind_col_error2, Irina_ind_col_error2, jive_col_ind_error2, slide_col_ind_error2, ajive_col_ind_error2)

ggplot_col_mat = data.frame(Signal_Error = Signal_Error_Col, Method, Section, Note)

levels(ggplot_col_mat$Method)<- c('AJIVE', 'Irina', 'DMMD', 'JIVE', 'SLIDE')

levels(ggplot_col_mat$Section)[levels(ggplot_col_mat$Section)=="1st Joint"] <- "'1st Joint'"
levels(ggplot_col_mat$Section)[levels(ggplot_col_mat$Section)=="1st Ind"] <- "'1st Ind'"
levels(ggplot_col_mat$Section)[levels(ggplot_col_mat$Section)=="2nd Joint"] <- "'2nd Joint'"
levels(ggplot_col_mat$Section)[levels(ggplot_col_mat$Section)=="2nd Ind"] <- "'2nd Ind'"

ggplot_col_mat$Note <- factor(ggplot_col_mat$Note,
                              levels=c("1-35","36-70","71-105", "106-140"),
                              labels=c(expression(paste(r[c] == 0, ",", r[r] == 0)),
                                       expression(paste(r[c] == 0, ",", r[r] == r[2])),
                                       expression(paste(r[c] == r[2], ",", r[r] == 0)),
                                       expression(paste(r[c], " = ", r[r] == r[2]))))

gg_col <- ggplot(ggplot_col_mat, aes(x=Method, y=Signal_Error)) + geom_boxplot(aes(fill=Method)) + facet_grid(rows = vars(Note), cols = vars(Section), labeller = label_parsed, scales = "free") + theme_bw() + ylab("Absolute Estimation Error")+ theme(strip.background=element_rect(fill="black")) + theme(strip.text=element_text(color="white", face="bold",size = 25), text=element_text(size = 25)) + theme(axis.title.x=element_blank(), axis.ticks.x=element_blank(), axis.text.x=element_blank(),  legend.position = "bottom")
gg_col <- gg_col + scale_fill_manual(values = c("#ffffcc", "#a1dab4", "#41b6c4", "#2c7fb8","#253494"))
gg_col

# Save
fig.path = "Simulations/Signal_Identification_Setting6_Special/FinalFigures/Draft/"
pdf(file = paste(fig.path,"Signal Identification for Row Decomposition_All_Setting6.pdf",sep=""), width = 9, height = 9)
print(gg_row)
dev.off()

pdf(file = paste(fig.path,"Signal Identification for Column Decomposition_All_Setting6.pdf",sep=""), width = 9, height = 9)
print(gg_col)
dev.off()