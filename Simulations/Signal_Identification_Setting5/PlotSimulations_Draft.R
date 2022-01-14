rm(list = ls())
# Load my Fnorm function
function_path = "DMMDFunctions/"
source(paste(function_path,"Preliminary_Functions.R",sep=''))
# Load the original generated data to get the signal information.
load("Simulations/Data_Setting5_FixRank_SNR0.5/Data.RData")

load("Simulations/Signal_Identification_Setting5/DMMD_i/output.RData")

# Size of each data matrix.
n = 240
p = 200

library(ggplot2)
# Here is the rank information used to generate the data.
nrep = 140
r1 = 20
r2 = 18
rc = 4
rr = 3
total_rank1 = rep(r1,nrep)
total_rank2 = rep(r2,nrep)
joint_rank_col = rep(rc,nrep)
joint_rank_row = rep(rr,nrep)

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
load("Simulations/Signal_Identification_Setting5/SLIDEoutput.RData")
output_slide = output
rm(output)
load("Simulations/Signal_Identification_Setting5/AJIVE_output.RData")
output_ajive = output
rm(output)
load("Simulations/Signal_Identification_Setting5/output.RData")

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
    
    slide_col_error1[i] = Fnorm(temp_c[,1:p] - signal1)^2/Fnorm(signal1)^2
    slide_col_error2[i] = Fnorm(temp_c[,(p+1):(2*p)] - signal2)^2/Fnorm(signal2)^2
    
    slide_row_error1[i] = Fnorm(temp_r[,1:n] - t(signal1))^2/Fnorm(signal1)^2
    slide_row_error2[i] = Fnorm(temp_r[,(n+1):(2*n)] - t(signal2))^2/Fnorm(signal2)^2
    
    slide_Jc = tcrossprod(Uc[,1:rc],Vc[,1:rc])
    slide_Jr = tcrossprod(Ur[,1:rr],Vr[,1:rr])
    slide_Jc1 = slide_Jc[,1:p]
    slide_Jc2 = slide_Jc[,(p+1):(2*p)]
    slide_Jr1 = t(slide_Jr[,1:n])
    slide_Jr2 = t(slide_Jr[,(n+1):(2*n)])
    
    slide_Ic1 = tcrossprod(Uc[,(rc+1):r1],Vc[,(rc+1):r1])[,1:p]
    slide_Ic2 = tcrossprod(Uc[,(r1+1):(r1+r2-rc)],Vc[,(r1+1):(r1+r2-rc)])[,(p+1):(2*p)]
    slide_Ir1 = t(tcrossprod(Ur[,(rr+1):r1],Vr[,(rr+1):r1])[,1:n])
    slide_Ir2 = t(tcrossprod(Ur[,(r1+1):(r1+r2-rr)],Vr[,(r1+1):(r1+r2-rr)])[,(n+1):(2*n)])
    
    slide_col_joint_error1[i] = Fnorm(slide_Jc1 - J1_c)^2/Fnorm(J1_c)^2
    slide_col_joint_error2[i] = Fnorm(slide_Jc2 - J2_c)^2/Fnorm(J2_c)^2
    slide_col_ind_error1[i] = Fnorm(slide_Ic1 - I1_c)^2/Fnorm(I1_c)^2
    slide_col_ind_error2[i] = Fnorm(slide_Ic2 - I2_c)^2/Fnorm(I2_c)^2
    
    slide_row_joint_error1[i] = Fnorm(slide_Jr1 - J1_r)^2/Fnorm(J1_r)^2
    slide_row_joint_error2[i] = Fnorm(slide_Jr2 - J2_r)^2/Fnorm(J2_r)^2
    slide_row_ind_error1[i] = Fnorm(slide_Ir1 - I1_r)^2/Fnorm(I1_r)^2
    slide_row_ind_error2[i] = Fnorm(slide_Ir2 - I2_r)^2/Fnorm(I2_r)^2
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
y_row = c(my_signal_error1, my_signal_error2, my_joint_row_error1, my_joint_row_error2, my_ind_row_error1, my_ind_row_error2,
          Irina_signal_error1, Irina_signal_error2, Irina_joint_row_error1, Irina_joint_row_error2, Irina_ind_row_error1, Irina_ind_row_error2,
          slide_row_error1, slide_row_error2, slide_row_joint_error1, slide_row_joint_error2, slide_row_ind_error1, slide_row_ind_error2,
          jive_row_error1, jive_row_error2, jive_row_joint_error1, jive_row_joint_error2, jive_row_ind_error1, jive_row_ind_error2,
          ajive_row_error1, ajive_row_error2, ajive_row_joint_error1, ajive_row_joint_error2, ajive_row_ind_error1, ajive_row_ind_error2)

y_col = c(my_signal_error1, my_signal_error2, my_joint_col_error1, my_joint_col_error2, my_ind_col_error1, my_ind_col_error2,
          Irina_signal_error1, Irina_signal_error2, Irina_joint_col_error1, Irina_joint_col_error2, Irina_ind_col_error1, Irina_ind_col_error2,
          slide_col_error1, slide_col_error2, slide_col_joint_error1, slide_col_joint_error2, slide_col_ind_error1, slide_col_ind_error2,
          jive_col_error1, jive_col_error2, jive_col_joint_error1, jive_col_joint_error2, jive_col_ind_error1, jive_col_ind_error2,
          ajive_col_error1, ajive_col_error2, ajive_col_joint_error1, ajive_col_joint_error2, ajive_col_ind_error1, ajive_col_ind_error2)

method = as.factor(c(rep("DMMD",6*nrep), rep("DMMD-i",6*nrep), rep("SLIDE",6*nrep), rep("JIVE",6*nrep), rep("AJIVE",6*nrep)))

section = as.factor(rep(c(rep("1st Signal", nrep), rep("2nd Signal", nrep), 
                          rep("1st Joint", nrep), rep("2nd Joint", nrep),
                          rep("1st Individual", nrep), rep("2nd Individual", nrep)),5))

row_mat = data.frame(Error = y_row, Method = method, Section = section)
col_mat = data.frame(Error = y_col, Method = method, Section = section)

# Version 1
gg_row <- ggplot(row_mat, aes(x=Method, y=Error))+ geom_boxplot(aes(fill=Method))+ facet_wrap(~Section)+ theme_bw() + ylab("Relative Error") + theme(strip.background=element_rect(fill="black")) + theme(strip.text=element_text(color="white", face="bold",size = 25), text=element_text(size = 25))
gg_row <- gg_row + scale_fill_manual(values = c("#ffffcc", "#a1dab4", "#41b6c4", "#2c7fb8","#253494"))
gg_row <- gg_row + theme(axis.text.x = element_text(angle = 45, size = 25, vjust = 1, hjust = 1)) + guides(fill="none") 
gg_row

# Version 2
gg_row <- ggplot(row_mat, aes(x=Method, y=Error))+ geom_boxplot(aes(fill=Method))+ facet_wrap(~Section, scales = "free")+ theme_bw() + ylab("Relative Error") + theme(strip.background=element_rect(fill="black")) + theme(strip.text=element_text(color="white", face="bold",size = 25), text=element_text(size = 25))
gg_row <- gg_row + theme(axis.title.x=element_blank(), axis.ticks.x=element_blank(), axis.text.x=element_blank(),  legend.position = "bottom") 
gg_row <- gg_row + scale_fill_manual(values = c("#ffffcc", "#a1dab4", "#41b6c4", "#2c7fb8","#253494"))
gg_row

# Version 1
gg_col <- ggplot(col_mat, aes(x=Method, y=Error)) + geom_boxplot(aes(fill=Method)) + facet_wrap(~Section) + theme_bw() + ylab("Relative Error") + theme(strip.background=element_rect(fill="black")) + theme(strip.text=element_text(color="white", face="bold",size = 25), text=element_text(size = 25))
gg_col <- gg_col + theme( axis.text = element_text(color = 'black', face="bold", size = 19))
gg_col <- gg_col + theme(axis.title.x=element_blank(), axis.ticks.x=element_blank(), legend.position = "none") 
gg_col <- gg_col + scale_fill_manual(values = c("#ffffcc", "#a1dab4", "#41b6c4", "#2c7fb8","#253494"))
gg_col

# Version 2
gg_col <- ggplot(col_mat, aes(x=Method, y=Error)) + geom_boxplot(aes(fill=Method)) + facet_wrap(~Section, scales = "free") + theme_bw() + ylab("Relative Error") + theme(strip.background=element_rect(fill="black")) + theme(strip.text=element_text(color="white", face="bold",size = 25), text=element_text(size = 25))
gg_col <- gg_col + theme(axis.title.x=element_blank(), axis.ticks.x=element_blank(), axis.text.x=element_blank(),  legend.position = "bottom") 
gg_col <- gg_col + scale_fill_manual(values = c("#ffffcc", "#a1dab4", "#41b6c4", "#2c7fb8","#253494"))
gg_col

# Draw boxplot
fig.path = "Simulations/Signal_Identification_Setting5/Figures/Draft/"
pdf(file = paste(fig.path,"Signal Identification for Row Decomposition_All_Setting5.pdf",sep=""), width = 12, height = 9)
print(gg_row)
dev.off()

pdf(file = paste(fig.path,"Signal Identification for Column Decomposition_All_Setting5.pdf",sep=""), width = 12, height = 9)
print(gg_col)
dev.off()
