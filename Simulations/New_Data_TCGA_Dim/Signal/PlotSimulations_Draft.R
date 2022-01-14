rm(list = ls())
load("Simulations/AdditionalSimulations/New_Data_TCGA_Dim/Signal/SLIDEoutput.RData")
output_slide = output
rm(output)
load("Simulations/AdditionalSimulations/New_Data_TCGA_Dim/Signal/AJIVE_output.RData.")
output_ajive = output
rm(output)
load("Simulations/AdditionalSimulations/New_Data_TCGA_Dim/Signal/Irina_output.RData.")
load("Simulations/AdditionalSimulations/New_Data_TCGA_Dim/Signal/output.RData.")
load("Simulations/AdditionalSimulations/New_Data_TCGA_Dim/Data/Data1.RData")
load("Simulations/AdditionalSimulations/New_Data_TCGA_Dim/Data/Data2.RData")
# Load my Fnorm function
source("Simulations/MyFunction/Preliminary_Functions.R")

n = 88
p = 736
nrep = 100
r1 = 8
r2 = 6
joint_rank_col = rep(0,nrep)
joint_rank_row = rep(2,nrep)

library(ggplot2)

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

gg_row <- ggplot(row_mat, aes(x=Method, y=Error))
gg_row <- gg_row + geom_boxplot(aes(color=Method))
gg_row <- gg_row + facet_wrap(~Section)
gg_row <- gg_row + theme_bw()
gg_row <- gg_row + ylab("Relative Error")
gg_row <- gg_row + theme(strip.background=element_rect(fill="black"))
gg_row <- gg_row + theme(strip.text=element_text(color="white", face="bold",size = 17))
gg_row <- gg_row + theme(text=element_text(size = 17))
gg_row <- gg_row + theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank()) 
gg_row <- gg_row + scale_colour_manual(values = c("blue", "black", "maroon", "red","green"))
gg_row

gg_col <- ggplot(col_mat, aes(x=Method, y=Error))
gg_col <- gg_col + geom_boxplot(aes(color=Method))
gg_col <- gg_col + facet_wrap(~Section)
gg_col <- gg_col + theme_bw()
gg_col <- gg_col + ylab("Relative Error")
gg_col <- gg_col + theme(strip.background=element_rect(fill="black"))
gg_col <- gg_col + theme(strip.text=element_text(color="white", face="bold",size = 17))
gg_col <- gg_col + theme(text=element_text(size = 17))
gg_col <- gg_col + theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank()) 
#gg_col_jive <- gg_col_jive + ggtitle("Signal Error for Column Decomposition Setting 4. DMMD vs JIVE") + ylab("Relative Error")
gg_col <- gg_col + scale_colour_manual(values = c("blue", "black", "maroon", "red","green"))
gg_col

# Draw boxplot
fig.path = "Simulations/AdditionalSimulations/New_Data_TCGA_Dim/Signal/Figures/"
pdf(file = paste(fig.path,"Signal Identification for Row Decomposition_TCGASetting.pdf",sep=""), width = 16, height = 9)
print(gg_row)
dev.off()

pdf(file = paste(fig.path,"Signal Identification for Column Decomposition_TCGASetting.pdf",sep=""), width = 16, height = 9)
print(gg_col)
dev.off()