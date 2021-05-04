rm(list = ls())
list.files()
# Specify path to save figures
# Load the result of DMMD and JIVE
load("Application/TCGA/ClinicalInformation/Subtype_Sex_Information.RData")
load('Application/TCGA/2021New/DoubleSTD_DMMD_JIVE_EqualVaraince_Result.RData')

cdata = dcs_cancer
ndata = dcs_normal
colnames(cdata)

# Index change: we need to find which row in the clinical dataset is corresponding 
# to the column (sample) in the cdata/ndata.
# The result will be in a vector form. e.g. if the first element is 18, that means
# the first column in cdata (TCGA.A7.A13E.01A.11R.A12O.13) appears in the 18th column 
# of clinical data.
index_change_vec = rep(NA,dim(cdata)[2])
for (i in 1:dim(cdata)[2]){
  tmp = substr(colnames(cdata)[i],1,12)
  # I need to change all the dots "." to "-" so that the format matches.
  new_tmp = gsub("[.]","-",tmp)
  ID_index = which(new_tmp == small_clinic$Patient.ID)
  # If there is no match, we mark the place as NA.
  if (length(ID_index) == 1){
    index_change_vec[i] = ID_index
  }

}
index_change_vec
# We found three samples in our double matched data do not have clinical information.
# colnames(cdata)[47]
# small_clinic$Patient.ID[520]
save(small_clinic,index_change_vec, file = "Application/TCGA/ClinicalInformation/IndexChange.RData")



