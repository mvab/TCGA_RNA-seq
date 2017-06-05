setwd("~/Bioinformatics MSc UCPH/0_MasterThesis/TCGAbiolinks/CBL_scripts/data")


######### data from GSE ######

data<-read.table("GSE62944_06_01_15_TCGA_24_548_Clinical_Variables_9264_Samples.txt",  sep = '\t', fill =T )
dim(data)
View(data[,c(1,4:24)])

data$V2<-NULL
data$V3<-NULL
data<-data[-c(191), ] #there is NA
View(data[,1:10])
dim(data)
cats<-as.vector(data$V1)
cats2<- append(c("full_barcode"), cats[2:315])
rownames(data)<-cats2
data$V1<-NULL
###



barcodes_all<- factor(unlist(unname(data[3,1:9264]))) 
colnames(data) <- barcodes_all

female_BRCA<-as.vector(read.table("female_onlyTandN.txt", header =F))$V1
length(female_BRCA)

dataBRCA<- data[,colnames(data) %in% female_BRCA]
dim(dataBRCA)
#dataBRCA_t<-as.data.frame(t(dataBRCA)) #for testing
View(dataBRCA[,1:7])

#replacing with NAs
to_remove <-c("[Not Available]", "[Not Applicable]", "")
for (i in 1:ncol(dataBRCA)){
  this_col<-dataBRCA[,i]
  this_col[this_col %in% to_remove]<-NA
  dataBRCA[,i]<-this_col
  print(i)  
}
## removing NAa
final <- dataBRCA[rowMeans(is.na(dataBRCA)) <= 0.2, ] #rows that have no more than 20% NAs
dim(final)
View(final[,1:7])

final_t<-as.data.frame(t(final))
View(final_t[1:10,])

personal_metadata<-subset(final_t[,c("bcr_patient_barcode", "race", "ethnicity")])
View(personal_metadata)
names(personal_metadata)[1]<-"patient"



#######################

#other data from Kristoffer -> this dataset is more complete -> use this + personal data from the other

brcadata<-read.table("clinical data/brca_clinical_data_Kristoffer.txt", sep = '\t', header = T)
dim(brcadata)
View(brcadata[1:10,1:10])

View(brcadata[1:10,c("sampleID", "AJCC_Stage_nature2012", "PAM50Call_RNAseq", "PAM50_mRNA_nature2012")])

brcadata_cancer <- brcadata[brcadata$sample_type == 'Primary Tumor',]
rownames(brcadata_cancer)<-brcadata_cancer$X_PATIENT

dataBRCA_other<- brcadata_cancer[rownames(brcadata_cancer) %in% female_BRCA,]
dim(dataBRCA_other)
View(dataBRCA_other[1:10,])
dataBRCA_other

#replacing with NAs


to_remove <-c("[Not Available]", "[Not Applicable]", "[Not Evaluated]", "")
for (i in 1:ncol(dataBRCA_other)){
  this_col<-dataBRCA_other[,i]
  this_col[this_col %in% to_remove]<-NA
  dataBRCA_other[,i]<-this_col
  print(i)  
}


final_other <- dataBRCA_other[, colMeans(is.na(dataBRCA_other)) <= 0.2 ] #cols that have no more than 20% NAs
dim(final_other) #1082 x 75
View(final_other)

#things to consider
table(final_other$menopause_status)
table(final_other$age_at_initial_pathologic_diagnosis)
table(final_other$form_completion_date) #can't use!
table(final_other$histological_type)
table(final_other$icd_o_3_histology)
table(final_other$pathologic_stage)
table(final_other$year_of_initial_pathologic_diagnosis)
table(final_other$tissue_source_site)
table(final_other$pathologic_M)
table(final_other$pathologic_T)
table(final_other$pathologic_N)


########### MENOPAUSE
#############
#[1] Pre (<6 months since LMP AND no prior bilateral ovariectomy AND not on estrogen replacement)
#[2] Post (prior bilateral ovariectomy OR >12 mo since LMP with no prior hysterectomy)           
#[3] <NA>                                                                                        
#[4] Indeterminate (neither Pre or Postmenopausal)                                               
#[5] Peri (6-12 months since last menstrual period)  
##########

menopause<- factor(final_other$menopause_status)
menop_renamed <- vector(mode="character", length=0)
for (i in 1:length(menopause)){
  if( is.na(menopause[i]) ){
    menop_renamed[i]<-"Unknown"
  }else if ( startsWith(as.character(menopause[i]), "Pre")){
    menop_renamed[i]<-"Pre menopause"
  }else if ( startsWith(as.character(menopause[i]), "Post")){
    menop_renamed[i]<-"Post menopause"
  }else if ( startsWith(as.character(menopause[i]), "Peri")){
    menop_renamed[i]<-"Peri menopause"
  }else if ( startsWith(as.character(menopause[i]), "Indeterminate")){
    menop_renamed[i]<-"Unknown"
  }
}
final_other$menopause <- menop_renamed

# AGE #
age_diag<- as.vector((final_other$age_at_initial_pathologic_diagnosis))
ageGroups <- vector(mode="character", length=0)
for (i in 1:length(age_diag)){
  print(age_diag[i])
  
  if (is.na(age_diag[i])){
    print ("error")
    ageGroups[i]<-  NA 
  } else if (age_diag[i] < 40) {
    ageGroups[i]<- "< 40"
  } else if ((age_diag[i] >= 40) & (age_diag[i] < 50)){
    ageGroups[i]<- "40-49"
  } else if  ((age_diag[i] >= 50) & (age_diag[i] < 60)){
    ageGroups[i]<- "50-59"
  } else if ((age_diag[i] >= 60) & (age_diag[i] < 70)){
    ageGroups[i] <- "60-69"
  } else if (age_diag[i] >= 70) {
    ageGroups[i]<- "70+"  
  }  
}
final_other$ageGroups<-ageGroups


# MORPHOLOGY #
table(final_other$histological_type)
histtype<- as.vector((final_other$histological_type))
for (i in 1:length(histtype)){
  if (!is.na(histtype[i])){
      if(histtype[i]=="Other  specify"){
        histtype[i]<-"Other"}
      else if (histtype[i]=="Mixed Histology (please specify)"){
        histtype[i]<-"Mixed"}
  }
}
final_other$morphology<-histtype


# STAGE #
stage<-final_other$pathologic_stage

for (i in 1:length(stage)){
  if( is.na(stage[i]) ){
    stage[i]<-"Stage X"
  }else if ( stage[i] == "[Discrepancy]"){
    stage[i]<-"Stage X"
  }else if ( stage[i] == "Stage Tis"){
    stage[i]<-"Stage I"
  }
}
final_other$pathologic_stage <- factor(stage)


# pathology #

# M #

metastasis<- factor(final_other$pathologic_M)
for (i in 1:length(metastasis)){
  if( is.na(metastasis[i]) ){
    metastasis[i]<-"MX"
  }else if ( metastasis[i] == "cM0 (i+)"){
    metastasis[i]<-"M0"
  }
}
final_other$metastasis <- factor(metastasis)

# N #
nodes<- factor(final_other$pathologic_N)
for (i in 1:length(nodes)){
  if( is.na(nodes[i]) ){
    nodes[i]<-"NX"
  }else if ( startsWith(as.character(nodes[i]), "N0")){
    nodes[i]<-"N0"
  }else if ( startsWith(as.character(nodes[i]), "N1")){
    nodes[i]<-"N1"
  }else if ( startsWith(as.character(nodes[i]), "N2")){
    nodes[i]<-"N2"
  }else if ( startsWith(as.character(nodes[i]), "N3")){
    nodes[i]<-"N3"
  }else if ( startsWith(as.character(nodes[i]), "NX")){
    nodes[i]<-"NX"
  }else{
    print (nodes[i])
  }
}
final_other$nodes <- factor(nodes)

# T #
tumour_size<- factor(final_other$pathologic_T)
for (i in 1:length(tumour_size)){
  if( is.na(tumour_size[i]) ){
    tumour_size[i]<-"TX"
  }else if ( startsWith(as.character(tumour_size[i]), "T1")){
    tumour_size[i]<-"T1"
  }else if ( startsWith(as.character(tumour_size[i]), "T2")){
    tumour_size[i]<-"T2"
  }else if ( startsWith(as.character(tumour_size[i]), "T3")){
    tumour_size[i]<-"T3"
  }else if ( startsWith(as.character(tumour_size[i]), "T4")){
    tumour_size[i]<-"T4"
  }else if ( startsWith(as.character(tumour_size[i]), "TX")){
    tumour_size[i]<-"TX"
  }else{
    print (tumour_size[i])
  }
}
final_other$tumour_size <- factor(tumour_size)


# get year data from Elenas data
clindata_elena<-get(load("clinData_prepared.rda"))
moredata<-subset(clindata_elena[, c("patient", "year_diagnosed")])
for (i in 1:length(moredata$year_diagnosed)){
    if (moredata$year_diagnosed[i] == 0){
      moredata$year_diagnosed[i] <- NA
    }else if (moredata$year_diagnosed[i] != 0){
      moredata$year_diagnosed[i] <- as.character(moredata$year_diagnosed[i])
    }
}

 

## new metadata

metadata<-subset(final_other[,c("X_PATIENT","menopause", "ageGroups",  "morphology","icd_o_3_histology",
                                "pathologic_stage", "tissue_source_site", "metastasis", "nodes", "tumour_size")])
View(metadata)
names(metadata)<- c("patient","menopause", "ageGroups", "morphology", "histology", "substage", "tissue source site","metastasis", "nodes", "tumour_size")

# adding personal data to main metadata
metadata_full<-merge(metadata, personal_metadata, by="patient", all.x=TRUE)
metadata_full<-merge(metadata_full, moredata, by="patient", all.x=TRUE)
View(metadata_full)


to_remove <-c("[Not Available]", "[Not Applicable]", "[Not Evaluated]", "")
for (i in 1:ncol(metadata_full)){
  this_col<-metadata_full[,i]
  this_col[this_col %in% to_remove]<-NA
  metadata_full[,i]<-this_col
  print(i)  
}


save(metadata_full, file="cleaned_metadata.rda")
