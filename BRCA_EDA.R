#! /usr/bin/env Rscript
rm(list=ls(all=TRUE))

#library(TCGAbiolinks)
library(SummarizedExperiment)
library(plyr)
library(edgeR)

library(ggplot2)
library(RColorBrewer)

library(limma)

suppressMessages(library(magrittr))
suppressMessages(library(tidyr))
suppressMessages(library(dplyr))


setwd("~/Bioinformatics MSc UCPH/0_MasterThesis/TCGAbiolinks/CBL_scripts/data")


#femaly only T/N data upload
#dataSE<-get(load("brcaExp_PreprocessedData_wo_batch_updatedSE_allFemale.rda"))
#samples.matrixA<-get(load("brcaExp_PreprocessedData_wo_batch_sampleMatrix_allFemale.rda"))
#dim(samples.matrixA)


## most recent all types and stages cancer and normal
dataSE<-get(load("BRCA_Illumina_HiSeqnew_updatedSE_allFemale_PreprocessedData_wo_batch_GC_010_updatedSE_allTypes_allStages_Female.rda"))
samples.matrix<-get(load("BRCA_Illumina_HiSeqnew_updatedSE_allFemale_PreprocessedData_wo_batch_GC_010_sampleMatrix_allTypes_allStages_Female.rda"))
dim(samples.matrix)
dim(dataSE)



#trying cancer only with subtypes (1081)
dataSE<-get(load("BRCA_Illumina_HiSeqnew_PreprocessedData_wo_batch_GC_010_updatedSE_allTypes_allStages_Female_cancerOnly_2017-04-02_16-22.rda"))
samples.matrix<-get(load("BRCA_Illumina_HiSeqnew_PreprocessedData_wo_batch_GC_010_sampleMatrix_allTypes_allStages_Female_cancerOnly_2017-04-02_16-23.rda"))
dim(samples.matrix)
dim(dataSE)



#### extract data only for autophagy genes

getAutophagyGenes <- function(dataSE){
  
  autophagy_genes<- as.vector(read.table("autopahagy_genes.txt", as.is = T, header = FALSE))$V1
  all_genes<-rownames(dataSE)
  shared <- intersect(autophagy_genes,all_genes)
  dataSE<- dataSE[c(shared),]
  
  print(paste0("Total number of genes: ", length(all_genes)))
  print(paste0("Autophagy genes: ", length(shared)))
  
  return(dataSE)
} 

dataSE <- getAutophagyGenes(dataSE)
dim(dataSE)




########## exploration of different matrices : DO NOT RUN for standard case ########

addClinData<-function(sample.matrix){
  
  clindata<-get(load("clinData_prepared.rda"))
  samples.matrix_clinincal <- merge(samples.matrix, clindata, by="patient", all.x=TRUE) 
  #samples.matrix_clinincal <- samples.matrix_clinincal[order(samples.matrix_clinincal$myorder), ]
  
  return(samples.matrix_clinincal)
}

# only do this if want to add clinical data
samples.matrix<-addClinData(samples.matrix)
dim(samples.matrix)
head(samples.matrix)
dataSE<- dataSE[,c(samples.matrix$barcode)]
dim(dataSE)


addOnlyTopTSS <- function(samples.matrix){
  
  #filtering out patient with neither Positive nor Negative
  #tss_above50samples <-c("BH","A2","E2","A8","D8","E9", "AR", "B6","AC" ) #660 patients
  tss_2 <-c("A2","E9") 
  
  #samples.matrix <- samples.matrix[samples.matrix$tss %in% tss_above50samples,] 
  samples.matrix <- samples.matrix[samples.matrix$tss %in% tss_2,] 
  
  return (samples.matrix)
}   
  # if used addOnlyTopTSS DO THIS:
  samples.matrixTSS<-addOnlyTopTSS(samples.matrix)
  dim(samples.matrixTSS)

  dataSE<- dataSE[,c(samples.matrixTSS$barcode)]
  dim(dataSE)


addBRCAReceptorStatus <- function(samples.matrix){

    # get information on subtype/pation details
    dataSubt <- TCGAquery_subtype(tumor = "BRCA") 
    
    #add extra patient information, but some have NAs; if FALSE- no NAs but lose ~200 patients #~735
    subdiagnosis <- merge(samples.matrix, dataSubt, by="patient", all.x=TRUE) 
    subdiagnosis <- subdiagnosis[order(subdiagnosis$myorder), ]
    
    # we only have female samples, so NA in gender can be replaced with FEMALE
    subdiagnosis[c("Gender")][is.na(subdiagnosis[c("Gender")])] <- "FEMALE"
    
    # note:  the are many thing not relevant - remove them and then do complete cases
    
    # testing how many pations how no NAs for ER PR HER cols
    complete_cases_receptors<- subdiagnosis[complete.cases(subdiagnosis[, c(1:12, 15:17)]),]# 704 patients 
    
    #filtering out patient with neither Positive nor Negative
    other_recp_status <-c("Indeterminate", "Not Performed", "Performed but Not Available", "Equivocal", "Not Available")
    complete_cases_receptors$ER.Status[complete_cases_receptors$ER.Status %in% other_recp_status] <- NA
    complete_cases_receptors$PR.Status[complete_cases_receptors$PR.Status %in% other_recp_status] <- NA
    complete_cases_receptors$HER2.Final.Status[complete_cases_receptors$HER2.Final.Status %in% other_recp_status] <- NA
    
    complete_cases_receptors<- complete_cases_receptors[complete.cases(complete_cases_receptors[, c(1:12, 15:17)]),]
    complete_cases_receptors<- subset(complete_cases_receptors[, c(1:12, 15:17)]) # 646 patients
    dim(complete_cases_receptors)
    
    
    # replacing Positive and Negative in receptor column with 0 and 1 
    
    complete_cases_receptors$ER.Status <- gsub('Positive', 1, complete_cases_receptors$ER.Status)
    complete_cases_receptors$ER.Status <- gsub('Negative', 0, complete_cases_receptors$ER.Status)
    complete_cases_receptors$PR.Status <- gsub('Positive', 1, complete_cases_receptors$PR.Status)
    complete_cases_receptors$PR.Status <- gsub('Negative', 0, complete_cases_receptors$PR.Status)
    complete_cases_receptors$HER2.Final.Status <- gsub('Positive', 1, complete_cases_receptors$HER2.Final.Status)
    complete_cases_receptors$HER2.Final.Status <- gsub('Negative', 0, complete_cases_receptors$HER2.Final.Status)
    
    receptor_statusNum = as.matrix(as.data.frame(lapply(subset(complete_cases_receptors[13:15]), as.numeric)))
    receptor_status$sum <- rowSums(receptor_statusNum)
    
    # adding subtype classification:
    
    ### Luminal     ###   ER+PR+HER+ or ER+PR+HER- , so 3 or 2
    ### HER2 over   ###   ER-PR-HER+  , so 1
    ### Triple Neg  ###   ER-PR-HER-  , so 0
    
    receptSubtype <- vector(mode="character", length=0)
    for (i in 1:length(receptor_status$sum)){
      if (receptor_status$sum[i] > 1) {
          receptSubtype[i]<-"Luminal"
      } else if (receptor_status$sum[i] == 1) {
          receptSubtype[i]<-"HER2"
      } else {
          receptSubtype[i]<-"TNR"
      }
    }
    receptor_status$receptSubtype <- receptSubtype
    # add receptSubtype to the main sample matrix and remove sepate receptors
    
    complete_cases_receptors<-subset(complete_cases_receptors[, c(1:12)])
    complete_cases_receptors$receptSubtype <- receptSubtype
    
    return (complete_cases_receptors)
}

  samples.matrixRecep<-addBRCAReceptorStatus(samples.matrix)
  # if used addBRCAReceptorStatus DO THIS:
  dataSE<- dataSE[,c(samples.matrixRecep$barcode)]
  dim(dataSE)
  
########################################################  

addXtraandNormal<-function(samples.matrix){
    
    extradata<-get(load("extradata_var_w_old.rda"))
    
    #add extra patient information #578
    subdiagnosis <- merge(samples.matrix, extradata, by="patient", all.x=T) 
    subdiagnosis <- subdiagnosis[order(subdiagnosis$myorder), ]
    
    # get barcodes of samples with PAM50 recorded
    diag <- subdiagnosis[!is.na(subdiagnosis$PAM50new), ] 
    diag <- diag$barcode #578
    
    # get barcodes of patients with paired samples
    
    #get patients with 2 entries == 224
    paired <- samples.matrix[ duplicated(samples.matrix$participant, fromLast=FALSE) | duplicated(samples.matrix$participant, fromLast=TRUE),] 
    paired <- as.character(paired$participant)
    #get normal samples from paired data == 112
    paired <- samples.matrix[samples.matrix$participant %in% paired & samples.matrix$condition == "normal", ]  
    paired <- as.character(paired$barcode)
    
    # concatnate barcodes to keep
    #join NA and paired cancer samples
    to_keep <- unique(sort(c(as.character(diag), as.character(paired)))) #578+112=626 unique (some normals have subtype label)

    # remove samples with no subdtype label and paired tumor sample from dataframe, subtype vector and ID-information.
    samples.matrix <- samples.matrix[samples.matrix$barcode %in% to_keep, ]    
    subdiagnosis <- subdiagnosis[subdiagnosis$barcode %in% to_keep, ]
    
    # make normal samples from pairs "normal" while retaining info on subtype.
    subdiagnosis$PAM50 <- ifelse(subdiagnosis$condition == "normal", "normal", as.character(subdiagnosis$PAM50new))
    
    print ("Keeping PAM50 samples and all normals:")
    print( table(subdiagnosis$PAM50))
    
    return(list(samples.matrix=subdiagnosis, samplesToKeep=to_keep, normal_barcodes = paired))
}  
  
PAMNPout<-addXtraandNormal(samples.matrix)
samples.matrix<-PAMNPout$samples.matrix 
dim(samples.matrix)

sampleTokeepSE <- PAMNPout$samplesToKeep
dataSE<- dataSE[,colnames(dataSE) %in% sampleTokeepSE]
dim(dataSE) 
  

source("../functions.R")  
  
  # if used addPAM50andNormal DO THIS:
PAMNPout<-addPAM50andNormal(samples.matrix)
samples.matrix<-PAMNPout$samples.matrix 
dim(samples.matrix)

sampleTokeepSE <- PAMNPout$samplesToKeep
dataSE<- dataSE[,colnames(dataSE) %in% sampleTokeepSE]
dim(dataSE) 


## renaming differing PAM labels
addAltPAM50<-function(samples.matrix){

    PAM50_upd<-get(load("potentially_wrong_PAM.rda"))
    names(PAM50_upd)[2]<-"PAM50upd"

    samples.matrix <- merge(samples.matrix, PAM50_upd, by="patient", all.x=TRUE) 
    samples.matrix <- samples.matrix[order(samples.matrix$myorder), ]
    samples.matrix$PAM50upd <- ifelse(samples.matrix$condition == "normal", "normal", as.character(samples.matrix$PAM50upd))
    
    for (i in 1:length(samples.matrix$patient)){
      if (is.na(samples.matrix$PAM50upd[i])){
        samples.matrix$PAM50upd[i]<-samples.matrix$PAM50[i]
      }
    }
    return(samples.matrix)
}

samples.matrix<-addAltPAM50(samples.matrix)
names(samples.matrix)

###### renaming the largest morphology group to NA for testing############

#first rename NA to normal

samples.matrix$tumourTypes <- as.character(samples.matrix$tumourTypes)
samples.matrix$tumourTypes[is.na(samples.matrix$tumourTypes)] <- 'normal'
samples.matrix$tumourTypes <- as.factor(samples.matrix$tumourTypes)

#then rename major group to NA so that it is not displayed
samples.matrix$tumourTypes <- as.character(samples.matrix$tumourTypes)
samples.matrix$tumourTypes[samples.matrix$tumourTypes == "female_85003"] <- NA
samples.matrix$tumourTypes <- as.factor(samples.matrix$tumourTypes)

##############



  
###################### EDA for (full) dataset #########################
cbPalette <- c( "#E69F00", "#0072B2", "#F0E442","#D55E00", "#009E73", "#CC79A7",   "#56B4E9", "black")

dim(dataSE)    
dge <- DGEList(dataSE)
dge <- calcNormFactors(object=dge, method="TMM")
EM <- cpm(x=dge, log=TRUE)

# Perfrom PCA
pca <- prcomp(x=t(EM), scale=TRUE, center=TRUE)

# Inspect components
summary(pca)
#type                               
qplot(data=as.data.frame(pca$x), x=PC1, y=PC2, geom=c("point"), color=samples.matrix$PAM50)

#PAM50 status
ggplot(data=as.data.frame(pca$x),aes(x=PC1,y=PC2,col=samples.matrix$PAM50upd))+ #, shape =samples.matrix$tumourType
  geom_point(size=2.5,alpha=0.9)+ #Size and alpha just for fun
  #scale_colour_brewer(palette = "Set2")+ #your colors 
  scale_colour_manual(values=cbPalette)+
  theme_classic()+
  theme(legend.position="bottom",
        legend.title=element_blank())

# tumour types
ggplot(data=as.data.frame(pca$x),aes(x=PC3,y=PC2,col=samples.matrix$tumourTypes))+ #, shape =samples.matrix$tumourType
  geom_point(size=2.5,alpha=0.9)+ #Size and alpha just for fun
  #scale_colour_brewer(palette = "Set2")+ #your colors here
  scale_colour_manual(values=cbPalette)+
  theme_classic()+
  theme(legend.position="bottom",
        legend.title=element_blank())
  


# tumour stages
ggplot(data=as.data.frame(pca$x),aes(x=PC1,y=PC3,col=samples.matrix$tumourStages))+ #, shape =samples.matrix$tumourType
  geom_point(size=2.5,alpha=0.9)+ #Size and alpha just for fun
  #scale_colour_brewer(palette = "Set2")+ #your colors here
  scale_colour_manual(values=cbPalette)+
  theme_classic()+
  theme(legend.position="bottom",
        legend.title=element_blank())

#age group
ggplot(data=as.data.frame(pca$x),aes(x=PC1,y=PC2,col=samples.matrix$ageGroup))+ #, shape =samples.matrix$tumourType
  geom_point(size=2.5,alpha=0.9)+ #Size and alpha just for fun
  scale_colour_brewer(palette = "Set2")+ #your colors here
  theme_classic()+
  theme(legend.position="bottom",
        legend.title=element_blank())

#year dignosed
ggplot(data=as.data.frame(pca$x),aes(x=PC1,y=PC3,col=samples.matrix$yearGroup))+ #, shape =samples.matrix$tumourType
  geom_point(size=2.5,alpha=0.9)+ #Size and alpha just for fun
  scale_colour_brewer(palette = "Set2")+ #your colors here
  theme_classic()+
  theme(legend.position="bottom",
        legend.title=element_blank())

#race
ggplot(data=as.data.frame(pca$x),aes(x=PC1,y=PC3,col=samples.matrix$race))+ #, shape =samples.matrix$tumourType
  geom_point(size=2.5,alpha=0.9)+ #Size and alpha just for fun
  scale_colour_brewer(palette = "Set2")+ #your colors here
  theme_classic()+
  theme(legend.position="bottom",
        legend.title=element_blank())


#### exploring PCs

# T/N
d <- data.frame(samples.matrix, pca$x)
p <- d %>%  tbl_df %>%  gather(key="PC", value="score", contains("PC")) %>% filter(PC %in% paste0("PC", 1:9))

ggplot(p, aes(x=condition, y=score, fill=condition)) + 
  geom_boxplot(outlier.colour=NaN) + 
  #geom_jitter(alpha=0.5) + 
  facet_wrap(~PC) + 
  scale_fill_brewer(palette="Set2")

#PAM50
d <- data.frame(samples.matrix, pca$x)
p <- d %>%  tbl_df %>%  gather(key="PC", value="score", contains("PC")) %>% filter(PC %in% paste0("PC", 1:9))

ggplot(p, aes(x=PAM50, y=score, fill=PAM50)) + 
  geom_boxplot(outlier.colour=NaN) + 
  #geom_jitter(alpha=0.5) + 
  facet_wrap(~PC) + 
  scale_fill_manual(values=cbPalette)
  #scale_fill_brewer(palette="Set2")


#tumour types 
d <- data.frame(samples.matrix, pca$x)
p <- d %>%  tbl_df %>%  gather(key="PC", value="score", contains("PC")) %>% filter(PC %in% paste0("PC", 1:9))

ggplot(p, aes(x=tumourTypes, y=score, fill=tumourTypes)) + 
  geom_boxplot(outlier.colour=NaN) + 
  #geom_jitter(alpha=0.5) + 
  facet_wrap(~PC) + 
  scale_fill_manual(values=cbPalette)
  #scale_fill_brewer(palette="Set2")

#tumour stages
d <- data.frame(samples.matrix, pca$x)
p <- d %>%  tbl_df %>%  gather(key="PC", value="score", contains("PC")) %>% filter(PC %in% paste0("PC", 1:9))

ggplot(p, aes(x=tumourStages, y=score, fill=tumourStages)) + 
  geom_boxplot(outlier.colour=NaN) + 
  #geom_jitter(alpha=0.5) + 
  facet_wrap(~PC) + 
  scale_fill_manual(values=cbPalette)
  #scale_fill_brewer(palette="Set2")


#age group
d <- data.frame(samples.matrix, pca$x)
p <- d %>%  tbl_df %>%  gather(key="PC", value="score", contains("PC")) %>% filter(PC %in% paste0("PC", 1:9))

ggplot(p, aes(x=ageGroups, y=score, fill=ageGroups)) + 
  geom_boxplot(outlier.colour=NaN) + 
  #geom_jitter(alpha=0.5) + 
  facet_wrap(~PC) + 
  scale_fill_brewer(palette="Set2")

#year group
d <- data.frame(samples.matrix, pca$x)
p <- d %>%  tbl_df %>%  gather(key="PC", value="score", contains("PC")) %>% filter(PC %in% paste0("PC", 1:9))

ggplot(p, aes(x=yearGroups, y=score, fill=yearGroups)) + 
  geom_boxplot(outlier.colour=NaN) + 
  #geom_jitter(alpha=0.5) + 
  facet_wrap(~PC) + 
  scale_fill_brewer(palette="Set2")

#race
d <- data.frame(samples.matrix, pca$x)
p <- d %>%  tbl_df %>%  gather(key="PC", value="score", contains("PC")) %>% filter(PC %in% paste0("PC", 1:9))

ggplot(p, aes(x=race, y=score, fill=race)) + 
  geom_boxplot(outlier.colour=NaN) + 
  #geom_jitter(alpha=0.5) + 
  facet_wrap(~PC) + 
  scale_fill_brewer(palette="Set2")



stop()
################## EDA for subtypes ###################

dge <- DGEList(dataSE)
dge <- calcNormFactors(object=dge, method="TMM")
EM <- cpm(x=dge, log=TRUE)


# Perfrom PCA
pca <- prcomp(x=t(EM), scale=TRUE, center=TRUE)

# Inspect components
summary(pca)




fortyColours<-c("#3579b4","#c8c049","#8996da","#ee5345","#c84297","#43ea3b","#50a376","#281340","#6e5c41","#94f5d2","#fd0d32","#f19832","#b1f555","#d727b1","#f27456","#4bfe9f","#61789b","#2896be","#db1453","#c7a233","#d9a5c8","#1e785f","#3183e5","#82117f","#e5cbb0","#2dc194","#8f2ccf","#4e8fec","#e7ad8a","#234220","#4cee30","#d7b51c","#c96629","#472134","#36d1c8","#9f6f63","#ac8d3c","#a63dbd","#1db9d9","#10c399")
thirteenColours<-c("#3579b4","#c8c049","#8996da","#ee5345","#c84297","#43ea3b","#50a376","#281340","#6e5c41","#fd0d32","#f19832","#94f5d2","#b1f555")
#TSS
qplot(data=as.data.frame(pca$x), x=PC1, y=PC2, geom=c("point"), color=out$tss)

qplot(data=as.data.frame(pca$x), x=PC1, y=PC2, geom=c("text"), color=out$tss, label = out$tss) +scale_color_manual(values=fortyColours)

ggplot(data=as.data.frame(pca$x),aes(x=PC1,y=PC2,col=out$tss))+
  geom_point(size=2,alpha=0.5)+ #Size and alpha just for fun
  scale_color_manual(values = fortyColours)+ #your colors here
  theme_classic()

#TSS after reduction 
qplot(data=as.data.frame(pca$x), x=PC1, y=PC3, geom=c("point"), color=samples.matrixTSS$tss)
ggplot(data=as.data.frame(pca$x),aes(x=PC1,y=PC5,col=samples.matrixTSS$tss))+
  geom_point(size=2,alpha=0.9)+ #Size and alpha just for fun
  #scale_color_manual(values = thirteenColours)+ #your colors here
  theme_classic()

#type                               #or 5
qplot(data=as.data.frame(pca$x), x=PC1, y=PC3, geom=c("point"), color=out$tumourType)

#portion (nothin)
qplot(data=as.data.frame(pca$x), x=PC3, y=PC2, geom=c("point"), color=out$portion)

# receptor
ggplot(data=as.data.frame(pca$x),aes(x=PC1,y=PC3,col=samples.matrixRecep$receptSubtype))+
  geom_point(size=2)+ #Size and alpha just for fun
  scale_fill_brewer(palette="Set3")+ #your colors here
  theme(legend.position="bottom")



ggplot(data=as.data.frame(pca$x),aes(x=PC3,y=PC7,col=out$portion))+
  geom_point(size=2,alpha=0.7)+ #Size and alpha just for fun
  scale_color_manual(values = thirteenColours)+ #your colors here
  theme_classic()

ggplot(data=as.data.frame(pca$x),aes(x=PC1,y=PC3,col=out$plate))+
  geom_point(size=2,alpha=0.7)+ #Size and alpha just for fun
  scale_color_manual(values = fortyColours)+ #your colors here
  theme_classic()



#### exploring PCs

suppressMessages(library(magrittr))
suppressMessages(library(tidyr))
suppressMessages(library(dplyr))

d <- data.frame(samples.matrixTSS, pca$x)
p <- d %>%  tbl_df %>%  gather(key="PC", value="score", contains("PC")) %>% filter(PC %in% paste0("PC", 1:9))

ggplot(p, aes(x=tumourType, y=score, fill=tumourType)) + 
  geom_boxplot(outlier.colour=NaN) + 
  #geom_jitter(alpha=0.5) + 
  facet_wrap(~PC) + 
  scale_fill_brewer(palette="Set2")

ggplot(p, aes(x=portion, y=score, fill=portion)) + 
  geom_boxplot(outlier.colour=NaN) + 
  #geom_jitter(alpha=0.05) +   #uncomment this to see most dense group
  facet_wrap(~PC) + 
  scale_color_manual(values = thirteenColours)
  #scale_fill_brewer(palette="Set3")

ggplot(p, aes(x=tss, y=score, fill=tss)) + 
  geom_boxplot(outlier.colour=NaN) + 
  #geom_jitter(alpha=0.05) +   #uncomment this to see most dense group
  facet_wrap(~PC) 

ggplot(p, aes(x=sample, y=score, fill=sample)) + 
  geom_boxplot(outlier.colour=NaN) + 
  #geom_jitter(alpha=0.5) + 
  facet_wrap(~PC) + 
  scale_fill_brewer(palette="Set3")

ggplot(p, aes(x=plate, y=score, fill=plate)) + 
  geom_boxplot(outlier.colour=NaN) + 
  #geom_jitter(alpha=0.5) + 
  facet_wrap(~PC) + 
  scale_fill_brewer(palette="Set3")

ggplot(p, aes(x=receptSubtype, y=score, fill=receptSubtype)) + 
  geom_boxplot(outlier.colour=NaN) + 
  #geom_jitter(alpha=0.5) + 
  facet_wrap(~PC) + 
  scale_fill_brewer(palette="Set3")

ggplot(p, aes(x=tss, y=score, fill=tss)) + 
  geom_boxplot(outlier.colour=NaN) + 
  #geom_jitter(alpha=0.5) + 
  facet_wrap(~PC) 
