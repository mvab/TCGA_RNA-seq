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
source("../functions.R") 


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



#### rename morphology
samples.matrix<-renameMorph(samples.matrix)

####  add clinical data
samples.matrix<-addClinData(samples.matrix)

#### extract data only for autophagy genes
dataSE <- getAutophagyGenes(dataSE)


# if used addOnlyTopTSS DO THIS:
samples.matrix<-addOnlyTopTSS(samples.matrix)
dim(samples.matrix)
dataSE<- dataSE[,c(samples.matrix$barcode)]
dim(dataSE)

# if used addPAM50 DO THIS:

#PAMNPout<-addPAM50andNormal(samples.matrix) #514
PAMNPout<-addXtraPAMandNormal(samples.matrix)# (807 +104 -39) + 112 = 984
samples.matrix<-PAMNPout$samples.matrix 
dim(samples.matrix)

sampleTokeepSE <- PAMNPout$samplesToKeep
dataSE<- dataSE[,colnames(dataSE) %in% sampleTokeepSE]
dim(dataSE) 



############# ad hoc testing 
################################
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


  





## renaming differing PAM labels (514, 39 renamed)
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
cbPalette <- c( "#E69F00", "#0072B2", "#F0E442","#D55E00",  "#009E73", "#CC79A7",   "#56B4E9", "black")
cbPaletteM <- c(  "#F0E442", "#D55E00","#E69F00", "#56B4E9", "#CC79A7",  "#0072B2", "#009E73",  "black","white")#
cbPaletteOther <- c(  "#F0E442", "#56B4E9", "#CC79A7", "#D55E00","#E69F00",  "#0072B2", "#009E73",  "black")
fortyColours<-c("#3579b4","#c8c049","#8996da","#ee5345","#c84297","#43ea3b","#50a376","#281340","#6e5c41","#94f5d2","#fd0d32","#f19832","#b1f555","#d727b1","#f27456","#4bfe9f","#61789b","#2896be","#db1453","#c7a233","#d9a5c8","#1e785f","#3183e5","#82117f","#e5cbb0","#2dc194","#8f2ccf","#4e8fec","#e7ad8a","#234220","#4cee30","#d7b51c","#c96629","#472134","#36d1c8","#9f6f63","#ac8d3c","#a63dbd","#1db9d9","#10c399")


dim(dataSE)    
dge <- DGEList(dataSE)
dge <- calcNormFactors(object=dge, method="TMM")
EM <- cpm(x=dge, log=TRUE)

# Perfrom PCA
pca <- prcomp(x=t(EM), scale=TRUE, center=TRUE)

# Inspect components
summary(pca)
#type                               
qplot(data=as.data.frame(pca$x), x=PC1, y=PC2, geom=c("point"), color=samples.matrix$condition)

#PAM50 status
ggplot(data=as.data.frame(pca$x),aes(x=PC1,y=PC2,col=samples.matrix$PAM50))+ #, shape =samples.matrix$tumourType
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
scale_colour_manual(values=cbPaletteM)+
theme_classic()+
theme(legend.position="bottom",
legend.title=element_blank())


########### testing other categories ###########
thirteenColours<-c("#3579b4","#c8c049","#8996da","#ee5345","#c84297","#43ea3b","#50a376","#281340","#6e5c41","#fd0d32","#f19832","#94f5d2","white")
fortyColours2<-c("#3579b4","#c8c049","#8996da","#ee5345","#c84297",
                "#43ea3b","#50a376","#281340","#6e5c41","#94f5d2",
                "#94f5d2","#94f5d2","#94f5d2","#d727b1","#f27456",
                "#ac8d3c","#61789b","#2896be","#db1453","#c7a233",
                "#d9a5c8","#1e785f","#3183e5","#82117f","#e5cbb0",
                "#9f6f63","#8f2ccf","#4e8fec","#e7ad8a","#234220",
                "#4cee30","#d7b51c","#c96629","#472134","#36d1c8",
                "#9f6f63","#ac8d3c","#a63dbd","#1db9d9","#10c399")


# "#3579b4", normal
substagesCols=c( "greenyellow", "green1","green4",
                  "deeppink", "deeppink2","deeppink4", 
                  "deepskyblue","deepskyblue2", "deepskyblue4", "dodgerblue4",
                  "orange","white")

names(samples.matrix)

ggplot(data=as.data.frame(pca$x),aes(x=PC1,y=PC2,col=samples.matrix$metastasis))+ #, shape =samples.matrix$tumourType
  geom_point(size=2.5,alpha=0.9)+ #Size and alpha just for fun
  #scale_colour_brewer(palette = "Set2")+ #your colors here
  scale_colour_manual(values=cbPaletteM)+# thirteenColours)+
  theme_classic()+
  theme(legend.position="bottom",
        legend.title=element_blank())

d <- data.frame(samples.matrix, pca$x)
p <- d %>%  tbl_df %>%  gather(key="PC", value="score", contains("PC")) %>% filter(PC %in% paste0("PC", 1:9))

ggplot(p, aes(x=tss, y=score, fill=tss)) + 
  geom_boxplot(outlier.colour=NaN) + 
  #geom_jitter(alpha=0.5) + 
  facet_wrap(~PC) + 
  scale_fill_manual(values=fortyColours)
  #scale_fill_brewer(palette="Set2")


##################################

  
# tumour stages
ggplot(data=as.data.frame(pca$x),aes(x=PC1,y=PC2,col=samples.matrix$tumourStages))+ #, shape =samples.matrix$tumourType
  geom_point(size=2.5,alpha=0.9)+ #Size and alpha just for fun
  #scale_colour_brewer(palette = "Set2")+ #your colors here
  scale_colour_manual(values=cbPaletteM)+
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
  scale_fill_manual(values=cbPaletteM)
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





stop()
