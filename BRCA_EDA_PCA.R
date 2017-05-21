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
#dataSE<-get(load("BRCA_Illumina_HiSeqnew_PreprocessedData_wo_batch_GC_010_updatedSE_allTypes_allStages_Female_cancerOnly_2017-04-02_16-22.rda"))
#samples.matrix<-get(load("BRCA_Illumina_HiSeqnew_PreprocessedData_wo_batch_GC_010_sampleMatrix_allTypes_allStages_Female_cancerOnly_2017-04-02_16-23.rda"))
#dim(samples.matrix)
#dim(dataSE)



#### rename morphology
samples.matrix<-renameMorph(samples.matrix)

####  add clinical data
samples.matrix<-addClinData(samples.matrix)

#### extract data only for autophagy genes
dataSE <- getAutophagyGenes(dataSE)


# if used addPAM50 DO THIS:
#PAMNPout<-addPAM50andNormal(samples.matrix) #514
PAMNPout<-addXtraPAMandNormal(samples.matrix)# (807 +104 -39) + 112 = 984
samples.matrix<-PAMNPout$samples.matrix 
dim(samples.matrix)

sampleTokeepSE <- PAMNPout$samplesToKeep
dataSE<- dataSE[,colnames(dataSE) %in% sampleTokeepSE]
dim(dataSE) 

# if used addOnlyTopTSS DO THIS:
samples.matrix<-addOnlyTopTSS(samples.matrix)
dim(samples.matrix)
dataSE<- dataSE[,c(samples.matrix$barcode)]
dim(dataSE)






###### renaming the largest morphology group to NA for testing -DO NOT RUN############

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

customCol=c("#4e8fec","#e7ad8a", "#8f2ccf",
            "#4cee30","#d7b51c","#c96629","#61789b","#36d1c8",
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
p <- d %>%  tbl_df %>%  gather(key="PC", value="score", contains("PC")) %>% filter(PC %in% paste0("PC", 1:6))

ggplot(p, aes(x=year_diagnosed, y=score, fill=year_diagnosed)) + 
  geom_boxplot(outlier.colour=NaN) + 
  geom_jitter(aes(color=year_diagnosed),alpha=0.3) + 
  facet_wrap(~PC) + 
  #scale_fill_manual(values=customCol)+
  theme(legend.position="bottom",
        axis.text.x = element_text(angle = 90, hjust = 1))
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
