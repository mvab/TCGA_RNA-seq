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


setwd("~/Bioinformatics MSc UCPH/0_MasterThesis/TCGAbiolinks/CBL_scripts/TCGA-BRCA_analysis_methods/data/")
source("../BRCA_functions.R") 


## most recent all types and stages cancer and normal

dataSE<-get(load("BRCA_Illumina_HiSeqnew_updatedSE_allFemale_PreprocessedData_wo_batch_GC_010_withSampleInfo_SE.rda"))
samples.matrix<-get(load("BRCA_Illumina_HiSeqnew_updatedSE_allFemale_PreprocessedData_wo_batch_GC_010_withSampleInfo_SM.rda"))
dim(samples.matrix)
dim(dataSE)


########## PRE-ANALYIS STEPS ################

#### rename morphology
samples.matrix<-renameMorph(samples.matrix)


####  add clinical data
samples.matrix<-addClinData(samples.matrix)


#### add PAM50 annotation to samples.matrix and update data SE
PAMNPout<-addPAM50annotation(samples.matrix)

samples.matrix<-PAMNPout$samples.matrix #add annotation
dim(samples.matrix) 

sampleTokeepSE <- PAMNPout$samplesToKeep #update SE
dataSE<- dataSE[,colnames(dataSE) %in% sampleTokeepSE]
dim(dataSE) 


# exclude mrphology samples Ductal mixed with Others -- it was added at the beginning, 
# but later it was decided not use it in the analysis. Therefore this step removes those samples 
# from dataSE and samples.matrix

removedsamples<-removeDuctalOther(dataSE, samples.matrix)
dataSE<-removedsamples$dataSE
samples.matrix<-removedsamples$samples.matrix
dim(dataSE)
dim(samples.matrix)




#  #  #  #  ONLY RUN THIS BLOCK IF DO PCA ON STAGES ~ excludes unknown stage 
known_stage<-samples.matrix[samples.matrix$tumourStages!='unknown',]$barcode
dataSE<- dataSE[,colnames(dataSE) %in% known_stage]
dim(dataSE)
samples.matrix<-samples.matrix[samples.matrix$tumourStages!='unknown',]
dim(samples.matrix)





## palettes (run all)
tnPalette=c("deeppink2", "dodgerblue2")
cbPalette <- c( "#E69F00", "#0072B2", "#F0E442", "#D55E00",  "#009E73", "#CC79A7",   "#56B4E9", "black")
tumourTypCol = c( "#56B4E9", "#F0E442", "#CC79A7", "#D55E00", "#0072B2", "#009E73", "#E69F00")
tumourStagesCol = c( "#009E73", "#F0E442", "#E69F00", "#CC79A7", "red3")
year_cols<-c("#3579b4", "#c8c049", "#8996da", "#ee5345", "#c84297",
             "#ac8d3c", "#50a376", "#9f6f63", "#36d1c8", "#a63dbd",
             "#8f2ccf", "#4e8fec", "#e7ad8a", "#d727b1", "#f27456",
             "#ac8d3c", "#61789b", "#2896be", "#db1453", "#c7a233",
             "#d9a5c8", "#1e785f", "#3183e5", "#82117f", "#e5cbb0",
             "#d7b51c")
source_cols<-c("#fc8785", "#3cb777", "#575a6c", "#c390ec", "#c5eebd",
               "#4dc9d5", "#d1f6f2", "#e4f90d", "#5c9134", "#4159e2",
               "#f3132d", "#d0abef", "#80eb6e", "#05dda1", "#af7fae",
               "#349ad6", "#f2e5a9", "#43a405", "#f05017", "#f655f5",   
               "#01e034", "#3f7f56", "#016ecd", "#e2523f", "#b270e4")




#### PERFORMING PCA ###


dim(dataSE)    
dge <- DGEList(dataSE)

dge <- calcNormFactors(object=dge, method="TMM")
EM <- cpm(x=dge, log=TRUE)

# Perfrom PCA
pca <- prcomp(x=t(EM), scale=TRUE, center=TRUE)



### this code chunck ca be used use the % variance described by each PC

#getting variance proportion of PCs
#std<-pca$sdev #compute standard deviation of each principal component
#var<-std^2  #compute variance
#var[1:10]#check variance of first 10 components
#prop_varex <- var*100/sum(var)#proportion of variance explained
#prop_varex[1:10] # get varinces in percntages



######## PLOTTING PCA 2D sctter plots #######


library(extrafont)
loadfonts(device = "win")

#T/N  
dev.new()
ggplot(data=as.data.frame(pca$x),aes(x=PC1,y=PC2,col=samples.matrix$condition))+ 
  geom_point(size=2.1,alpha=0.9)+ 
  scale_colour_manual(values=tnPalette)+
  ggtitle("PCA plot of cancer and normal samples")+
  labs(x="PC1 (11.04%)", y="PC2 (8.63%)")+
  theme_classic(base_size = 12)+
  theme(#legend.position="bottom",
    legend.title=element_blank(),
    plot.title = element_text(hjust = 0.1),
    text=element_text(size=14,  family="serif"))

#PAM50   
dev.new()
ggplot(data=as.data.frame(pca$x),aes(x=PC1,y=PC2,col=samples.matrix$PAM50))+ 
  geom_point(size=2.1,alpha=0.9)+ 
  scale_colour_manual(values=cbPalette)+
  ggtitle("PCA plot of PAM50 subtypes and normal samples")+
  labs(x="PC1 (11.04%)", y="PC2 (8.63%)")+
  theme_classic(base_size = 12)+
  theme(#legend.position="bottom",
        legend.title=element_blank(),
        plot.title = element_text(hjust = 0.1),
        text=element_text(size=14,  family="serif"))

#morphology
dev.new()
ggplot(data=as.data.frame(pca$x),aes(x=PC3,y=PC2,col=samples.matrix$tumourTypes))+ 
  geom_point(size=2.1,alpha=1)+ 
  scale_colour_manual(values=tumourTypCol)+
  ggtitle("PCA plot of morphology groups and normal samples")+
  labs(x="PC3 (5.43%)", y="PC2 (8.63%)")+
  theme_classic(base_size = 10)+
  theme(legend.position="bottom",
        legend.title=element_blank(),
        plot.title = element_text(hjust = 0.6),
        text=element_text(size=14,  family="serif"))


#stage
ggplot(data=as.data.frame(pca$x),aes(x=PC1,y=PC2,col=samples.matrix$tumourStages))+ 
  geom_point(size=2.1,alpha=0.9)+ 
  scale_colour_manual(values=tumourStagesCol)+
  ggtitle("PCA plot of cancer stages and normal samples")+
  labs(x="PC1 (11.04%)", y="PC2 (8.63%)")+
  theme_classic(base_size = 10)+
  theme(legend.position="bottom",
        legend.title=element_blank(),
        plot.title = element_text(hjust = 0.7),
        text=element_text(size=14,  family="serif"))


######## PLOTTING 1D PLOTS ##########



#t/n
d <- data.frame(samples.matrix, pca$x)
p <- d %>%  tbl_df %>%  gather(key="PC", value="score", contains("PC")) %>% filter(PC %in% paste0("PC", 1:9))

ggplot(p, aes(x=condition, y=score, fill=condition)) + 
  geom_boxplot(outlier.colour=NaN) + 
  #geom_jitter(alpha=0.5) + 
  facet_wrap(~PC) + 
  scale_fill_manual(values=tnPalette)+
  ggtitle(" 1D PCA plot showing varince along PC1-PC9 \nfor cancer and normal samples")+
  labs(x=" ", y=" ")+
  theme(#legend.position="bottom",
        legend.title=element_blank(),
        plot.title = element_text(hjust = 0.5),
        text=element_text(size=14,  family="serif"))


#PAM50
d <- data.frame(samples.matrix, pca$x)
p <- d %>%  tbl_df %>%  gather(key="PC", value="score", contains("PC")) %>% filter(PC %in% paste0("PC", 1:9))

ggplot(p, aes(x=PAM50, y=score, fill=PAM50)) + 
  geom_boxplot(outlier.colour=NaN) + 
  facet_wrap(~PC) + 
  scale_fill_manual(values=cbPalette)+
  ggtitle("1D PCA plot showing varince along PC1-PC9 \nfor PAM50 subtypes and normal samples")+
  theme(axis.text.x=element_blank(),
      axis.ticks.x=element_blank(),
      legend.title=element_blank(),
      plot.title = element_text(hjust = 0.5),
      text=element_text(size=14,  family="serif"))




#tumour morphology
d <- data.frame(samples.matrix, pca$x)
p <- d %>%  tbl_df %>%  gather(key="PC", value="score", contains("PC")) %>% filter(PC %in% paste0("PC", 1:9))

ggplot(p, aes(x=tumourTypes, y=score, fill=tumourTypes)) + 
  geom_boxplot(outlier.colour=NaN) + 
  facet_wrap(~PC) + 
  scale_fill_manual(values=tumourTypCol)+
  ggtitle("1D PCA plot showing varince along PC1-PC9 \nfor morphology groups and normal samples")+
  labs(x="Morphology")+
  guides(fill=guide_legend(nrow=2))+
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        legend.title=element_blank(),
        legend.position="bottom",
        plot.title = element_text(hjust = 0.5),
        text=element_text(size=14,  family="serif"))




#tumour stages
d <- data.frame(samples.matrix, pca$x)
p <- d %>%  tbl_df %>%  gather(key="PC", value="score", contains("PC")) %>% filter(PC %in% paste0("PC", 1:9))

dev.new()
ggplot(p, aes(x=tumourStages, y=score, fill=tumourStages)) + 
  geom_boxplot(outlier.colour=NaN) + 
  facet_wrap(~PC) + 
  scale_fill_manual(values=tumourStagesCol)+
  ggtitle("1D PCA plot showing varince along PC1-PC9 \nfor cancer stages and normal samples")+
  labs(x="Cancer stages")+
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        legend.title=element_blank(),
        legend.position="bottom",
        plot.title = element_text(hjust = 0.5),
        text=element_text(size=14,  family="serif"))


#year taken
d <- data.frame(samples.matrix, pca$x)
p <- d %>%  tbl_df %>%  gather(key="PC", value="score", contains("PC")) %>% filter(PC %in% paste0("PC", 1:6))

dev.new()
ggplot(p, aes(x=year_diagnosed, y=score, fill=year_diagnosed)) + 
  geom_boxplot(outlier.colour=NaN) + 
  scale_fill_manual(values=year_cols)+
  geom_jitter(alpha=0.2, aes(color=year_diagnosed)) +
  scale_color_manual(values=year_cols)+
  facet_wrap(~PC) + 
  theme(legend.position="bottom",
        axis.text.x = element_text(angle = 90, hjust = 1, size=7))+
  guides(fill=guide_legend(nrow=3))+
  ggtitle("1D PCA plot showing varince along PC1-PC6 for years samples were taken")+
  labs(x="Year sample taken")+
  theme(legend.title=element_blank(),
        plot.title = element_text(hjust = 0.5),
        text=element_text(size=14,  family="serif"))


# source site
names(samples.matrix)[19]<-'source_site'
d <- data.frame(samples.matrix, pca$x)
p <- d %>%  tbl_df %>%  gather(key="PC", value="score", contains("PC")) %>% filter(PC %in% paste0("PC", 1:6))
dev.new()
ggplot(p, aes(x=source_site, y=score, fill=source_site)) + 
  geom_boxplot(outlier.colour=NaN) + 
  scale_fill_manual(values=source_cols)+
  geom_jitter(alpha=0.1, aes(color=source_site)) +
  scale_color_manual(values=source_cols)+
  facet_wrap(~PC) + 
  theme(legend.position="bottom",
        axis.text.x = element_text(angle =0, hjust = 0.5, size=5))+
  guides(fill=guide_legend(nrow=2))+
  ggtitle("1D PCA plot showing varince along PC1-PC6 for sample source site ")+ # (only cancer samples)
  labs(x="Sample source site codes")+
  theme(legend.title=element_blank(),
        plot.title = element_text(hjust = 0.5),
        text=element_text(size=14,  family="serif"))


library(scales)
show_col(source_cols)

#age groups!

ggplot(p, aes(x=ageGroups, y=score, fill=ageGroups)) + 
  geom_boxplot(outlier.colour=NaN) + 
  scale_fill_manual(values=source_cols)+
  scale_color_manual(values=source_cols)+
  facet_wrap(~PC) + 
  theme(legend.position="bottom",
        axis.text.x = element_text(angle =45, hjust = 0.5, size=8))+
  guides(fill=guide_legend(nrow=1))+
  ggtitle("1D PCA plot showing varince along PC1-PC6 \nfor age groups")+
  labs(x="Patients age groups")+
  theme(legend.title=element_blank(),
        plot.title = element_text(hjust = 0.5))





