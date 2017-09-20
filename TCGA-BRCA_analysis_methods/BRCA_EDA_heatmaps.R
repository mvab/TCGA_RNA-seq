#! /usr/bin/env Rscript
rm(list=ls(all=TRUE))

#library(TCGAbiolinks)
library(SummarizedExperiment)
library(plyr)
library(edgeR)

library(ggplot2)
library(RColorBrewer)

library(limma)

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

normal <- PAMNPout$normal_barcodes #get normal barcodes for later


# exclude mrphology samples Ductal mixed with Others -- it was added at the beginning, 
# but later it was decided not use it in the analysis. Therefore this step removes those samples 
# from dataSE and samples.matrix

removedsamples<-removeDuctalOther(dataSE, samples.matrix)
dataSE<-removedsamples$dataSE
samples.matrix<-removedsamples$samples.matrix
dim(dataSE)
dim(samples.matrix)





###### SETUP for COLOURBARS of PATIENT CLASSES #####

#get pationt classes info
names(samples.matrix)
design<-subset(samples.matrix[,c('tumourTypes','tumourStages','PAM50')])
row.names(design)<-samples.matrix$barcode
colnames(design)<-c("Morphology", "Stages", "PAM50")

#reorder annotations
design$PAM50 = factor(design$PAM50, levels = c("Basal-like", "HER2-enriched", "Luminal A", "Luminal B", "Normal-like", "Normal"))
design$Stages = factor(design$Stages, levels = c("stage1", "stage2", "stage3", "stage4", "unknown")) 
design$Morphology = factor(design$Morphology,  
              levels = c("Ductal and lobular mixed","Ductal carcinoma","Lobular carcinoma", 
                         "Metaplastic carcinoma" ,"Mucinous carcinoma", "Normal",  "Other")) #other is a mix of few samples of different ones

# colours for each subgroup
PAM50Col <- c( "#E69F00", "#0072B2", "#F0E442","red3", "#CC79A7","#009E73")##D55E00->altred
tumourStgCol = c( "red","green","blue", "yellow2", "white")
tumourTypCol = c( "#F0E442", "#56B4E9", "#CC79A7", "#D55E00", "#0072B2", "#009E73", "#E69F00")

#assign colours to subgroups
names(PAM50Col)<-levels(design$PAM50)
names(tumourTypCol)<-levels(design$Morphology)
names(tumourStgCol)<-levels(design$Stages)

# create list
annColour <-list(
  PAM50=PAM50Col,
  Morphology=tumourTypCol,
  Stages=tumourStgCol
)
######



##### SETUP FOR HEATMAP ########


# normalise with log
dge <- DGEList(dataSE)
dge <- calcNormFactors(object=dge, method="TMM")
EM <- cpm(x=dge, log=TRUE)


## as there are 17k genes, we are only going to look at the ones that are have the most change in expression
## selecting top 1000 genes:

# calculate variance for each gene (log2 data), take top 1000 highest ->
#look at the stuff that is actually changing
#use only cancer samples for varince calculation


EMc <- EM[, !(colnames(EM) %in% normal)] #select only cancer samples
EM_var<- apply(EMc, MARGIN=1, var) #calc var across the row
length(EM_var)
top_genes<- sort(EM_var, decreasing = TRUE)[1:1000] # sort data with hightest variance genes at the top, select top 1000

# cancer and normal
EM2<- subset(EM[c(names(top_genes)),]) #select top 1000 genes in the original EM object (that has both T and N samples)
dim(EM2) # should have 1000 rows



###### PLOTTING the HEATMAP


pheatmap::pheatmap(mat = as.matrix(EM2), color = colorRampPalette(brewer.pal(9,"YlGnBu"))(20),
                   clustering_distance_rows = 'euclidean', 
                   clustering_distance_cols = 'euclidean',
                   annotation_col=design,
                   annotation_colors = annColour,
                   cluster_cols = T, cluster_rows = T, 
                   show_rownames = F,show_colnames = F,
                   fontsize = 5)


