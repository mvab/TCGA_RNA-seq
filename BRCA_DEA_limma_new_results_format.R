#! /usr/bin/env Rscript
rm(list=ls(all=TRUE))

#library(TCGAbiolinks)
library(SummarizedExperiment)
library(plyr)
library(edgeR)
library(limma)
library(ggplot2)
library(RColorBrewer)
library(dplyr)


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

# only PAM50 samples 
PAMNPout<-addXtraPAMandNormal(samples.matrix)# (807 +104 -39) + 112 = 984
samples.matrix<-PAMNPout$samples.matrix 
#save(samples.matrix, file="samplesMatrix_full.rda")
dim(samples.matrix)
sampleTokeepSE <- PAMNPout$samplesToKeep
dataSE<- dataSE[,colnames(dataSE) %in% sampleTokeepSE]
dim(dataSE) 

## 

#get gene lenghts for genes in the analysis : NB this 
#getGeneLenght_out<-getGeneLenght(dataSE)
#gene_lengths<-getGeneLenght_out$gene_lengths
#dataSE<-getGeneLenght_out$dataSE
#dim(dataSE)
#dim(gene_lengths)

# edgeR object
dge <- DGEList(counts=dataSE)#, genes = data.frame(Length= as.numeric(gene_lengths$gene_length)) ) #here will be adding genes: genes = Ann)
head(dge$genes)

#### IF want to RPKM
#yr <- calcNormFactors(dge)
# rpkm will use the normalized effective library sizes to compute rpkm instead of the raw library sizes. 
#RPKM <- rpkm(yr)
#dim(RPKM) #17368   984 matrix
#keep.exprs <- rowSums(RPKM > 1) >=50     
#dge <- dge[keep.exprs,, keep.lib.sizes=FALSE]
#dim(dge) # 17098 


# filtering low expression by cpm

#how many genes have at least 50 samples with 0 expression
#Keep genes with total counts more than 50.
cpm <- cpm(dge, log=FALSE)
#keep.exprs <- rowSums(cpm > 1) >= 50        #A CPM value of 1 is equivalent to a log-CPM value of 0.
keep.exprs <- rowSums(cpm > 2) >= 19
dge <- dge[keep.exprs,, keep.lib.sizes=FALSE]
dim(dge) # 

#check autophagy in the current geneset
#autophagy_genes<- as.vector(read.table("autopahagy_genes.txt", as.is = T, header = FALSE))$V1
#shared <- intersect(autophagy_genes,rownames(dge))
#print(paste0("Total number of genes in analysis: ", length(rownames(dge))))
#print(paste0("Autophagy genes: ", length(shared)))




addSampleData<-function(y, samples.matrix){
  
  # adding samples information
  y$samples$condition <- as.factor(samples.matrix$condition)
  y$samples$PAM50 <- as.factor(samples.matrix$PAM50)
  y$samples$morphology <- as.factor(samples.matrix$tumourTypes) # from sample lists!
  y$samples$stages <- as.factor(samples.matrix$tumourStages) # from sample lists!
  y$samples$year <- as.factor(samples.matrix$year_diagnosed)
  y$samples$tss <- as.factor(samples.matrix$tss)
  y$samples$age <- as.factor(samples.matrix$ageGroups)
  
  #stage + PAM50
  y$samples$Group1 <- factor(paste( gsub(" ", "", samples.matrix$PAM50) , samples.matrix$tumourStages,sep="."))
  #morphology +PAM50
  y$samples$Group2 <- factor(paste( gsub(" ", "", samples.matrix$PAM50) , samples.matrix$tumourTypes,sep="."))
  
  
  #fixing NAs
  y$samples$age <- as.character(y$samples$age)
  y$samples[is.na(y$samples$age),]$age<-"Unknown"
  y$samples[y$samples$age=="70+",]$age<-"70andUp"
  y$samples$age <- as.factor(y$samples$age)
  
  y$samples$year <- as.character(y$samples$year)
  y$samples[is.na(y$samples$year),]$year<-"Unknown"
  y$samples$year <- as.factor(y$samples$year)
  
  
  ## making normal the baselayer
  y$samples$condition = relevel(y$samples$condition, ref="normal")
  y$samples$PAM50 = relevel(y$samples$PAM50, ref="Normal")
  y$samples$morphology = relevel(y$samples$morphology, ref="Normal")
  y$samples$stages = relevel(y$samples$stages, ref="Normal")
  y$samples$Group1 = relevel(y$samples$Group1, ref = "Normal.Normal")
  y$samples$Group1 = relevel(y$samples$Group2, ref = "Normal.Normal")
  
  
  return (y)
}
dge<-addSampleData(dge,samples.matrix)
# Scale normalisation: correct for library size
#to match the between-sample distributions of gene counts
#in terms of parameters such as quantiles
y <- calcNormFactors(dge)



#trying MDS condition 
#condition <- substring(y$samples$condition,1,6) 

#plotMDS(y, labels=condition, top=1000, col=ifelse(condition=="cancer","blue","red"),        gene.selection="common", prior.count = 5)

#emergency PCA
#pca <- prcomp(x=t(log.cpm), scale=TRUE, center=TRUE)
#qplot(data=as.data.frame(pca$x), x=PC1, y=PC2, geom=c("point"), color=y$samples$morphology)

#trying MDS PAM50
#pam50 <- substring(samples.matrix$PAM50,1,9) 
#plotMDS(y, labels=pam50, top=50, gene.selection="common", prior.count = 5, pch = 'o')
#plotMDS(y, labels=pam50, col=col.concent[y$samples$PAM50], prior.count=5)





#design <- model.matrix(~0+PAM50, data=y$samples) #setting up model contrasts is more straight forward in the absence of an intercept for
#design <- model.matrix(~tss+PAM50  +morphology + stages + age +year , data=y$samples) #setting up model contrasts is more straight forward in the absence of an intercept for

design <- model.matrix(~0 + PAM50 + age + stages  + morphology + tss + year , data=y$samples) #nested interaction # makes all possible pairings (alt to manual contarsts)
design <- model.matrix(~0 +  Group1 + PAM50 +   age + stages  + morphology + tss + year , data=y$samples) #nested interaction # makes all possible pairings (alt to manual contarsts)

#design <- model.matrix(~0+Group1, data=y$samples) #setting up model contrasts is more straight forward in the absence of an intercept for

colnames(design) <- gsub("Group1", "", colnames(design))
colnames(design) <- gsub("PAM50", "", colnames(design))
colnames(design) <- gsub("morphology", "", colnames(design))
colnames(design) <- gsub("stages", "", colnames(design))
colnames(design) <- gsub(":age", ":", colnames(design))


colnames(design) <- gsub(" ", "", colnames(design))
colnames(design) <- gsub("-l", "L", colnames(design))
colnames(design) <- gsub("-", "", colnames(design))

#colnames(design) <- gsub(":t", "t", colnames(design))

colnames(design) 

# contr matrices
#######       


contr.matrix <- makeContrasts(
  ######## BY PAM50  
  LuminalA.stage1vsNorm = LuminalA.stage1 - Normal.Normal,  
  LuminalA.stage2vsNorm = LuminalA.stage2 - Normal.Normal,
  LuminalA.stage3vsNorm = LuminalA.stage3 - Normal.Normal,
  LuminalA.stage4vsNorm = LuminalA.stage4- Normal.Normal,
  
  LuminalA.stage1vsLuminalA.stage2 = LuminalA.stage1 - LuminalA.stage2,
  LuminalA.stage1vsLuminalA.stage3 = LuminalA.stage1 - LuminalA.stage3,
  LuminalA.stage1vsLuminalA.stage4 = LuminalA.stage1 - LuminalA.stage4,
  
  LuminalA.stage2vsLuminalA.stage1 = LuminalA.stage2 - LuminalA.stage1,
  LuminalA.stage2vsLuminalA.stage3 = LuminalA.stage2 - LuminalA.stage3,
  LuminalA.stage2vsLuminalA.stage4 = LuminalA.stage2 - LuminalA.stage4,
  
  LuminalA.stage3vsLuminalA.stage4 = LuminalA.stage3 - LuminalA.stage4,
  LuminalA.stage3vsLuminalA.stage2 = LuminalA.stage3 - LuminalA.stage2,
  LuminalA.stage3vsLuminalA.stage1 = LuminalA.stage3 - LuminalA.stage1,
  
  
  
  LuminalB.stage1vsNorm = LuminalB.stage1 - Normal.Normal,  
  LuminalB.stage2vsNorm = LuminalB.stage2 - Normal.Normal,
  LuminalB.stage3vsNorm = LuminalB.stage3 - Normal.Normal,
  LuminalB.stage4vsNorm = LuminalB.stage4- Normal.Normal,
  
  LuminalB.stage1vsLuminalB.stage2 = LuminalB.stage1 - LuminalB.stage2,
  LuminalB.stage1vsLuminalB.stage3 = LuminalB.stage1 - LuminalB.stage3,
  LuminalB.stage1vsLuminalB.stage4 = LuminalB.stage1 - LuminalB.stage4,
  
  LuminalB.stage2vsLuminalB.stage1 = LuminalB.stage2 - LuminalB.stage1,
  LuminalB.stage2vsLuminalB.stage3 = LuminalB.stage2 - LuminalB.stage3,
  LuminalB.stage2vsLuminalB.stage4 = LuminalB.stage2 - LuminalB.stage4,
  
  LuminalB.stage3vsLuminalB.stage4 = LuminalB.stage3 - LuminalB.stage4,
  LuminalB.stage3vsLuminalB.stage2 = LuminalB.stage3 - LuminalB.stage2,
  LuminalB.stage3vsLuminalB.stage1 = LuminalB.stage3 - LuminalB.stage1,
  
  
  
  BasalLike.stage1vsNorm = BasalLike.stage1 - Normal.Normal,  
  BasalLike.stage2vsNorm = BasalLike.stage2 - Normal.Normal,
  BasalLike.stage3vsNorm = BasalLike.stage3 - Normal.Normal,
  BasalLike.stage4vsNorm = BasalLike.stage4- Normal.Normal,
  
  BasalLike.stage1vsBasalLike.stage2 = BasalLike.stage1 - BasalLike.stage2,
  BasalLike.stage1vsBasalLike.stage3 = BasalLike.stage1 - BasalLike.stage3,
  BasalLike.stage1vsBasalLike.stage4 = BasalLike.stage1 - BasalLike.stage4,
  
  BasalLike.stage2vsBasalLike.stage1 = BasalLike.stage2 - BasalLike.stage1,
  BasalLike.stage2vsBasalLike.stage3 = BasalLike.stage2 - BasalLike.stage3,
  BasalLike.stage2vsBasalLike.stage4 = BasalLike.stage2 - BasalLike.stage4,
  
  BasalLike.stage3vsBasalLike.stage4 = BasalLike.stage3 - BasalLike.stage4,
  BasalLike.stage3vsBasalLike.stage2 = BasalLike.stage3 - BasalLike.stage2,
  BasalLike.stage3vsBasalLike.stage1 = BasalLike.stage3 - BasalLike.stage1,
  
  
  
  HER2enriched.stage1vsNorm = HER2enriched.stage1 - Normal.Normal,  
  HER2enriched.stage2vsNorm = HER2enriched.stage2 - Normal.Normal,
  HER2enriched.stage3vsNorm = HER2enriched.stage3 - Normal.Normal,
  HER2enriched.stage4vsNorm = HER2enriched.stage4- Normal.Normal,
  
  HER2enriched.stage1vsHER2enriched.stage2 = HER2enriched.stage1 - HER2enriched.stage2,
  HER2enriched.stage1vsHER2enriched.stage3 = HER2enriched.stage1 - HER2enriched.stage3,
  HER2enriched.stage1vsHER2enriched.stage4 = HER2enriched.stage1 - HER2enriched.stage4,
  
  HER2enriched.stage2vsHER2enriched.stage1 = HER2enriched.stage2 - HER2enriched.stage1,
  HER2enriched.stage2vsHER2enriched.stage3 = HER2enriched.stage2 - HER2enriched.stage3,
  HER2enriched.stage2vsHER2enriched.stage4 = HER2enriched.stage2 - HER2enriched.stage4,
  
  HER2enriched.stage3vsHER2enriched.stage4 = HER2enriched.stage3 - HER2enriched.stage4,
  HER2enriched.stage3vsHER2enriched.stage2 = HER2enriched.stage3 - HER2enriched.stage2,
  HER2enriched.stage3vsHER2enriched.stage1 = HER2enriched.stage3 - HER2enriched.stage1,
  
  
  NormalLike.stage1vsNorm = NormalLike.stage1 - Normal.Normal,  
  NormalLike.stage2vsNorm = NormalLike.stage2 - Normal.Normal,
  NormalLike.stage3vsNorm = NormalLike.stage3 - Normal.Normal,
  
  NormalLike.stage1vsNormalLike.stage2 = NormalLike.stage1 - NormalLike.stage2,
  NormalLike.stage1vsNormalLike.stage3 = NormalLike.stage1 - NormalLike.stage3,
  
  NormalLike.stage2vsNormalLike.stage1 = NormalLike.stage2 - NormalLike.stage1,
  NormalLike.stage2vsNormalLike.stage3 = NormalLike.stage2 - NormalLike.stage3,
  
  NormalLike.stage3vsNormalLike.stage2 = NormalLike.stage3 - NormalLike.stage2,
  NormalLike.stage3vsNormalLike.stage1 = NormalLike.stage3 - NormalLike.stage1,
  
  # normallike.stage4 does not exist!
  
  #######BY STAGE
  
  stage1.LumAvsLumB = LuminalA.stage1 - LuminalB.stage1,     
  stage1.LumAvsBasal = LuminalA.stage1 - BasalLike.stage1,
  stage1.LumAvsHER2 = LuminalA.stage1 - HER2enriched.stage1,
  stage1.LumAvsNormLike = LuminalA.stage1 - NormalLike.stage1,
  stage1.LumAvsNormal = LuminalA.stage1 - Normal.Normal,
  
  stage1.LumBvsLumA = LuminalB.stage1 - LuminalA.stage1,
  stage1.LumBvsBasal = LuminalB.stage1 - BasalLike.stage1,
  stage1.LumBvsHER2 = LuminalB.stage1 - HER2enriched.stage1,
  stage1.LumBvsNormLike = LuminalB.stage1 - NormalLike.stage1,
  stage1.LumBvsNormal = LuminalB.stage1 - Normal.Normal,
  
  stage1.BasalvsLumA = BasalLike.stage1- LuminalA.stage1,
  stage1.BasalvsLumB = BasalLike.stage1- LuminalB.stage1,
  stage1.BasalvsHER2 = BasalLike.stage1- HER2enriched.stage1,
  stage1.BasalvsNormLike = BasalLike.stage1 - NormalLike.stage1,
  stage1.BasalvsNormal = BasalLike.stage1 - Normal.Normal,
  
  stage1.HER2vsLumA = HER2enriched.stage1- LuminalA.stage1,
  stage1.HER2vsLumB = HER2enriched.stage1- LuminalB.stage1,
  stage1.HER2vsBasal = HER2enriched.stage1- BasalLike.stage1,
  stage1.HER2vsNormLike = HER2enriched.stage1 - NormalLike.stage1,
  stage1.HER2vsNormal = HER2enriched.stage1 - Normal.Normal,
  
  stage1.NormalLikevsLumA = NormalLike.stage1- LuminalA.stage1,
  stage1.NormalLikevsLumB = NormalLike.stage1- LuminalB.stage1,
  stage1.NormalLikevsBasal = NormalLike.stage1- BasalLike.stage1,
  stage1.NormalLikevsHER2 = NormalLike.stage1- HER2enriched.stage1,
  stage1.NormalLikevsNormal =NormalLike.stage1 - Normal.Normal,     
  
  
  
  
  stage2.LumAvsLumB = LuminalA.stage2 - LuminalB.stage2,     
  stage2.LumAvsBasal = LuminalA.stage2 - BasalLike.stage2,
  stage2.LumAvsHER2 = LuminalA.stage2 - HER2enriched.stage2,
  stage2.LumAvsNormLike = LuminalA.stage2 - NormalLike.stage2,
  stage2.LumAvsNormal = LuminalA.stage2 - Normal.Normal,
  
  stage2.LumBvsLumA = LuminalB.stage2 - LuminalA.stage2,
  stage2.LumBvsBasal = LuminalB.stage2 - BasalLike.stage2,
  stage2.LumBvsHER2 = LuminalB.stage2 - HER2enriched.stage2,
  stage2.LumBvsNormLike = LuminalB.stage2 - NormalLike.stage2,
  stage2.LumBvsNormal = LuminalB.stage2 - Normal.Normal,
  
  stage2.BasalvsLumA = BasalLike.stage2- LuminalA.stage2,
  stage2.BasalvsLumB = BasalLike.stage2- LuminalB.stage2,
  stage2.BasalvsHER2 = BasalLike.stage2- HER2enriched.stage2,
  stage2.BasalvsNormLike = BasalLike.stage2 - NormalLike.stage2,
  stage2.BasalvsNormal = BasalLike.stage2 - Normal.Normal,
  
  stage2.HER2vsLumA = HER2enriched.stage2- LuminalA.stage2,
  stage2.HER2vsLumB = HER2enriched.stage2- LuminalB.stage2,
  stage2.HER2vsBasal = HER2enriched.stage2- BasalLike.stage2,
  stage2.HER2vsNormLike = HER2enriched.stage2 - NormalLike.stage2,
  stage2.HER2vsNormal = HER2enriched.stage2 - Normal.Normal,
  
  stage2.NormalLikevsLumA = NormalLike.stage2- LuminalA.stage2,
  stage2.NormalLikevsLumB = NormalLike.stage2- LuminalB.stage2,
  stage2.NormalLikevsBasal = NormalLike.stage2- BasalLike.stage2,
  stage2.NormalLikevsHER2 = NormalLike.stage2- HER2enriched.stage2,
  stage2.NormalLikevsNormal =NormalLike.stage2 - Normal.Normal,  
  
  
  
  stage3.LumAvsLumB = LuminalA.stage3 - LuminalB.stage3,     
  stage3.LumAvsBasal = LuminalA.stage3 - BasalLike.stage3,
  stage3.LumAvsHER2 = LuminalA.stage3 - HER2enriched.stage3,
  stage3.LumAvsNormLike = LuminalA.stage3 - NormalLike.stage3,
  stage3.LumAvsNormal = LuminalA.stage3 - Normal.Normal,
  
  stage3.LumBvsLumA = LuminalB.stage3 - LuminalA.stage3,
  stage3.LumBvsBasal = LuminalB.stage3 - BasalLike.stage3,
  stage3.LumBvsHER2 = LuminalB.stage3 - HER2enriched.stage3,
  stage3.LumBvsNormLike = LuminalB.stage3 - NormalLike.stage3,
  stage3.LumBvsNormal = LuminalB.stage3 - Normal.Normal,
  
  stage3.BasalvsLumA = BasalLike.stage3- LuminalA.stage3,
  stage3.BasalvsLumB = BasalLike.stage3- LuminalB.stage3,
  stage3.BasalvsHER2 = BasalLike.stage3- HER2enriched.stage3,
  stage3.BasalvsNormLike = BasalLike.stage3 - NormalLike.stage3,
  stage3.BasalvsNormal = BasalLike.stage3 - Normal.Normal,
  
  stage3.HER2vsLumA = HER2enriched.stage3- LuminalA.stage3,
  stage3.HER2vsLumB = HER2enriched.stage3- LuminalB.stage3,
  stage3.HER2vsBasal = HER2enriched.stage3- BasalLike.stage3,
  stage3.HER2vsNormLike = HER2enriched.stage3 - NormalLike.stage3,
  stage3.HER2vsNormal = HER2enriched.stage3 - Normal.Normal,
  
  stage3.NormalLikevsLumA = NormalLike.stage3- LuminalA.stage3,
  stage3.NormalLikevsLumB = NormalLike.stage3- LuminalB.stage3,
  stage3.NormalLikevsBasal = NormalLike.stage3- BasalLike.stage3,
  stage3.NormalLikevsHER2 = NormalLike.stage3- HER2enriched.stage3,
  stage3.NormalLikevsNormal =NormalLike.stage3 - Normal.Normal,         
  
  
  stage4.LumAvsLumB = LuminalA.stage4 - LuminalB.stage4,     
  stage4.LumAvsBasal = LuminalA.stage4 - BasalLike.stage4,
  stage4.LumAvsHER2 = LuminalA.stage4 - HER2enriched.stage4,
  stage4.LumAvsNormal = LuminalA.stage4 - Normal.Normal,
  
  stage4.LumBvsLumA = LuminalB.stage4 - LuminalA.stage4,
  stage4.LumBvsBasal = LuminalB.stage4 - BasalLike.stage4,
  stage4.LumBvsHER2 = LuminalB.stage4 - HER2enriched.stage4,
  stage4.LumBvsNormal = LuminalB.stage4 - Normal.Normal,
  
  stage4.BasalvsLumA = BasalLike.stage4- LuminalA.stage4,
  stage4.BasalvsLumB = BasalLike.stage4- LuminalB.stage4,
  stage4.BasalvsHER2 = BasalLike.stage4- HER2enriched.stage4,
  stage4.BasalvsNormal = BasalLike.stage4 - Normal.Normal,
  
  stage4.HER2vsLumA = HER2enriched.stage4- LuminalA.stage4,
  stage4.HER2vsLumB = HER2enriched.stage4- LuminalB.stage4,
  stage4.HER2vsBasal = HER2enriched.stage4- BasalLike.stage4,
  stage4.HER2vsNormal = HER2enriched.stage4 - Normal.Normal,
  
  levels = colnames(design)) 


#group1
contr.matrix <- makeContrasts( LuminalA.stage1vsNorm = LuminalA.stage1 - Normal.Normal,  
                               LuminalA.stage2vsNorm = LuminalA.stage2 - Normal.Normal,
                               LuminalA.stage3vsNorm = LuminalA.stage3 - Normal.Normal,
                               LuminalA.stage4vsNorm = LuminalA.stage4- Normal.Normal,
                               
                               LuminalA.stage1vsLuminalA.stage2 = LuminalA.stage1 - LuminalA.stage2,
                               LuminalA.stage1vsLuminalA.stage3 = LuminalA.stage1 - LuminalA.stage3,
                               LuminalA.stage1vsLuminalA.stage4 = LuminalA.stage1 - LuminalA.stage4,
                               
                               LuminalA.stage2vsLuminalA.stage1 = LuminalA.stage2 - LuminalA.stage1,
                               LuminalA.stage2vsLuminalA.stage3 = LuminalA.stage2 - LuminalA.stage3,
                               LuminalA.stage2vsLuminalA.stage4 = LuminalA.stage2 - LuminalA.stage4,
                               
                               LuminalA.stage3vsLuminalA.stage4 = LuminalA.stage3 - LuminalA.stage4,
                               LuminalA.stage3vsLuminalA.stage2 = LuminalA.stage3 - LuminalA.stage2,
                               LuminalA.stage3vsLuminalA.stage1 = LuminalA.stage3 - LuminalA.stage1,
                               
                               
                               
                               LuminalB.stage1vsNorm = LuminalB.stage1 - Normal.Normal,  
                               LuminalB.stage2vsNorm = LuminalB.stage2 - Normal.Normal,
                               LuminalB.stage3vsNorm = LuminalB.stage3 - Normal.Normal,
                               LuminalB.stage4vsNorm = LuminalB.stage4- Normal.Normal,
                               
                               LuminalB.stage1vsLuminalB.stage2 = LuminalB.stage1 - LuminalB.stage2,
                               LuminalB.stage1vsLuminalB.stage3 = LuminalB.stage1 - LuminalB.stage3,
                               LuminalB.stage1vsLuminalB.stage4 = LuminalB.stage1 - LuminalB.stage4,
                               
                               LuminalB.stage2vsLuminalB.stage1 = LuminalB.stage2 - LuminalB.stage1,
                               LuminalB.stage2vsLuminalB.stage3 = LuminalB.stage2 - LuminalB.stage3,
                               LuminalB.stage2vsLuminalB.stage4 = LuminalB.stage2 - LuminalB.stage4,
                               
                               LuminalB.stage3vsLuminalB.stage4 = LuminalB.stage3 - LuminalB.stage4,
                               LuminalB.stage3vsLuminalB.stage2 = LuminalB.stage3 - LuminalB.stage2,
                               LuminalB.stage3vsLuminalB.stage1 = LuminalB.stage3 - LuminalB.stage1,
                               
                               
                               
                               BasalLike.stage1vsNorm = BasalLike.stage1 - Normal.Normal,  
                               BasalLike.stage2vsNorm = BasalLike.stage2 - Normal.Normal,
                               BasalLike.stage3vsNorm = BasalLike.stage3 - Normal.Normal,
                               BasalLike.stage4vsNorm = BasalLike.stage4- Normal.Normal,
                               
                               BasalLike.stage1vsBasalLike.stage2 = BasalLike.stage1 - BasalLike.stage2,
                               BasalLike.stage1vsBasalLike.stage3 = BasalLike.stage1 - BasalLike.stage3,
                               BasalLike.stage1vsBasalLike.stage4 = BasalLike.stage1 - BasalLike.stage4,
                               
                               BasalLike.stage2vsBasalLike.stage1 = BasalLike.stage2 - BasalLike.stage1,
                               BasalLike.stage2vsBasalLike.stage3 = BasalLike.stage2 - BasalLike.stage3,
                               BasalLike.stage2vsBasalLike.stage4 = BasalLike.stage2 - BasalLike.stage4,
                               
                               BasalLike.stage3vsBasalLike.stage4 = BasalLike.stage3 - BasalLike.stage4,
                               BasalLike.stage3vsBasalLike.stage2 = BasalLike.stage3 - BasalLike.stage2,
                               BasalLike.stage3vsBasalLike.stage1 = BasalLike.stage3 - BasalLike.stage1,
                               
                               
                               
                               HER2enriched.stage1vsNorm = HER2enriched.stage1 - Normal.Normal,  
                               HER2enriched.stage2vsNorm = HER2enriched.stage2 - Normal.Normal,
                               HER2enriched.stage3vsNorm = HER2enriched.stage3 - Normal.Normal,
                               HER2enriched.stage4vsNorm = HER2enriched.stage4- Normal.Normal,
                               
                               HER2enriched.stage1vsHER2enriched.stage2 = HER2enriched.stage1 - HER2enriched.stage2,
                               HER2enriched.stage1vsHER2enriched.stage3 = HER2enriched.stage1 - HER2enriched.stage3,
                               HER2enriched.stage1vsHER2enriched.stage4 = HER2enriched.stage1 - HER2enriched.stage4,
                               
                               HER2enriched.stage2vsHER2enriched.stage1 = HER2enriched.stage2 - HER2enriched.stage1,
                               HER2enriched.stage2vsHER2enriched.stage3 = HER2enriched.stage2 - HER2enriched.stage3,
                               HER2enriched.stage2vsHER2enriched.stage4 = HER2enriched.stage2 - HER2enriched.stage4,
                               
                               HER2enriched.stage3vsHER2enriched.stage4 = HER2enriched.stage3 - HER2enriched.stage4,
                               HER2enriched.stage3vsHER2enriched.stage2 = HER2enriched.stage3 - HER2enriched.stage2,
                               HER2enriched.stage3vsHER2enriched.stage1 = HER2enriched.stage3 - HER2enriched.stage1,
                               
                               
                               NormalLike.stage1vsNorm = NormalLike.stage1 - Normal.Normal,  
                               NormalLike.stage2vsNorm = NormalLike.stage2 - Normal.Normal,
                               NormalLike.stage3vsNorm = NormalLike.stage3 - Normal.Normal,

                               NormalLike.stage1vsNormalLike.stage2 = NormalLike.stage1 - NormalLike.stage2,
                               NormalLike.stage1vsNormalLike.stage3 = NormalLike.stage1 - NormalLike.stage3,

                               NormalLike.stage2vsNormalLike.stage1 = NormalLike.stage2 - NormalLike.stage1,
                               NormalLike.stage2vsNormalLike.stage3 = NormalLike.stage2 - NormalLike.stage3,

                               NormalLike.stage3vsNormalLike.stage2 = NormalLike.stage3 - NormalLike.stage2,
                               NormalLike.stage3vsNormalLike.stage1 = NormalLike.stage3 - NormalLike.stage1,
                               
                               # normallike.stage4 does not exist!
                               
                               levels = colnames(design)) 

#group1 other 
contr.matrix <- makeContrasts( stage1.LumAvsLumB = LuminalA.stage1 - LuminalB.stage1,     
                               stage1.LumAvsBasal = LuminalA.stage1 - BasalLike.stage1,
                               stage1.LumAvsHER2 = LuminalA.stage1 - HER2enriched.stage1,
                               stage1.LumAvsNormLike = LuminalA.stage1 - NormalLike.stage1,
                               stage1.LumAvsNormal = LuminalA.stage1 - Normal.Normal,
                               
                               stage1.LumBvsLumA = LuminalB.stage1 - LuminalA.stage1,
                               stage1.LumBvsBasal = LuminalB.stage1 - BasalLike.stage1,
                               stage1.LumBvsHER2 = LuminalB.stage1 - HER2enriched.stage1,
                               stage1.LumBvsNormLike = LuminalB.stage1 - NormalLike.stage1,
                               stage1.LumBvsNormal = LuminalB.stage1 - Normal.Normal,
                               
                               stage1.BasalvsLumA = BasalLike.stage1- LuminalA.stage1,
                               stage1.BasalvsLumB = BasalLike.stage1- LuminalB.stage1,
                               stage1.BasalvsHER2 = BasalLike.stage1- HER2enriched.stage1,
                               stage1.BasalvsNormLike = BasalLike.stage1 - NormalLike.stage1,
                               stage1.BasalvsNormal = BasalLike.stage1 - Normal.Normal,
                               
                               stage1.HER2vsLumA = HER2enriched.stage1- LuminalA.stage1,
                               stage1.HER2vsLumB = HER2enriched.stage1- LuminalB.stage1,
                               stage1.HER2vsBasal = HER2enriched.stage1- BasalLike.stage1,
                               stage1.HER2vsNormLike = HER2enriched.stage1 - NormalLike.stage1,
                               stage1.HER2vsNormal = HER2enriched.stage1 - Normal.Normal,
                               
                               stage1.NormalLikevsLumA = NormalLike.stage1- LuminalA.stage1,
                               stage1.NormalLikevsLumB = NormalLike.stage1- LuminalB.stage1,
                               stage1.NormalLikevsBasal = NormalLike.stage1- BasalLike.stage1,
                               stage1.NormalLikevsHER2 = NormalLike.stage1- HER2enriched.stage1,
                               stage1.NormalLikevsNormal =NormalLike.stage1 - Normal.Normal,     
                               
                               
                               
                               
                               stage2.LumAvsLumB = LuminalA.stage2 - LuminalB.stage2,     
                               stage2.LumAvsBasal = LuminalA.stage2 - BasalLike.stage2,
                               stage2.LumAvsHER2 = LuminalA.stage2 - HER2enriched.stage2,
                               stage2.LumAvsNormLike = LuminalA.stage2 - NormalLike.stage2,
                               stage2.LumAvsNormal = LuminalA.stage2 - Normal.Normal,
                               
                               stage2.LumBvsLumA = LuminalB.stage2 - LuminalA.stage2,
                               stage2.LumBvsBasal = LuminalB.stage2 - BasalLike.stage2,
                               stage2.LumBvsHER2 = LuminalB.stage2 - HER2enriched.stage2,
                               stage2.LumBvsNormLike = LuminalB.stage2 - NormalLike.stage2,
                               stage2.LumBvsNormal = LuminalB.stage2 - Normal.Normal,
                               
                               stage2.BasalvsLumA = BasalLike.stage2- LuminalA.stage2,
                               stage2.BasalvsLumB = BasalLike.stage2- LuminalB.stage2,
                               stage2.BasalvsHER2 = BasalLike.stage2- HER2enriched.stage2,
                               stage2.BasalvsNormLike = BasalLike.stage2 - NormalLike.stage2,
                               stage2.BasalvsNormal = BasalLike.stage2 - Normal.Normal,
                               
                               stage2.HER2vsLumA = HER2enriched.stage2- LuminalA.stage2,
                               stage2.HER2vsLumB = HER2enriched.stage2- LuminalB.stage2,
                               stage2.HER2vsBasal = HER2enriched.stage2- BasalLike.stage2,
                               stage2.HER2vsNormLike = HER2enriched.stage2 - NormalLike.stage2,
                               stage2.HER2vsNormal = HER2enriched.stage2 - Normal.Normal,
                               
                               stage2.NormalLikevsLumA = NormalLike.stage2- LuminalA.stage2,
                               stage2.NormalLikevsLumB = NormalLike.stage2- LuminalB.stage2,
                               stage2.NormalLikevsBasal = NormalLike.stage2- BasalLike.stage2,
                               stage2.NormalLikevsHER2 = NormalLike.stage2- HER2enriched.stage2,
                               stage2.NormalLikevsNormal =NormalLike.stage2 - Normal.Normal,  
                               
                               
                               
                               stage3.LumAvsLumB = LuminalA.stage3 - LuminalB.stage3,     
                               stage3.LumAvsBasal = LuminalA.stage3 - BasalLike.stage3,
                               stage3.LumAvsHER2 = LuminalA.stage3 - HER2enriched.stage3,
                               stage3.LumAvsNormLike = LuminalA.stage3 - NormalLike.stage3,
                               stage3.LumAvsNormal = LuminalA.stage3 - Normal.Normal,
                               
                               stage3.LumBvsLumA = LuminalB.stage3 - LuminalA.stage3,
                               stage3.LumBvsBasal = LuminalB.stage3 - BasalLike.stage3,
                               stage3.LumBvsHER2 = LuminalB.stage3 - HER2enriched.stage3,
                               stage3.LumBvsNormLike = LuminalB.stage3 - NormalLike.stage3,
                               stage3.LumBvsNormal = LuminalB.stage3 - Normal.Normal,
                               
                               stage3.BasalvsLumA = BasalLike.stage3- LuminalA.stage3,
                               stage3.BasalvsLumB = BasalLike.stage3- LuminalB.stage3,
                               stage3.BasalvsHER2 = BasalLike.stage3- HER2enriched.stage3,
                               stage3.BasalvsNormLike = BasalLike.stage3 - NormalLike.stage3,
                               stage3.BasalvsNormal = BasalLike.stage3 - Normal.Normal,
                               
                               stage3.HER2vsLumA = HER2enriched.stage3- LuminalA.stage3,
                               stage3.HER2vsLumB = HER2enriched.stage3- LuminalB.stage3,
                               stage3.HER2vsBasal = HER2enriched.stage3- BasalLike.stage3,
                               stage3.HER2vsNormLike = HER2enriched.stage3 - NormalLike.stage3,
                               stage3.HER2vsNormal = HER2enriched.stage3 - Normal.Normal,
                               
                               stage3.NormalLikevsLumA = NormalLike.stage3- LuminalA.stage3,
                               stage3.NormalLikevsLumB = NormalLike.stage3- LuminalB.stage3,
                               stage3.NormalLikevsBasal = NormalLike.stage3- BasalLike.stage3,
                               stage3.NormalLikevsHER2 = NormalLike.stage3- HER2enriched.stage3,
                               stage3.NormalLikevsNormal =NormalLike.stage3 - Normal.Normal,         
                               
                               
                               stage4.LumAvsLumB = LuminalA.stage4 - LuminalB.stage4,     
                               stage4.LumAvsBasal = LuminalA.stage4 - BasalLike.stage4,
                               stage4.LumAvsHER2 = LuminalA.stage4 - HER2enriched.stage4,
                               stage4.LumAvsNormal = LuminalA.stage4 - Normal.Normal,
                               
                               stage4.LumBvsLumA = LuminalB.stage4 - LuminalA.stage4,
                               stage4.LumBvsBasal = LuminalB.stage4 - BasalLike.stage4,
                               stage4.LumBvsHER2 = LuminalB.stage4 - HER2enriched.stage4,
                               stage4.LumBvsNormal = LuminalB.stage4 - Normal.Normal,
                               
                               stage4.BasalvsLumA = BasalLike.stage4- LuminalA.stage4,
                               stage4.BasalvsLumB = BasalLike.stage4- LuminalB.stage4,
                               stage4.BasalvsHER2 = BasalLike.stage4- HER2enriched.stage4,
                               stage4.BasalvsNormal = BasalLike.stage4 - Normal.Normal,
                               
                               stage4.HER2vsLumA = HER2enriched.stage4- LuminalA.stage4,
                               stage4.HER2vsLumB = HER2enriched.stage4- LuminalB.stage4,
                               stage4.HER2vsBasal = HER2enriched.stage4- BasalLike.stage4,
                               stage4.HER2vsNormal = HER2enriched.stage4 - Normal.Normal,
                               
                               levels = colnames(design)) 
                               


contr.matrix <- makeContrasts(LobularvsNormal = Lobularcarcinoma - Normal,
                              DuctalvsNormal = Ductalcarcinoma - Normal,
                              DuctLobvsNormal = Ductalandlobularmixed - Normal,
                              DuctOthersvsNormal = Ductualmixedwithothers - Normal,
                              MetaplastvsNormal = Metaplasticcarcinoma - Normal,
                              MucinousvsNormal = Mucinousadenocarcinoma - Normal,
                              
                              LobularvsDuctal = Lobularcarcinoma - Ductalcarcinoma,
                              LobularvsDuctLob = Lobularcarcinoma - Ductalandlobularmixed,
                              LobularvsDuctOther = Lobularcarcinoma - Ductualmixedwithothers,
                              LobularvsMetaplast = Lobularcarcinoma - Metaplasticcarcinoma,
                              LobularvsMucinous = Lobularcarcinoma - Mucinousadenocarcinoma,
                              
                              DuctalvsLobular = Ductalcarcinoma - Lobularcarcinoma,
                              DuctalvsDuctLob = Ductalcarcinoma - Ductalandlobularmixed,
                              DuctalvsDuctOthers = Ductalcarcinoma - Ductualmixedwithothers,
                              DuctalvsMetaplast = Ductalcarcinoma - Metaplasticcarcinoma,
                              DuctalvsMucinous = Ductalcarcinoma - Mucinousadenocarcinoma,
                              
                              DuctLobvsLobular = Ductalandlobularmixed - Lobularcarcinoma,
                              DuctLobvsDuctal  = Ductalandlobularmixed - Ductalcarcinoma,
                              DuctLobvsDuctOther =Ductalandlobularmixed - Ductualmixedwithothers,
                              DuctLobvsMetaplast =Ductalandlobularmixed - Metaplasticcarcinoma,
                              DuctLobvsMucinous  =Ductalandlobularmixed - Mucinousadenocarcinoma,
                              
                              DuctOthersvsLobular = Ductualmixedwithothers - Lobularcarcinoma,
                              DuctOthersvsDuctal = Ductualmixedwithothers - Ductalcarcinoma,
                              DuctOthersvsDuctLob = Ductualmixedwithothers -  - Ductalandlobularmixed,
                              DuctOthersvsMetaplast = Ductualmixedwithothers - - Metaplasticcarcinoma,
                              DuctOthersvsMucinous  = Ductualmixedwithothers -  Mucinousadenocarcinoma,
                              
                              MetaplastvsLobular = Metaplasticcarcinoma -  Lobularcarcinoma,
                              MetaplastvsDuctal = Metaplasticcarcinoma - Ductalcarcinoma,
                              MetaplastvsDuctLob = Metaplasticcarcinoma -  Ductalandlobularmixed,
                              MetaplastvsDuctOther = Metaplasticcarcinoma -  Ductualmixedwithothers,
                              MetaplastvsMucinous = Metaplasticcarcinoma  -   Mucinousadenocarcinoma,
                              
                              MucinousvsLobular = Mucinousadenocarcinoma -  Lobularcarcinoma,
                              MucinousvsDuctal = Mucinousadenocarcinoma -Ductalcarcinoma,
                              MucinousvsDuctLob = Mucinousadenocarcinoma -  Ductalandlobularmixed,
                              MucinousvsDuctOther = Mucinousadenocarcinoma - Ductualmixedwithothers,
                              MucinousvsMetaplast = Mucinousadenocarcinoma - Metaplasticcarcinoma,
                              
                              levels = colnames(design)) 
  

contr.matrix <- makeContrasts( Stage1vsNorm = stage1 - Normal,  
                               Stage2vsNorm = stage2 - Normal,
                               Stage3vsNorm = stage3 - Normal,
                               Stage4vsNorm = stage4- Normal,
                               
                               Stage1vsStage2 = stage1 - stage2,
                               Stage1vsStage3 = stage1 - stage3,
                               Stage1vsStage4 = stage1 - stage4,
                               
                               Stage2vsStage1 = stage2 - stage1,
                               Stage2vsStage3 = stage2 - stage3,
                               Stage2vsStage4 = stage2 - stage4,
                               
                               Stage3vsStage4 = stage3 - stage4,
                               Stage3vsStage2 = stage3 - stage2,
                               Stage3vsStage1 = stage3 - stage1,
                               
                               levels = colnames(design)) 

######
contr.matrix <- makeContrasts( LumAvsLumB = LuminalA - LuminalB,     #unique comparisons 15
                               LumAvsBasal = LuminalA - BasalLike,
                               LumAvsHER2 = LuminalA - HER2enriched,
                               LumAvsNormLike = LuminalA - NormalLike,
                               LumAvsNormal = LuminalA - Normal,
                               #
                               LumBvsLumA = LuminalB - LuminalA,
                               LumBvsBasal = LuminalB - BasalLike,
                               LumBvsHER2 = LuminalB - HER2enriched,
                               LumBvsNormLike = LuminalB - NormalLike,
                               LumBvsNormal = LuminalB - Normal,
                               #
                               BasalvsLumA = BasalLike- LuminalA,
                               BasalvsLumB = BasalLike- LuminalB,
                               BasalvsHER2 = BasalLike- HER2enriched,
                               BasalvsNormLike = BasalLike - NormalLike,
                               BasalvsNormal = BasalLike - Normal,
                               #
                               HER2vsLumA = HER2enriched- LuminalA,
                               HER2vsLumB = HER2enriched- LuminalB,
                               HER2vsBasal = HER2enriched- BasalLike,
                               HER2vsNormLike = HER2enriched - NormalLike,
                               HER2vsNormal = HER2enriched - Normal,
                               #
                               NormalLikevsLumA = NormalLike- LuminalA,
                               NormalLikevsLumB = NormalLike- LuminalB,
                               NormalLikevsBasal = NormalLike- BasalLike,
                               NormalLikevsHER2 = NormalLike- HER2enriched,
                               NormalLikevsNormal =NormalLike - Normal,
                               
                               levels = colnames(design)) 


#####


# eBayes
DEA_limmaVoom <- function(y, design, contr.matrix=NULL) {
  #the voom transformation is applied to the normalized and filtered DGEList object:
  #Use voom() to convert the read counts to log2-cpm, with associated weights, 
  #ready for linear modelling:
  #par(mfrow=c(2,2))
  v <- voom(y, design, plot=F); print ("done Voom")
  
  #After this, the usual limma pipelines for differential expression can be applied
  #fit a separate model to the expression values for each gene
  vfit <- lmFit(v, design); print ("done lmFit")
  if (!is.null(contr.matrix)){
    print ("Using contrast matrix... ")
    vfit <- contrasts.fit(vfit, contrasts=contr.matrix)  
  }
  #empirical Bayes moderation is carried out by borrowing information across all genes
  #to obtain more precise estimates of gene-wise variability
  
  efit <- eBayes(vfit); print ("done eBayes")
  
  #The model's residual variances are plotted against average expression values
  # plotSA(efit,main="plotSA()")
  
  #testing FDR p-val distribution #always check p-val disrubution to make sure FDR worked
  #hist(as.vector(as.matrix(efit$p.value[,c(1:ncol(efit$p.value))])), col="red", main="P-values before MTC") #NBBBB not sure if this is the right p-vals
  par(mfrow=c(1,1))
  
  return(list(fit=efit, v=v)) #change if choose stricter option
}  
DEA_MTC_save <-function(fit, my_coef, logFC, FDR){
  
  #For a quick look at differential expression levels, the number of significantly up-
  #and down-regulated genes can be summarised in a table.
  dt<-decideTests(fit, p.value=FDR, lfc= logFC, adjust.method="BH")
  
  
  print (my_coef)
  print (summary(dt)[,c(my_coef)])
  
  
  tt <- topTable(fit, coef=my_coef, adjust='BH', n=Inf)
  #tt <- topTreat(fit, coef=my_coef, adjust='BH', n=Inf)
  
  #if toy want acceess to df of genes
  #up <- data.frame(rownames(tt[tt$logFC >= logFC & tt$adj.P.Val < FDR, ]))
  #down <- data.frame(rownames(tt[tt$logFC <= -logFC & tt$adj.P.Val < FDR, ]))
  #colnames(up) <- as.character("up")
  #colnames(down) <- as.character("down")
  
  #tt <- subset(tt,abs(tt$logFC) >= logFC & tt$adj.P.Val < FDR)
  index.up <- which(tt$logFC >= logFC & tt$adj.P.Val < FDR)
  index.down <- which(tt$logFC <= -logFC & tt$adj.P.Val < FDR)
  direction <- c()
  direction[index.up] <- "up"
  direction[index.down] <- "down"
  direction[!(1:nrow(tt) %in% union(index.up,index.down))] <- "no DE"
  tt <- cbind(tt,direction)
  #final <- list(up, down)
  
  #save for the specified coefficient
  #write.csv(tt, "both_PAM_DEG_23-04.csv", quote = FALSE)
  up <- data.frame(rownames(tt[tt$direction == "up", ]))
  down <- data.frame(rownames(tt[tt$direction == "down", ]))
  
  setwd("~/Bioinformatics MSc UCPH/0_MasterThesis/TCGAbiolinks/CBL_scripts/data/DEA/")
  #write.table(up, paste0(my_coef,"_up.txt"), sep = "\t",col.names = FALSE,row.names = FALSE, quote = FALSE)
  #write.table(down, paste0(my_coef,"_down.txt"),sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)
  
  return(list(tt=tt, dt=dt))
}

DEA_limmaVoom_out<-DEA_limmaVoom(y,design,contr.matrix)
fit<-DEA_limmaVoom_out$fit


#Treat
#DEA_limmaVoom <- function(y, design, contr.matrix=NULL) {
  #the voom transformation is applied to the normalized and filtered DGEList object:
  #Use voom() to convert the read counts to log2-cpm, with associated weights, 
  #ready for linear modelling:
  #par(mfrow=c(2,2))
  v <- voom(y, design, plot=FALSE); print ("done Voom")
  
  #After this, the usual limma pipelines for differential expression can be applied
  #fit a separate model to the expression values for each gene
  vfit <- lmFit(v, design); print ("done lmFit")
  if (!is.null(contr.matrix)){
    print ("Using contrast matrix... ")
    vfit <- contrasts.fit(vfit, contrasts=contr.matrix)  
  }
  #empirical Bayes moderation is carried out by borrowing information across all genes
  #to obtain more precise estimates of gene-wise variability
  
  #efit <- eBayes(vfit); print ("done eBayes")
  
  # for stricter filtering on significance
  
  # The treat method can be used to calculate p-values from
  # empirical Bayes moderated t-statistics with a minimum log-FC requirement
  #when testing requires genes to have a log-FC that is significantly 
  #greater than 1 (equivalent to a 2-fold difference between cell types on the original scale)
  tfit <- treat(vfit, lfc=1) ; print ("done treat")
  #dt <- decideTests(tfit)
  
  
  #The model's residual variances are plotted against average expression values
  # plotSA(efit,main="plotSA()")
  
  #testing FDR p-val distribution #always check p-val disrubution to make sure FDR worked
  #hist(as.vector(as.matrix(efit$p.value[,c(1:ncol(efit$p.value))])), col="red", main="P-values before MTC") #NBBBB not sure if this is the right p-vals
  par(mfrow=c(1,1))
  
  return(list(fit=tfit, v=v)) #change if choose stricter option
}  
#DEA_MTC_save <-function(fit, my_coef, logFC, FDR){
  
  #For a quick look at differential expression levels, the number of significantly up-
  #and down-regulated genes can be summarised in a table.
  dt<-decideTests(fit, p.value=FDR, lfc= logFC, adjust.method="BH")
  
  
  print (my_coef)
  print (summary(dt)[,c(my_coef)])
  
  
  #tt <- topTable(fit, coef=my_coef, adjust='BH', n=Inf)
  tt <- topTreat(fit, coef=my_coef, adjust='BH', n=Inf)
  
  #if toy want acceess to df of genes
  #up <- data.frame(rownames(tt[tt$logFC >= logFC & tt$adj.P.Val < FDR, ]))
  #down <- data.frame(rownames(tt[tt$logFC <= -logFC & tt$adj.P.Val < FDR, ]))
  #colnames(up) <- as.character("up")
  #colnames(down) <- as.character("down")
  
  #tt <- subset(tt,abs(tt$logFC) >= logFC & tt$adj.P.Val < FDR)
  index.up <- which(tt$logFC >= logFC & tt$adj.P.Val < FDR)
  index.down <- which(tt$logFC <= -logFC & tt$adj.P.Val < FDR)
  direction <- c()
  direction[index.up] <- "up"
  direction[index.down] <- "down"
  direction[!(1:nrow(tt) %in% union(index.up,index.down))] <- "no DE"
  tt <- cbind(tt,direction)
  #final <- list(up, down)
  
  #save for the specified coefficient
  #write.csv(tt, "both_PAM_DEG_23-04.csv", quote = FALSE)
  up <- data.frame(rownames(tt[tt$direction == "up", ]))
  down <- data.frame(rownames(tt[tt$direction == "down", ]))
  
  setwd("~/Bioinformatics MSc UCPH/0_MasterThesis/TCGAbiolinks/CBL_scripts/data/DEA/testing")
  write.table(up, paste0(my_coef,"_up.txt"), sep = "\t",col.names = FALSE,row.names = FALSE, quote = FALSE)
  write.table(down, paste0(my_coef,"_down.txt"),sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)
  
  return(list(tt=tt, dt=dt))
}
#DEA_limmaVoom_out<-DEA_limmaVoom(y,design,contr.matrix)
#fit<-DEA_limmaVoom_out$fit



#set parameters!
set_logFC=1
set_FDR=0.05

## decideTests        # here coef does not matter
dt <-DEA_MTC_save(fit, 2, set_logFC , set_FDR)$dt #coef can be anything, report dt for all anyway
summary(dt)

contast_DE_list<-list()
for (i in colnames(dt)){
  this_coef<-DEA_MTC_save(fit, i, set_logFC , set_FDR)$tt
  contast_DE_list[[i]]<-this_coef
}

names(contast_DE_list)



#### visualisation of limma results####

# PVALUES AFTER MTC

getVarname<-function(var){
  deparse(substitute(var))
}

### STAGES
plot1 <- qplot(data=Stage1vsNorm_tt  , main = "Stage1vsNorm_tt  ", x=P.Value, geom="histogram", color=I("black"), fill=I("hotpink"), binwidth=0.05)
plot2 <- qplot(data=Stage2vsNorm_tt  , main = "Stage2vsNorm_tt  ", x=P.Value, geom="histogram", color=I("black"), fill=I("hotpink"), binwidth=0.05)
plot3 <- qplot(data=Stage3vsNorm_tt  , main = "Stage3vsNorm_tt  ", x=P.Value, geom="histogram", color=I("black"), fill=I("hotpink"), binwidth=0.05)
plot4 <- qplot(data=Stage4vsNorm_tt  , main = "Stage4vsNorm_tt  ", x=P.Value, geom="histogram", color=I("black"), fill=I("hotpink"), binwidth=0.05)
plot5 <- qplot(data=Stage1vsStage2_tt, main = "Stage1vsStage2_tt", x=P.Value, geom="histogram", color=I("black"), fill=I("hotpink"), binwidth=0.05)
plot6 <- qplot(data=Stage1vsStage3_tt, main = "Stage1vsStage3_tt", x=P.Value, geom="histogram", color=I("black"), fill=I("hotpink"), binwidth=0.05)
plot7 <- qplot(data=Stage1vsStage4_tt, main = "Stage1vsStage4_tt", x=P.Value, geom="histogram", color=I("black"), fill=I("hotpink"), binwidth=0.05)
plot8 <- qplot(data=Stage2vsStage3_tt, main = "Stage2vsStage3_tt", x=P.Value, geom="histogram", color=I("black"), fill=I("hotpink"), binwidth=0.05)
plot9 <- qplot(data=Stage2vsStage4_tt, main = "Stage2vsStage4_tt", x=P.Value, geom="histogram", color=I("black"), fill=I("hotpink"), binwidth=0.05)
plot10<- qplot(data=Stage3vsStage4_tt, main = "Stage3vsStage4_tt", x=P.Value, geom="histogram", color=I("black"), fill=I("hotpink"), binwidth=0.05)

require(gridExtra)
grid.arrange(plot1, plot2,plot3,plot4,plot5,plot6,plot7,plot8,plot9,plot10, ncol=3)
#A flat p-value histogram does not indicate any problem of the methods that operate 
#on p-values but suggests the disappointing result that very few, if any, genes are
#differentially expressed. Depletion of small p-values can indicate the presence of 
#confounding hidden variables - "batch effect2

#GROUP1
plot1 <- qplot(data=Basal_st1vs4_tt      , main = "Basal_st1vs4_tt", x=P.Value, geom="histogram", color=I("black"), fill=I("hotpink"), binwidth=0.05)
plot2 <- qplot(data=LumA_st1vs4_tt       , main = "LumA_st1vs4_tt ", x=P.Value, geom="histogram", color=I("black"), fill=I("hotpink"), binwidth=0.05)
plot3 <- qplot(data=allBasal_tt          , main = "allBasal_tt", x=P.Value, geom="histogram", color=I("black"), fill=I("hotpink"), binwidth=0.05)
plot4 <- qplot(data=allLumA_tt           , main = "allLumA_tt ", x=P.Value, geom="histogram", color=I("black"), fill=I("hotpink"), binwidth=0.05)
plot5 <- qplot(data=Basal_stage1_tt      , main = "Basal_stage1_tt", x=P.Value, geom="histogram", color=I("black"), fill=I("hotpink"), binwidth=0.05)
plot6 <- qplot(data=Basal_stage4_tt      , main = "Basal_stage4_tt", x=P.Value, geom="histogram", color=I("black"), fill=I("hotpink"), binwidth=0.05)
plot7 <- qplot(data=LumA_stage1_tt       , main = "LumA_stage1_tt", x=P.Value, geom="histogram", color=I("black"), fill=I("hotpink"), binwidth=0.05)
plot8 <- qplot(data=LumA_stage4_tt       , main = "LumA_stage4_tt", x=P.Value, geom="histogram", color=I("black"), fill=I("hotpink"), binwidth=0.05)
plot9 <- qplot(data=BasalVSLumA_stage1_tt, main = "BasalVSLumA_stage1_tt", x=P.Value, geom="histogram", color=I("black"), fill=I("hotpink"), binwidth=0.05)
plot10<- qplot(data=BasalVSLumA_stage4_tt, main = "BasalVSLumA_stage4_tt", x=P.Value, geom="histogram", color=I("black"), fill=I("hotpink"), binwidth=0.05)
plot11 <- qplot(data=average_stage4_tt, main = "average_stage4_tt", x=P.Value, geom="histogram", color=I("black"), fill=I("hotpink"), binwidth=0.05)
plot12<- qplot(data=LumA_Basal_diff_tt, main = "LumA_Basal_diff_tt", x=P.Value, geom="histogram", color=I("black"), fill=I("hotpink"), binwidth=0.05)


require(gridExtra)
grid.arrange(plot1, plot2,plot3,plot4,plot5,plot6,plot7,plot8,plot9,plot10, plot11, plot12, ncol=3)


################



# Get all genes
getGeneNames<-function(comparison, direction){
  #comparison<-comparison[1:100,] #take top 100 lowest Pval (maybe logFC is btter??)
  if (direction != 'both'){
    genes<-rownames(comparison[which(comparison$direction==direction),])
  } else{
    genes_up<-rownames(comparison[which(comparison$direction=='up'),])
    genes_down<-rownames(comparison[which(comparison$direction=='down'),])
    genes<-append(genes_up,genes_down)
  }    
  return(as.character(genes))
}

## getting autophagy genes DE in a supplies gene list
sharedWithAuto <-function (gene_list){
  autophagy_genes<- as.vector(read.table("../autopahagy_genes.txt", as.is = T, header = FALSE))$V1
  shared <- intersect(autophagy_genes,gene_list)
  print(paste0("Total number of genes in the list : ", length(gene_list)))
  print(paste0("Autophagy genes among these: ", length(shared)))
  return(shared)
  
}
sharedWithAutoCORE <-function (gene_list){
  autophagy_genes<- as.vector(read.table("../autophagic_core.txt", as.is = T, header = FALSE))$V1
  shared <- intersect(autophagy_genes,gene_list)
  print(paste0("Total number of genes in the list : ", length(gene_list)))
  print(paste0("Autophagy CORE genes among these: ", length(shared)))
  return(shared)
  
}
sharedWithAutoTF <-function (gene_list){
  autophagy_genes<- as.vector(read.table("../transcription_factors.txt", as.is = T, header = FALSE))$V1
  shared <- intersect(autophagy_genes,gene_list)
  print(paste0("Total number of genes in the list : ", length(gene_list)))
  print(paste0("Autophagy TF genes among these: ", length(shared)))
  return(shared)
  
}



#create a df to hold DE genes counts 
genesDEdf<-data.frame(matrix( NA,  nrow = length(names(contast_DE_list)), ncol = 3))
rownames(genesDEdf)<-names(contast_DE_list)
colnames(genesDEdf)<-c("up","down","both")
head(genesDEdf)


#create a df to hold DE AUTOPHAGY genes counts 
genesDE_AUTOdf<-data.frame(matrix( NA,  nrow = length(names(contast_DE_list)), ncol = 3))
rownames(genesDE_AUTOdf)<-names(contast_DE_list)
colnames(genesDE_AUTOdf)<-c("up","down","both")
head(genesDE_AUTOdf)


#create a df to hold DE AUTOPHAGY CORE genes counts 
genesDE_AUTOCOREdf<-data.frame(matrix( NA,  nrow = length(names(contast_DE_list)), ncol = 3))
rownames(genesDE_AUTOCOREdf)<-names(contast_DE_list)
colnames(genesDE_AUTOCOREdf)<-c("up","down","both")
head(genesDE_AUTOCOREdf)

#create a df to hold DE AUTOPHAGY CORE genes counts 
genesDE_AUTOTFdf<-data.frame(matrix( NA,  nrow = length(names(contast_DE_list)), ncol = 3))
rownames(genesDE_AUTOTFdf)<-names(contast_DE_list)
colnames(genesDE_AUTOTFdf)<-c("up","down","both")
head(genesDE_AUTOTFdf)


all_genesDE <- c() # for reference, count how many genes were DE in all comparisons
autophagy_genesDE <- c()
autophagyCORE_genesDE <- c()
autophagyTF_genesDE <- c()

# fill in df with numbers
for (i in 1:length(rownames(genesDEdf))){
  this_contrast<- rownames(genesDEdf)[i]
  
  # get numbers for all DE
  genesDEdf[i,"up"]<- length(getGeneNames(contast_DE_list[[this_contrast]],'up'))
  genesDEdf[i,"down"]<- length(getGeneNames(contast_DE_list[[this_contrast]],'down'))
  genesDEdf[i,"both"]<- length(getGeneNames(contast_DE_list[[this_contrast]],'both'))
  
  #get numbers for autophagy DE only 
  genesDE_AUTOdf[i,"up"]<- length(sharedWithAuto(getGeneNames(contast_DE_list[[this_contrast]],'up')))
  genesDE_AUTOdf[i,"down"]<- length(sharedWithAuto(getGeneNames(contast_DE_list[[this_contrast]],'down')))
  genesDE_AUTOdf[i,"both"]<- length(sharedWithAuto(getGeneNames(contast_DE_list[[this_contrast]],'both')))
  
  #get numbers for autophagy CORE DE only 
  genesDE_AUTOCOREdf[i,"up"]<- length(sharedWithAutoCORE(getGeneNames(contast_DE_list[[this_contrast]],'up')))
  genesDE_AUTOCOREdf[i,"down"]<- length(sharedWithAutoCORE(getGeneNames(contast_DE_list[[this_contrast]],'down')))
  genesDE_AUTOCOREdf[i,"both"]<- length(sharedWithAutoCORE(getGeneNames(contast_DE_list[[this_contrast]],'both')))
  
  #get numbers for autophagy TF DE only 
  genesDE_AUTOTFdf[i,"up"]<- length(sharedWithAutoTF(getGeneNames(contast_DE_list[[this_contrast]],'up')))
  genesDE_AUTOTFdf[i,"down"]<- length(sharedWithAutoTF(getGeneNames(contast_DE_list[[this_contrast]],'down')))
  genesDE_AUTOTFdf[i,"both"]<- length(sharedWithAutoTF(getGeneNames(contast_DE_list[[this_contrast]],'both')))
  
  all_genesDE <- unique(sort(append(all_genesDE, getGeneNames(contast_DE_list[[this_contrast]],'both'))))
  autophagy_genesDE <- unique(sort(append(autophagy_genesDE, sharedWithAuto(getGeneNames(contast_DE_list[[this_contrast]],'both')))))
  autophagyCORE_genesDE <- unique(sort(append(autophagyCORE_genesDE, sharedWithAutoCORE(getGeneNames(contast_DE_list[[this_contrast]],'both')))))
  autophagyTF_genesDE <- unique(sort(append(autophagyTF_genesDE, sharedWithAutoTF(getGeneNames(contast_DE_list[[this_contrast]],'both')))))
  
  
}
length(all_genesDE)
length(autophagy_genesDE)
length(autophagyCORE_genesDE)
length(autophagyTF_genesDE)


# when group1 is in the model, for PAM50: 6818 / 412/ 52 /48

save(autophagy_genesDE, file="autophagy_genesDE.rda")



###### saving data for Enrichment Analysis

###PAM50
#  x   is cpm>1>=50 Treat
# 2 is with cpm>2>=19 + Treat
# 3 is with cpm>2>=19 + ebayes ----- this is the one I use now!

setwd("~/Bioinformatics MSc UCPH/0_MasterThesis/TCGAbiolinks/CBL_scripts/data/DEA/EnrichmentAnalysis")

#PAM50
save(genesDEdf, file="PAM50_DE_genes_numbers3.rda")
save(genesDE_AUTOdf, file="PAM50_DE_AUTO_genes_numbers3.rda")
save(genesDE_AUTOCOREdf, file="PAM50_DE_AUTOCORE_genes_numbers3.rda")
save(genesDE_AUTOTFdf, file="PAM50_DE_AUTOTF_genes_numbers3.rda")


### STAGES
save(genesDEdf, file="stages_DE_genes_numbers3.rda")
save(genesDE_AUTOdf, file="stages_DE_AUTO_genes_numbers3.rda")
save(genesDE_AUTOCOREdf, file="stages_DE_AUTOCORE_genes_numbers3.rda")
save(genesDE_AUTOTFdf, file="stages_DE_AUTOTF_genes_numbers3.rda")

### morphology
save(genesDEdf, file="morph_DE_genes_numbers3.rda")
save(genesDE_AUTOdf, file="morph_DE_AUTO_genes_numbers3.rda")
save(genesDE_AUTOCOREdf, file="morph_DE_AUTOCORE_genes_numbers3.rda")
save(genesDE_AUTOTFdf, file="morph_DE_AUTOTF_genes_numbers3.rda")


#Group1
save(genesDEdf, file="Group1_DE_genes_numbers3.rda")
save(genesDE_AUTOdf, file="Group1_DE_AUTO_genes_numbers3.rda")
save(genesDE_AUTOCOREdf, file="Group1_DE_AUTOCORE_genes_numbers3.rda")
save(genesDE_AUTOTFdf, file="Group1_DE_AUTOTF_genes_numbers3.rda")

#Group2
save(genesDEdf, file="Group2_DE_genes_numbers3.rda")
save(genesDE_AUTOdf, file="Group2_DE_AUTO_genes_numbers3.rda")
save(genesDE_AUTOCOREdf, file="Group2_DE_AUTOCORE_genes_numbers3.rda")
save(genesDE_AUTOTFdf, file="Group2_DE_AUTOTF_genes_numbers3.rda")



#checking
allDE<-get(load("PAM50_DE_genes_numbers3.rda"))
dim(allDE)

stop()




####### HEATMAPS OF DEA #########



#to_show_in_hm <- all_genesDE
to_show_in_hm <- autophagy_genesDE
to_show_in_hm <- autophagyCORE_genesDE
to_show_in_hm <- autophagyTF_genesDE

# Get expression values
v<-DEA_limmaVoom_out$v  
save(v, file="v_expression_for_hm.rda")

all_geneExp <- v$E[rownames(v$E) %in% to_show_in_hm , ]
dim(all_geneExp)
#View(all_geneExp)

#write.table(all_L, paste0(my.dir, my.pattern, "_DE_Limma.txt"), sp = "\t", quote = FALSE)



#### genes classes setup #####

g_functions <- get(load("../autophagy_functions.rda"))       #### add comments!
rownames(g_functions)<-g_functions$genes 

ok_autop_genes<-rownames(dataSE)
ok_g_functions<-g_functions[c(ok_autop_genes),]

thirteen_cols<- c("#3579b4","#c8c049","#8996da","#ee5345",
                  "#c84297","#43ea3b","#50a376","#281340",
                  "#6e5c41","#fd0d32","#f19832","#94f5d2","#b1f555")
ok_g_functions$gene_functions=factor(ok_g_functions$gene_functions,
                                     levels=c("multifunction","lipid","phosphatidyl","endo_exosomes",    
                                              "transport","rabs","docking_fusion","mito","autoph_core",
                                              "transcr_factors","receptors_ligands", "mTOR_induction","lysosome"))
names(thirteen_cols)<-levels(ok_g_functions$gene_functions)
ok_g_functions$genes <-NULL

###### patient classes setup #####

#get only interesting cols
names(samplesMatrix)
hm.design<-subset(samples.matrix[,c('tumourTypes','tumourStages','PAM50')])#, 'tss'
row.names(hm.design)<-samples.matrix$barcode

#recorder annotations
hm.design$PAM50 = factor(hm.design$PAM50, levels = c("Basal-like", "HER2-enriched", "Luminal A", "Luminal B", "Normal-like", "Normal"))
hm.design$tumourStages = factor(hm.design$tumourStages, levels = c("stage1", "stage2", "stage3", "stage4", "unknown")) 

hm.design$tumourTypes = factor(hm.design$tumourTypes,  
                            levels = c("Normal",  "Mucinous adenocarcinoma", 
                                       "Ductal carcinoma","Lobular carcinoma", "Ductual mixed with others ",
                                       "Ductal and lobular mixed", "Metaplastic carcinoma" ,"Other")) #other is a mix of few samples of different ones

# colours
PAM50Col <- c( "#E69F00", "#0072B2", "#F0E442","red3", "#CC79A7","#009E73")##D55E00->altred
#PAM50Col= c( "turquoise4", "darkorange2", "slateblue2", "palevioletred2", "gold2", "olivedrab3")

#prev tumourTypCol = c( "#bebada", "#8dd3c7","#b3de69",  "#fb8072", "#80b1d3", "#fdb462","white")
tumourTypCol = c("#009E73",  "#0072B2", "#D55E00","#56B4E9","#E69F00","#F0E442","#CC79A7", "white")

tumourStgCol = c( "red","green","blue", "yellow2", "white")

names(PAM50Col)<-levels(hm.design$PAM50)
names(tumourTypCol)<-levels(hm.design$tumourTypes)
names(tumourStgCol)<-levels(hm.design$tumourStages)

annColour <-list(
  PAM50=PAM50Col,
  tumourTypes=tumourTypCol,
  tumourStages=tumourStgCol,
  gene_functions=thirteen_cols
)



######
pheatmap::pheatmap(mat = as.matrix(all_geneExp), color = brewer.pal(name = "YlGnBu", n = 9),
                   clustering_distance_rows = 'manhattan', 
                   clustering_distance_cols = 'manhattan', 
                   #scale="row",
                   annotation_col=hm.design,
                   annotation_colors = annColour,
                   annotation_row=ok_g_functions,
                   cluster_cols = T, cluster_rows = T, 
                   show_rownames = T,show_colnames = F,
                   fontsize = 5)










#Genes that are DE in multiple comparisons can be extracted using the results from decideTests

topTable(efit, coef = c(1,2))
topTable(efit, coef = "conditioncancer")
topTable(efit, coef = "conditioncancer", n = 10, sort = "p", p = 0.05)[,-1] 

## see commom in 2 conditions 
de.common <- which(dt[,5]!=0 & dt[,10]!=0) 
length(de.common)



library(VennDiagram)
par(mfrow=c(1,1))
# with model 0+PAM50+condition all dt are the same
colnames(summary(dt))

plot(venneuler(dt[,c(5,9,12,14,15)]))

vennDiagram(dt[,c(5,9,12,14,15)], circle.col=c("turquoise", "salmon", "blue", "yellow", "green"))
vennDiagram(dt[,c(5,9,15)], circle.col=c("turquoise", "salmon", "yellow"),cex=0.8)

vennDiagram(dt[,c(5:8)], circle.col=c("turquoise", "salmon", "yellow" , "blue"),cex=0.8)


#with intercept : get same count as without intercept for VS normal comparison
vennDiagram(dt[,c(2:6)], circle.col=c("turquoise", "salmon", "blue", "yellow", "green"))


#group1
vennDiagram(dt[,c(5,7)], circle.col=c("turquoise", "salmon"), cex=0.9) 
vennDiagram(dt[,c(6,8)], circle.col=c("turquoise", "salmon"), cex=0.9)
vennDiagram(dt[,c(5,7,9)], circle.col=c("turquoise", "salmon", "yellow"),cex=0.8)
vennDiagram(dt[,c(6,8,10)], circle.col=c("turquoise", "salmon", "yellow"),cex=0.8)

vennDiagram(dt[,c(5,6)], circle.col=c("turquoise", "salmon"), cex=0.9)
vennDiagram(dt[,c(7,8)], circle.col=c("turquoise", "salmon"), cex=0.9)
vennDiagram(dt[,c(5,6,1)], circle.col=c("turquoise", "salmon", "yellow"),cex=0.8)
vennDiagram(dt[,c(5,7,9)], circle.col=c("turquoise", "salmon", "yellow"),cex=0.8)

vennDiagram(dt[,c(1,2,11)], circle.col=c("turquoise", "salmon", "yellow"), cex=0.9)
vennDiagram(dt[,c(3,4,12)], circle.col=c("blue","green","pink"), cex=0.9) 


head(efit$genes$SYMBOL[de.common])

### examining DE genes 

head(Stage1vsNorm_tt[Stage1vsNorm_tt$direction != "no DE", ], n=10)
head(Stage2vsNorm_tt[Stage2vsNorm_tt$direction != "no DE", ], n=10)
head(Stage3vsNorm_tt[Stage3vsNorm_tt$direction != "no DE", ], n=10)
head(Stage4vsNorm_tt[Stage4vsNorm_tt$direction != "no DE", ], n=20)
head(Stage1vsStage2_tt[Stage1vsStage2_tt$direction != "no DE", ])
head(Stage1vsStage3_tt[Stage1vsStage3_tt$direction != "no DE", ])
head(Stage1vsStage4_tt[Stage1vsStage4_tt$direction != "no DE", ])
head(Stage2vsStage3_tt[Stage2vsStage3_tt$direction != "no DE", ])
head(Stage2vsStage4_tt[Stage2vsStage4_tt$direction != "no DE", ])
head(Stage3vsStage4_tt[Stage3vsStage4_tt$direction!="no DE", ])



## Useful graphical representations of differential expression results
#To summarise results for all genes visually, mean-difference plots, which display log-FCs
#from the linear model fit against the average log-CPM values can be generated using 
#the plotMD function, with the differenti3ally expressed genes highlighted.

colnames(dt)
summary(dt)
par(mfrow=c(1,5))

#vs Normal
plotMD(fit, column=5, status=dt[,5], main=colnames(fit)[5], ylim=c(-10,10))
plotMD(fit, column=9, status=dt[,9], main=colnames(fit)[9], ylim=c(-10,10))
plotMD(fit, column=12, status=dt[,12], main=colnames(fit)[12], ylim=c(-10,10))
plotMD(fit, column=14, status=dt[,14], main=colnames(fit)[14], ylim=c(-10,10))
plotMD(fit, column=15, status=dt[,15], main=colnames(fit)[15], ylim=c(-10,10))

#vs LumA

plotMD(fit, column=1, status=dt[,1], main=colnames(fit)[1], ylim=c(-10,10))
plotMD(fit, column=2, status=dt[,2], main=colnames(fit)[2], ylim=c(-10,10))
plotMD(fit, column=3, status=dt[,3], main=colnames(fit)[3], ylim=c(-10,10))
plotMD(fit, column=4, status=dt[,4], main=colnames(fit)[4], ylim=c(-10,10))
plotMD(fit, column=5, status=dt[,5], main=colnames(fit)[5], ylim=c(-10,10))

#vs LumB
plotMD(fit, column=1, status=dt[,1], main=colnames(fit)[1], ylim=c(-10,10))
plotMD(fit, column=6, status=dt[,6], main=colnames(fit)[6], ylim=c(-10,10))
plotMD(fit, column=7, status=dt[,7], main=colnames(fit)[7], ylim=c(-10,10))
plotMD(fit, column=8, status=dt[,8], main=colnames(fit)[8], ylim=c(-10,10))
plotMD(fit, column=9, status=dt[,9], main=colnames(fit)[9], ylim=c(-10,10))

#vs Basal
plotMD(fit, column=2, status=dt[,2], main=colnames(fit)[2], ylim=c(-10,10))
plotMD(fit, column=6, status=dt[,6], main=colnames(fit)[6], ylim=c(-10,10))
plotMD(fit, column=10, status=dt[,10], main=colnames(fit)[10], ylim=c(-10,10))
plotMD(fit, column=11, status=dt[,11], main=colnames(fit)[11], ylim=c(-10,10))
plotMD(fit, column=12, status=dt[,12], main=colnames(fit)[12], ylim=c(-10,10))

#vs HER2
plotMD(fit, column=3, status=dt[,3], main=colnames(fit)[3], ylim=c(-10,10))
plotMD(fit, column=7, status=dt[,7], main=colnames(fit)[7], ylim=c(-10,10))
plotMD(fit, column=9, status=dt[,9], main=colnames(fit)[9], ylim=c(-10,10))
plotMD(fit, column=13, status=dt[,13], main=colnames(fit)[13], ylim=c(-10,10))
plotMD(fit, column=14, status=dt[,14], main=colnames(fit)[14], ylim=c(-10,10))

#vs NormalLike
plotMD(fit, column=4, status=dt[,4], main=colnames(fit)[4], ylim=c(-10,10))
plotMD(fit, column=8, status=dt[,8], main=colnames(fit)[8], ylim=c(-10,10))
plotMD(fit, column=11, status=dt[,11], main=colnames(fit)[11], ylim=c(-10,10))
plotMD(fit, column=13, status=dt[,13], main=colnames(fit)[13], ylim=c(-10,10))
plotMD(fit, column=15, status=dt[,15], main=colnames(fit)[15], ylim=c(-10,10))



#MAplot
ggplot(BasalvsNormal_tt, aes(x=AveExpr, y=logFC, col=adj.P.Val<0.05)) +
  geom_point(alpha=0.25) + 
  geom_hline(aes(yintercept=2), col="blue") + 
  geom_hline(aes(yintercept=-2), col="blue") + 
  scale_color_manual(values=c("black", "red")) +
  ylim(-10,10)+
  theme(legend.position="none")+
  ggtitle(paste0(getVarname(BasalvsNormal_tt)))+
  xlab("Average Log Expression")

#volcano version 1
ggplot(LumBvsNormal_tt, aes(x=logFC, y=-log10(P.Value), col=adj.P.Val<0.05)) +
  geom_point(alpha=0.25) + 
  geom_vline(aes(xintercept=2), col="blue") + 
  geom_vline(aes(xintercept=-2), col="blue") + 
  scale_color_manual(values=c("black", "red")) + 
  theme(legend.position="none")


#volcano version 2
par(mfrow=c(1,1))
this_data =LumBvsNormal_tt
head(this_data)
this_data$Genes <-rownames(this_data)

for (i in 1:length(this_data$Genes)){
  if (this_data$Genes[i] %in% autophagy_genes){
    this_data$autophagy[i] <- this_data$Genes[i]
  }else{
    this_data$autophagy[i] <-'x'
  }
}


# Make a basic volcano plot
with(this_data, plot(logFC, -log10(P.Value), pch=20, main="Volcano plot: Luminal B vs Normal" ))#, xlim=c(-3,4), ylim=c(0,5)) )

# Add colored points: red if padj<0.05, orange of log2FC>1, green if both)
with(subset(this_data, adj.P.Val<0.05 ), points(logFC, -log10(P.Value), pch=20, col="blue"))
with(subset(this_data, abs(logFC)>2), points(logFC, -log10(P.Value), pch=20, col="orange"))
with(subset(this_data, adj.P.Val<.05 & abs(logFC)>2), points(logFC, -log10(P.Value), pch=20, col="red"))
with(subset(this_data, adj.P.Val<.05 & abs(logFC)>1 & autophagy != 'x'), points(logFC, -log10(P.Value), pch=20, col="cyan"))
abline(h = -log10(0.005), col = "green3", lty = 2) # adj Pvalue
abline(v = c(-2, 2), col = "blue", lty = 2) #logGC
mtext("adj pval\n = 0.05", side = 2, at = -log10(0.01), cex = 0.6, line = 0.5, las = 1)
mtext(c(paste("-2 fold"), paste("+2 fold")), side = 3, at = c(-2, 2), cex = 0.6, line = 0.2)

mtext(c(paste("DE Autophagy genes: "), paste(dim(subset(this_data, adj.P.Val<0.05 & abs(logFC)>1 & autophagy != 'x'))[1])), side = 3, las =1, at = c(-10,-5), cex = 1)

# Label points with the textxy function from the calibrate plot
#library(calibrate)
with(subset(this_data, adj.P.Val<0.05 & abs(logFC)>1 & autophagy != 'x'), textxy(logFC, -log10(P.Value), labs=autophagy, cex=.6))



