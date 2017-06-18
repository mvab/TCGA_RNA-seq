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

# exclude unknown satge and other mrphology samples
#removedsamples<-removeUnknownOther(dataSE, samples.matrix)
#dataSE<-removedsamples$dataSE
#samples.matrix<-removedsamples$samples.matrix
#dim(dataSE)
#dim(samples.matrix)

# exclude mrphology samples Ductal mixed with Others
removeUnknownOther <- function(dataSE, samples.matrix){
  samples.matrix[samples.matrix$tumourTypes=="Ductual mixed with others",]$barcode -> to_remove
  to_remove<- as.character(to_remove)
  
  samples.matrix<-samples.matrix[!samples.matrix$barcode %in% to_remove,]
  dim(samples.matrix)
  dataSE<-dataSE[, !colnames(dataSE) %in% to_remove ]
  dim(dataSE)
  
  return(list(dataSE=dataSE, samples.matrix=samples.matrix))
}

removedsamples<-removeUnknownOther(dataSE, samples.matrix)
dataSE<-removedsamples$dataSE
samples.matrix<-removedsamples$samples.matrix
dim(dataSE)
dim(samples.matrix)


## 

#get gene lenghts for genes in the analysis : NB this 
getGeneLenght_out<-getGeneLenght(dataSE)
gene_lengths<-getGeneLenght_out$gene_lengths
dataSE<-getGeneLenght_out$dataSE
dim(dataSE)
#dim(gene_lengths)

# edgeR object
dge <- DGEList(counts=dataSE, genes = data.frame(SYMBOL= as.character(gene_lengths$gene_name), Length= as.numeric(gene_lengths$gene_length)) ) #here will be adding genes: genes = Ann)
head(dge$genes)

#save(dge, file="dge_with_gene_length_prior_filtering.rda")

#get  entrezIDs - NB ~1649 genes get removed because can't find annotation
library("AnnotationDbi")
library("org.Hs.eg.db")
columns(org.Hs.eg.db)


head(dge$genes)
dge$genes$entrez <- mapIds(org.Hs.eg.db, keys=row.names(dge), column="ENTREZID", keytype="SYMBOL", multiVals="first")
genes_w_entrez<-dge$genes[complete.cases(dge$genes),]$SYMBOL
dge<- dge[rownames(dge) %in% genes_w_entrez,]
dim(dge) 
rownames(dge)<-dge$genes$entrez


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
genes_for_DEA<-rownames(dge)

#save(genes_for_DEA, file = "genes_for_DEA.rda")

#check autophagy in the current geneset
#autophagy_genes<- as.vector(read.table("autopahagy_genes.txt", as.is = T, header = FALSE))$V1
#shared <- intersect(autophagy_genes,rownames(dge))
#print(paste0("Total number of genes in analysis: ", length(rownames(dge))))
#print(paste0("Autophagy genes: ", length(shared)))

m_1_50 # "ASGR1" "RAB9B" "REP15" "KCNE2" unique compared to 3_12
        # "ASGR1" "RAB9B"  unique compared to m_2_25

m_2_25  # "ASGR1"  "ROS1"   "IRS4"   "PRKAG3" unique comapred to 2_19
        #  "RAB9B" "REP15" "KCNE2" unique compared to 3_12

m_3_12 # "RPH3A" "PDX1" uniwue here 
        # "RPH3A"  "PDX1"   "ROS1"   "IRS4"   "PRKAG3" but not in 1_50

m_2_19 # "ASGR1" "RAB9B" "REP15" "KCNE2" unique here



addSampleData<-function(y, samples.matrix) {
  
  # adding samples information
  #y$samples$condition <- as.factor(samples.matrix$condition)
  y$samples$PAM50 <- as.factor(samples.matrix$PAM50)
  y$samples$morphology <- as.factor(samples.matrix$tumourTypes) # from sample lists!
  y$samples$stages <- as.factor(samples.matrix$tumourStages) # from sample lists!
  y$samples$year <- as.factor(samples.matrix$year_diagnosed)
  y$samples$tss <- factor(samples.matrix$tss)
  y$samples$age <- as.factor(samples.matrix$ageGroups)
  

  #fixing NAs
  y$samples$age <- as.character(y$samples$age)
  y$samples[is.na(y$samples$age),]$age<-"UnknownAge"
  y$samples[y$samples$age=="70+",]$age<-"70andUp"
  y$samples[y$samples$age=="< 40",]$age<-"40andDown"
  y$samples$age <- as.factor(y$samples$age)
  
  y$samples$year <- as.character(y$samples$year)
  y$samples[is.na(y$samples$year),]$year<-"UnknownYear"
  y$samples$year <- as.factor(y$samples$year)
  
  
  #stage + PAM50
  y$samples$Group1 <- factor(paste( gsub(" ", "", samples.matrix$PAM50) , samples.matrix$tumourStages,sep="."))
  #morphology +stage
  y$samples$Group2 <- factor(paste( gsub(" ", "", samples.matrix$tumourTypes) , samples.matrix$tumourStages,sep="."))
  #pam50+morphology
  y$samples$Group3 <- factor(paste( gsub(" ", "", samples.matrix$PAM50) , samples.matrix$tumourTypes,sep="."))
  #pam50+age
  y$samples$Group4 <- factor(paste( gsub(" ", "", samples.matrix$PAM50) , y$samples$age,sep="."))
  
  
  
  ## making normal the baselayer
  #y$samples$condition = relevel(y$samples$condition, ref="normal")
  y$samples$PAM50 = relevel(y$samples$PAM50, ref="Normal")
  y$samples$morphology = relevel(y$samples$morphology, ref="Normal")
  y$samples$stages = relevel(y$samples$stages, ref="Normal")
  y$samples$age = relevel(y$samples$age, ref="40andDown") ###########--> should i make one for normal??
  y$samples$Group1 = relevel(y$samples$Group1, ref = "Normal.Normal")
  y$samples$Group2 = relevel(y$samples$Group2, ref = "Normal.Normal")
  y$samples$Group3 = relevel(y$samples$Group3, ref = "Normal.Normal")
  
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

design <- model.matrix(~ 0 + stages + age  + tss + year, data=y$samples)
design <- model.matrix(~ 0 + PAM50 + age  + tss + year , data=y$samples) 
design <- model.matrix(~ 0 + morphology + age + tss + year , data=y$samples) 
design <- model.matrix(~ 0 + Group1 + age + tss + year , data=y$samples) 
design <- model.matrix(~ 0 + Group2 + age + tss + year , data=y$samples) 
design <- model.matrix(~ 0 + Group3 + age + tss + year , data=y$samples) 


is.fullrank(design)
nonEstimable(design)

# check for linearly dependent columns
#rankifremoved <- sapply(1:ncol(design), function (x) qr(design[,-x])$rank)
#which(rankifremoved == max(rankifremoved))
#colnames(design)[2:18]

#design <- model.matrix(~0+Group1, data=y$samples) #setting up model contrasts is more straight forward in the absence of an intercept for

colnames(design) <- gsub("Group1", "", colnames(design))
colnames(design) <- gsub("Group2", "", colnames(design))
colnames(design) <- gsub("Group3", "", colnames(design))
colnames(design) <- gsub("PAM50", "", colnames(design))
colnames(design) <- gsub("morphology", "", colnames(design))
colnames(design) <- gsub("stages", "", colnames(design))


colnames(design) <- gsub(":age", ":", colnames(design))
colnames(design) <- gsub(" ", "", colnames(design))
colnames(design) <- gsub("-l", "L", colnames(design))
colnames(design) <- gsub("-", "", colnames(design))

#colnames(design) <- gsub(":t", "t", colnames(design))

colnames(design) 
#write.csv(design, file="testdesign.csv")

# contr matrices
#######       

# group1

#########
contr.matrix <- makeContrasts(
  ######## BY PAM50  
  LuminalA.stage1vsNorm = LuminalA.stage1 - Normal.Normal,  
  LuminalA.stage2vsNorm = LuminalA.stage2 - Normal.Normal,
  LuminalA.stage3vsNorm = LuminalA.stage3 - Normal.Normal,
  LuminalA.stage4vsNorm = LuminalA.stage4- Normal.Normal,
  
  LuminalA.stage1vsstage2 = LuminalA.stage1 - LuminalA.stage2,
  LuminalA.stage1vsstage3 = LuminalA.stage1 - LuminalA.stage3,
  LuminalA.stage1vsstage4 = LuminalA.stage1 - LuminalA.stage4,
  
  LuminalA.stage2vsstage1 = LuminalA.stage2 - LuminalA.stage1,
  LuminalA.stage2vsstage3 = LuminalA.stage2 - LuminalA.stage3,
  LuminalA.stage2vsstage4 = LuminalA.stage2 - LuminalA.stage4,
  
  LuminalA.stage3vsstage4 = LuminalA.stage3 - LuminalA.stage4,
  LuminalA.stage3vsstage2 = LuminalA.stage3 - LuminalA.stage2,
  LuminalA.stage3vsstage1 = LuminalA.stage3 - LuminalA.stage1,
  
  
  
  LuminalB.stage1vsNorm = LuminalB.stage1 - Normal.Normal,  
  LuminalB.stage2vsNorm = LuminalB.stage2 - Normal.Normal,
  LuminalB.stage3vsNorm = LuminalB.stage3 - Normal.Normal,
  LuminalB.stage4vsNorm = LuminalB.stage4- Normal.Normal,
  
  LuminalB.stage1vsstage2 = LuminalB.stage1 - LuminalB.stage2,
  LuminalB.stage1vsstage3 = LuminalB.stage1 - LuminalB.stage3,
  LuminalB.stage1vsstage4 = LuminalB.stage1 - LuminalB.stage4,
  
  LuminalB.stage2vsstage1 = LuminalB.stage2 - LuminalB.stage1,
  LuminalB.stage2vsstage3 = LuminalB.stage2 - LuminalB.stage3,
  LuminalB.stage2vsstage4 = LuminalB.stage2 - LuminalB.stage4,
  
  LuminalB.stage3vsstage4 = LuminalB.stage3 - LuminalB.stage4,
  LuminalB.stage3vsstage2 = LuminalB.stage3 - LuminalB.stage2,
  LuminalB.stage3vsstage1 = LuminalB.stage3 - LuminalB.stage1,
  
  
  
  BasalLike.stage1vsNorm = BasalLike.stage1 - Normal.Normal,  
  BasalLike.stage2vsNorm = BasalLike.stage2 - Normal.Normal,
  BasalLike.stage3vsNorm = BasalLike.stage3 - Normal.Normal,
  BasalLike.stage4vsNorm = BasalLike.stage4- Normal.Normal,
  
  BasalLike.stage1vsstage2 = BasalLike.stage1 - BasalLike.stage2,
  BasalLike.stage1vsstage3 = BasalLike.stage1 - BasalLike.stage3,
  BasalLike.stage1vsstage4 = BasalLike.stage1 - BasalLike.stage4,
  
  BasalLike.stage2vsstage1 = BasalLike.stage2 - BasalLike.stage1,
  BasalLike.stage2vsstage3 = BasalLike.stage2 - BasalLike.stage3,
  BasalLike.stage2vsstage4 = BasalLike.stage2 - BasalLike.stage4,
  
  BasalLike.stage3vsstage4 = BasalLike.stage3 - BasalLike.stage4,
  BasalLike.stage3vsstage2 = BasalLike.stage3 - BasalLike.stage2,
  BasalLike.stage3vsstage1 = BasalLike.stage3 - BasalLike.stage1,
  
  
  
  HER2enriched.stage1vsNorm = HER2enriched.stage1 - Normal.Normal,  
  HER2enriched.stage2vsNorm = HER2enriched.stage2 - Normal.Normal,
  HER2enriched.stage3vsNorm = HER2enriched.stage3 - Normal.Normal,
  HER2enriched.stage4vsNorm = HER2enriched.stage4- Normal.Normal,
  
  HER2enriched.stage1vsstage2 = HER2enriched.stage1 - HER2enriched.stage2,
  HER2enriched.stage1vsstage3 = HER2enriched.stage1 - HER2enriched.stage3,
  HER2enriched.stage1vsstage4 = HER2enriched.stage1 - HER2enriched.stage4,
  
  HER2enriched.stage2vsstage1 = HER2enriched.stage2 - HER2enriched.stage1,
  HER2enriched.stage2vsstage3 = HER2enriched.stage2 - HER2enriched.stage3,
  HER2enriched.stage2vsstage4 = HER2enriched.stage2 - HER2enriched.stage4,
  
  HER2enriched.stage3vsstage4 = HER2enriched.stage3 - HER2enriched.stage4,
  HER2enriched.stage3vsstage2 = HER2enriched.stage3 - HER2enriched.stage2,
  HER2enriched.stage3vsstage1 = HER2enriched.stage3 - HER2enriched.stage1,
  
  
  NormalLike.stage1vsNorm = NormalLike.stage1 - Normal.Normal,  
  NormalLike.stage2vsNorm = NormalLike.stage2 - Normal.Normal,
  NormalLike.stage3vsNorm = NormalLike.stage3 - Normal.Normal,
  
  NormalLike.stage1vsstage2 = NormalLike.stage1 - NormalLike.stage2,
  NormalLike.stage1vsstage3 = NormalLike.stage1 - NormalLike.stage3,
  
  NormalLike.stage2vsstage1 = NormalLike.stage2 - NormalLike.stage1,
  NormalLike.stage2vsstage3 = NormalLike.stage2 - NormalLike.stage3,
  
  NormalLike.stage3vsstage2 = NormalLike.stage3 - NormalLike.stage2,
  NormalLike.stage3vsstage1 = NormalLike.stage3 - NormalLike.stage1,
  
  # normallike.stage4 does not exist!
  
  #######BY STAGE
  
  stage1.LumAvsLumB = LuminalA.stage1 - LuminalB.stage1,     
  stage1.LumAvsBasal = LuminalA.stage1 - BasalLike.stage1,
  stage1.LumAvsHER2 = LuminalA.stage1 - HER2enriched.stage1,
  stage1.LumAvsNormLike = LuminalA.stage1 - NormalLike.stage1,

  stage1.LumBvsLumA = LuminalB.stage1 - LuminalA.stage1,
  stage1.LumBvsBasal = LuminalB.stage1 - BasalLike.stage1,
  stage1.LumBvsHER2 = LuminalB.stage1 - HER2enriched.stage1,
  stage1.LumBvsNormLike = LuminalB.stage1 - NormalLike.stage1,

  stage1.BasalvsLumA = BasalLike.stage1- LuminalA.stage1,
  stage1.BasalvsLumB = BasalLike.stage1- LuminalB.stage1,
  stage1.BasalvsHER2 = BasalLike.stage1- HER2enriched.stage1,
  stage1.BasalvsNormLike = BasalLike.stage1 - NormalLike.stage1,

  stage1.HER2vsLumA = HER2enriched.stage1- LuminalA.stage1,
  stage1.HER2vsLumB = HER2enriched.stage1- LuminalB.stage1,
  stage1.HER2vsBasal = HER2enriched.stage1- BasalLike.stage1,
  stage1.HER2vsNormLike = HER2enriched.stage1 - NormalLike.stage1,

  stage1.NormalLikevsLumA = NormalLike.stage1- LuminalA.stage1,
  stage1.NormalLikevsLumB = NormalLike.stage1- LuminalB.stage1,
  stage1.NormalLikevsBasal = NormalLike.stage1- BasalLike.stage1,
  stage1.NormalLikevsHER2 = NormalLike.stage1- HER2enriched.stage1,

  
  
  stage2.LumAvsLumB = LuminalA.stage2 - LuminalB.stage2,     
  stage2.LumAvsBasal = LuminalA.stage2 - BasalLike.stage2,
  stage2.LumAvsHER2 = LuminalA.stage2 - HER2enriched.stage2,
  stage2.LumAvsNormLike = LuminalA.stage2 - NormalLike.stage2,

  stage2.LumBvsLumA = LuminalB.stage2 - LuminalA.stage2,
  stage2.LumBvsBasal = LuminalB.stage2 - BasalLike.stage2,
  stage2.LumBvsHER2 = LuminalB.stage2 - HER2enriched.stage2,
  stage2.LumBvsNormLike = LuminalB.stage2 - NormalLike.stage2,

  stage2.BasalvsLumA = BasalLike.stage2- LuminalA.stage2,
  stage2.BasalvsLumB = BasalLike.stage2- LuminalB.stage2,
  stage2.BasalvsHER2 = BasalLike.stage2- HER2enriched.stage2,
  stage2.BasalvsNormLike = BasalLike.stage2 - NormalLike.stage2,

  stage2.HER2vsLumA = HER2enriched.stage2- LuminalA.stage2,
  stage2.HER2vsLumB = HER2enriched.stage2- LuminalB.stage2,
  stage2.HER2vsBasal = HER2enriched.stage2- BasalLike.stage2,
  stage2.HER2vsNormLike = HER2enriched.stage2 - NormalLike.stage2,

  stage2.NormalLikevsLumA = NormalLike.stage2- LuminalA.stage2,
  stage2.NormalLikevsLumB = NormalLike.stage2- LuminalB.stage2,
  stage2.NormalLikevsBasal = NormalLike.stage2- BasalLike.stage2,
  stage2.NormalLikevsHER2 = NormalLike.stage2- HER2enriched.stage2,

  
  
  stage3.LumAvsLumB = LuminalA.stage3 - LuminalB.stage3,     
  stage3.LumAvsBasal = LuminalA.stage3 - BasalLike.stage3,
  stage3.LumAvsHER2 = LuminalA.stage3 - HER2enriched.stage3,
  stage3.LumAvsNormLike = LuminalA.stage3 - NormalLike.stage3,

  stage3.LumBvsLumA = LuminalB.stage3 - LuminalA.stage3,
  stage3.LumBvsBasal = LuminalB.stage3 - BasalLike.stage3,
  stage3.LumBvsHER2 = LuminalB.stage3 - HER2enriched.stage3,
  stage3.LumBvsNormLike = LuminalB.stage3 - NormalLike.stage3,

  stage3.BasalvsLumA = BasalLike.stage3- LuminalA.stage3,
  stage3.BasalvsLumB = BasalLike.stage3- LuminalB.stage3,
  stage3.BasalvsHER2 = BasalLike.stage3- HER2enriched.stage3,
  stage3.BasalvsNormLike = BasalLike.stage3 - NormalLike.stage3,

  stage3.HER2vsLumA = HER2enriched.stage3- LuminalA.stage3,
  stage3.HER2vsLumB = HER2enriched.stage3- LuminalB.stage3,
  stage3.HER2vsBasal = HER2enriched.stage3- BasalLike.stage3,
  stage3.HER2vsNormLike = HER2enriched.stage3 - NormalLike.stage3,

  stage3.NormalLikevsLumA = NormalLike.stage3- LuminalA.stage3,
  stage3.NormalLikevsLumB = NormalLike.stage3- LuminalB.stage3,
  stage3.NormalLikevsBasal = NormalLike.stage3- BasalLike.stage3,
  stage3.NormalLikevsHER2 = NormalLike.stage3- HER2enriched.stage3,

  
  stage4.LumAvsLumB = LuminalA.stage4 - LuminalB.stage4,     
  stage4.LumAvsBasal = LuminalA.stage4 - BasalLike.stage4,
  stage4.LumAvsHER2 = LuminalA.stage4 - HER2enriched.stage4,

  stage4.LumBvsLumA = LuminalB.stage4 - LuminalA.stage4,
  stage4.LumBvsBasal = LuminalB.stage4 - BasalLike.stage4,
  stage4.LumBvsHER2 = LuminalB.stage4 - HER2enriched.stage4,

  stage4.BasalvsLumA = BasalLike.stage4- LuminalA.stage4,
  stage4.BasalvsLumB = BasalLike.stage4- LuminalB.stage4,
  stage4.BasalvsHER2 = BasalLike.stage4- HER2enriched.stage4,

  stage4.HER2vsLumA = HER2enriched.stage4- LuminalA.stage4,
  stage4.HER2vsLumB = HER2enriched.stage4- LuminalB.stage4,
  stage4.HER2vsBasal = HER2enriched.stage4- BasalLike.stage4,

  levels = colnames(design)) 

##########


### group2

###########
contr.matrix <- makeContrasts(stage1.LobularvsNormal = Lobularcarcinoma.stage1 - Normal.Normal,
                              stage1.DuctalvsNormal = Ductalcarcinoma.stage1 - Normal.Normal,
                              stage1.DuctLobvsNormal = Ductalandlobularmixed.stage1 - Normal.Normal,
                              #stage1.DuctOthersvsNormal = Ductualmixedwithothers.stage1 - Normal.Normal,
                              stage1.MetaplastvsNormal = Metaplasticcarcinoma.stage1 - Normal.Normal,
                              stage1.MucinousvsNormal = Mucinousadenocarcinoma.stage1 - Normal.Normal,
                              
                              stage1.LobularvsDuctal = Lobularcarcinoma.stage1 - Ductalcarcinoma.stage1,
                              stage1.LobularvsDuctLob = Lobularcarcinoma.stage1 - Ductalandlobularmixed.stage1,
                              #stage1.LobularvsDuctOther = Lobularcarcinoma.stage1 - Ductualmixedwithothers.stage1,
                              stage1.LobularvsMetaplast = Lobularcarcinoma.stage1 - Metaplasticcarcinoma.stage1,
                              stage1.LobularvsMucinous = Lobularcarcinoma.stage1 - Mucinousadenocarcinoma.stage1,
                              
                              stage1.DuctalvsLobular = Ductalcarcinoma.stage1 - Lobularcarcinoma.stage1,
                              stage1.DuctalvsDuctLob = Ductalcarcinoma.stage1 - Ductalandlobularmixed.stage1,
                              #stage1.DuctalvsDuctOthers = Ductalcarcinoma.stage1 - Ductualmixedwithothers.stage1,
                              stage1.DuctalvsMetaplast = Ductalcarcinoma.stage1 - Metaplasticcarcinoma.stage1,
                              stage1.DuctalvsMucinous = Ductalcarcinoma.stage1 - Mucinousadenocarcinoma.stage1,
                              
                              stage1.DuctLobvsLobular = Ductalandlobularmixed.stage1 - Lobularcarcinoma.stage1,
                              stage1.DuctLobvsDuctal  = Ductalandlobularmixed.stage1 - Ductalcarcinoma.stage1,
                              #stage1.DuctLobvsDuctOther =Ductalandlobularmixed.stage1 - Ductualmixedwithothers.stage1,
                              stage1.DuctLobvsMetaplast =Ductalandlobularmixed.stage1 - Metaplasticcarcinoma.stage1,
                              stage1.DuctLobvsMucinous  =Ductalandlobularmixed.stage1 - Mucinousadenocarcinoma.stage1,
                              
                              #stage1.DuctOthersvsLobular = Ductualmixedwithothers.stage1 - Lobularcarcinoma.stage1,
                              #stage1.DuctOthersvsDuctal = Ductualmixedwithothers.stage1 - Ductalcarcinoma.stage1,
                              #stage1.DuctOthersvsDuctLob = Ductualmixedwithothers.stage1 -   Ductalandlobularmixed.stage1,
                              #stage1.DuctOthersvsMetaplast = Ductualmixedwithothers.stage1 - Metaplasticcarcinoma.stage1,
                              #stage1.DuctOthersvsMucinous  = Ductualmixedwithothers.stage1 -  Mucinousadenocarcinoma.stage1,
                              
                              stage1.MetaplastvsLobular = Metaplasticcarcinoma.stage1 -  Lobularcarcinoma.stage1,
                              stage1.MetaplastvsDuctal = Metaplasticcarcinoma.stage1 - Ductalcarcinoma.stage1,
                              stage1.MetaplastvsDuctLob = Metaplasticcarcinoma.stage1 -  Ductalandlobularmixed.stage1,
                              #stage1.MetaplastvsDuctOther = Metaplasticcarcinoma.stage1 -  Ductualmixedwithothers.stage1,
                              stage1.MetaplastvsMucinous = Metaplasticcarcinoma.stage1  -   Mucinousadenocarcinoma.stage1,
                              
                              stage1.MucinousvsLobular = Mucinousadenocarcinoma.stage1 -  Lobularcarcinoma.stage1,
                              stage1.MucinousvsDuctal = Mucinousadenocarcinoma.stage1 -Ductalcarcinoma.stage1,
                              stage1.MucinousvsDuctLob = Mucinousadenocarcinoma.stage1 -  Ductalandlobularmixed.stage1,
                              #stage1.MucinousvsDuctOther = Mucinousadenocarcinoma.stage1 - Ductualmixedwithothers.stage1,
                              stage1.MucinousvsMetaplast = Mucinousadenocarcinoma.stage1 - Metaplasticcarcinoma.stage1,
                              
                              stage2.LobularvsNormal = Lobularcarcinoma.stage2 - Normal.Normal,
                              stage2.DuctalvsNormal = Ductalcarcinoma.stage2 - Normal.Normal,
                              stage2.DuctLobvsNormal = Ductalandlobularmixed.stage2 - Normal.Normal,
                              #stage2.DuctOthersvsNormal = Ductualmixedwithothers.stage2 - Normal.Normal,
                              stage2.MetaplastvsNormal = Metaplasticcarcinoma.stage2 - Normal.Normal,
                              stage2.MucinousvsNormal = Mucinousadenocarcinoma.stage2 - Normal.Normal,
                              
                              stage2.LobularvsDuctal = Lobularcarcinoma.stage2 - Ductalcarcinoma.stage2,
                              stage2.LobularvsDuctLob = Lobularcarcinoma.stage2 - Ductalandlobularmixed.stage2,
                              #stage2.LobularvsDuctOther = Lobularcarcinoma.stage2 - Ductualmixedwithothers.stage2,
                              stage2.LobularvsMetaplast = Lobularcarcinoma.stage2 - Metaplasticcarcinoma.stage2,
                              stage2.LobularvsMucinous = Lobularcarcinoma.stage2 - Mucinousadenocarcinoma.stage2,
                              
                              stage2.DuctalvsLobular = Ductalcarcinoma.stage2 - Lobularcarcinoma.stage2,
                              stage2.DuctalvsDuctLob = Ductalcarcinoma.stage2 - Ductalandlobularmixed.stage2,
                              #stage2.DuctalvsDuctOthers = Ductalcarcinoma.stage2 - Ductualmixedwithothers.stage2,
                              stage2.DuctalvsMetaplast = Ductalcarcinoma.stage2 - Metaplasticcarcinoma.stage2,
                              stage2.DuctalvsMucinous = Ductalcarcinoma.stage2 - Mucinousadenocarcinoma.stage2,
                              
                              stage2.DuctLobvsLobular = Ductalandlobularmixed.stage2 - Lobularcarcinoma.stage2,
                              stage2.DuctLobvsDuctal  = Ductalandlobularmixed.stage2 - Ductalcarcinoma.stage2,
                              #stage2.DuctLobvsDuctOther =Ductalandlobularmixed.stage2 - Ductualmixedwithothers.stage2,
                              stage2.DuctLobvsMetaplast =Ductalandlobularmixed.stage2 - Metaplasticcarcinoma.stage2,
                              stage2.DuctLobvsMucinous  =Ductalandlobularmixed.stage2 - Mucinousadenocarcinoma.stage2,
                              
                              #stage2.DuctOthersvsLobular = Ductualmixedwithothers.stage2 - Lobularcarcinoma.stage2,
                              #stage2.DuctOthersvsDuctal = Ductualmixedwithothers.stage2 - Ductalcarcinoma.stage2,
                              #stage2.DuctOthersvsDuctLob = Ductualmixedwithothers.stage2 -   Ductalandlobularmixed.stage2,
                              #stage2.DuctOthersvsMetaplast = Ductualmixedwithothers.stage2 - Metaplasticcarcinoma.stage2,
                              #stage2.DuctOthersvsMucinous  = Ductualmixedwithothers.stage2 -  Mucinousadenocarcinoma.stage2,
                              
                              stage2.MetaplastvsLobular = Metaplasticcarcinoma.stage2 -  Lobularcarcinoma.stage2,
                              stage2.MetaplastvsDuctal = Metaplasticcarcinoma.stage2 - Ductalcarcinoma.stage2,
                              stage2.MetaplastvsDuctLob = Metaplasticcarcinoma.stage2 -  Ductalandlobularmixed.stage2,
                              #stage2.MetaplastvsDuctOther = Metaplasticcarcinoma.stage2 -  Ductualmixedwithothers.stage2,
                              stage2.MetaplastvsMucinous = Metaplasticcarcinoma.stage2  -   Mucinousadenocarcinoma.stage2,
                              
                              stage2.MucinousvsLobular = Mucinousadenocarcinoma.stage2 -  Lobularcarcinoma.stage2,
                              stage2.MucinousvsDuctal = Mucinousadenocarcinoma.stage2 -Ductalcarcinoma.stage2,
                              stage2.MucinousvsDuctLob = Mucinousadenocarcinoma.stage2 -  Ductalandlobularmixed.stage2,
                              #stage2.MucinousvsDuctOther = Mucinousadenocarcinoma.stage2 - Ductualmixedwithothers.stage2,
                              stage2.MucinousvsMetaplast = Mucinousadenocarcinoma.stage2 - Metaplasticcarcinoma.stage2,
                              
                              
                              stage3.LobularvsNormal = Lobularcarcinoma.stage3 - Normal.Normal,
                              stage3.DuctalvsNormal = Ductalcarcinoma.stage3 - Normal.Normal,
                              stage3.DuctLobvsNormal = Ductalandlobularmixed.stage3 - Normal.Normal,
                              #stage3.DuctOthersvsNormal = Ductualmixedwithothers.stage3 - Normal.Normal,
                              stage3.MetaplastvsNormal = Metaplasticcarcinoma.stage3 - Normal.Normal,
                              stage3.MucinousvsNormal = Mucinousadenocarcinoma.stage3 - Normal.Normal,
                              
                              stage3.LobularvsDuctal = Lobularcarcinoma.stage3 - Ductalcarcinoma.stage3,
                              stage3.LobularvsDuctLob = Lobularcarcinoma.stage3 - Ductalandlobularmixed.stage3,
                              #stage3.LobularvsDuctOther = Lobularcarcinoma.stage3 - Ductualmixedwithothers.stage3,
                              stage3.LobularvsMetaplast = Lobularcarcinoma.stage3 - Metaplasticcarcinoma.stage3,
                              stage3.LobularvsMucinous = Lobularcarcinoma.stage3 - Mucinousadenocarcinoma.stage3,
                              
                              stage3.DuctalvsLobular = Ductalcarcinoma.stage3 - Lobularcarcinoma.stage3,
                              stage3.DuctalvsDuctLob = Ductalcarcinoma.stage3 - Ductalandlobularmixed.stage3,
                              #stage3.DuctalvsDuctOthers = Ductalcarcinoma.stage3 - Ductualmixedwithothers.stage3,
                              stage3.DuctalvsMetaplast = Ductalcarcinoma.stage3 - Metaplasticcarcinoma.stage3,
                              stage3.DuctalvsMucinous = Ductalcarcinoma.stage3 - Mucinousadenocarcinoma.stage3,
                              
                              stage3.DuctLobvsLobular = Ductalandlobularmixed.stage3 - Lobularcarcinoma.stage3,
                              stage3.DuctLobvsDuctal  = Ductalandlobularmixed.stage3 - Ductalcarcinoma.stage3,
                              #stage3.DuctLobvsDuctOther =Ductalandlobularmixed.stage3 - Ductualmixedwithothers.stage3,
                              stage3.DuctLobvsMetaplast =Ductalandlobularmixed.stage3 - Metaplasticcarcinoma.stage3,
                              stage3.DuctLobvsMucinous  =Ductalandlobularmixed.stage3 - Mucinousadenocarcinoma.stage3,
                              
                              #stage3.DuctOthersvsLobular = Ductualmixedwithothers.stage3 - Lobularcarcinoma.stage3,
                              #stage3.DuctOthersvsDuctal = Ductualmixedwithothers.stage3 - Ductalcarcinoma.stage3,
                              #stage3.DuctOthersvsDuctLob = Ductualmixedwithothers.stage3 -   Ductalandlobularmixed.stage3,
                              #stage3.DuctOthersvsMetaplast = Ductualmixedwithothers.stage3 - Metaplasticcarcinoma.stage3,
                              #stage3.DuctOthersvsMucinous  = Ductualmixedwithothers.stage3 -  Mucinousadenocarcinoma.stage3,
                              
                              stage3.MetaplastvsLobular = Metaplasticcarcinoma.stage3 -  Lobularcarcinoma.stage3,
                              stage3.MetaplastvsDuctal = Metaplasticcarcinoma.stage3 - Ductalcarcinoma.stage3,
                              stage3.MetaplastvsDuctLob = Metaplasticcarcinoma.stage3 -  Ductalandlobularmixed.stage3,
                              #stage3.MetaplastvsDuctOther = Metaplasticcarcinoma.stage3 -  Ductualmixedwithothers.stage3,
                              stage3.MetaplastvsMucinous = Metaplasticcarcinoma.stage3  -   Mucinousadenocarcinoma.stage3,
                              
                              stage3.MucinousvsLobular = Mucinousadenocarcinoma.stage3 -  Lobularcarcinoma.stage3,
                              stage3.MucinousvsDuctal = Mucinousadenocarcinoma.stage3 -Ductalcarcinoma.stage3,
                              stage3.MucinousvsDuctLob = Mucinousadenocarcinoma.stage3 -  Ductalandlobularmixed.stage3,
                              #stage3.MucinousvsDuctOther = Mucinousadenocarcinoma.stage3 - Ductualmixedwithothers.stage3,
                              stage3.MucinousvsMetaplast = Mucinousadenocarcinoma.stage3 - Metaplasticcarcinoma.stage3,
                              
                              
                              
                              stage4.DuctalvsNormal = Ductalcarcinoma.stage4 - Normal.Normal,
                              
                              # stage 4 is only Ductal!
                              
                              
                              Lobularcarcinoma.Stage1vsStage2 = Lobularcarcinoma.stage1 - Lobularcarcinoma.stage2,
                              Lobularcarcinoma.Stage1vsStage3 = Lobularcarcinoma.stage1 - Lobularcarcinoma.stage3,
                              Lobularcarcinoma.Stage2vsStage1 = Lobularcarcinoma.stage2 - Lobularcarcinoma.stage1,
                              Lobularcarcinoma.Stage2vsStage3 = Lobularcarcinoma.stage2 - Lobularcarcinoma.stage3,
                              Lobularcarcinoma.Stage3vsStage2 = Lobularcarcinoma.stage3 - Lobularcarcinoma.stage2,
                              Lobularcarcinoma.Stage3vsStage1 = Lobularcarcinoma.stage3 - Lobularcarcinoma.stage1,
                              
        
                              Ductalcarcinoma.Stage1vsStage2 = Ductalcarcinoma.stage1 - Ductalcarcinoma.stage2,
                              Ductalcarcinoma.Stage1vsStage3 = Ductalcarcinoma.stage1 - Ductalcarcinoma.stage3,
                              Ductalcarcinoma.Stage1vsStage4 = Ductalcarcinoma.stage1 - Ductalcarcinoma.stage4,
                              Ductalcarcinoma.Stage2vsStage1 = Ductalcarcinoma.stage2 - Ductalcarcinoma.stage1,
                              Ductalcarcinoma.Stage2vsStage3 = Ductalcarcinoma.stage2 - Ductalcarcinoma.stage3,
                              Ductalcarcinoma.Stage2vsStage4 = Ductalcarcinoma.stage2 - Ductalcarcinoma.stage4,
                              Ductalcarcinoma.Stage3vsStage4 = Ductalcarcinoma.stage3 - Ductalcarcinoma.stage4,
                              Ductalcarcinoma.Stage3vsStage2 = Ductalcarcinoma.stage3 - Ductalcarcinoma.stage2,
                              Ductalcarcinoma.Stage3vsStage1 = Ductalcarcinoma.stage3 - Ductalcarcinoma.stage1,
                              
                              
               
                              Ductalandlobularmixed.Stage1vsStage2 = Ductalandlobularmixed.stage1 - Ductalandlobularmixed.stage2,
                              Ductalandlobularmixed.Stage1vsStage3 = Ductalandlobularmixed.stage1 - Ductalandlobularmixed.stage3,
                              Ductalandlobularmixed.Stage2vsStage1 = Ductalandlobularmixed.stage2 - Ductalandlobularmixed.stage1,
                              Ductalandlobularmixed.Stage2vsStage3 = Ductalandlobularmixed.stage2 - Ductalandlobularmixed.stage3,
                              Ductalandlobularmixed.Stage3vsStage2 = Ductalandlobularmixed.stage3 - Ductalandlobularmixed.stage2,
                              Ductalandlobularmixed.Stage3vsStage1 = Ductalandlobularmixed.stage3 - Ductalandlobularmixed.stage1,
                              
                              
                              #Ductualmixedwithothers.Stage1vsNorm =   Ductualmixedwithothers.stage1 - Normal.Normal,  
                              #Ductualmixedwithothers.Stage2vsNorm =   Ductualmixedwithothers.stage2 - Normal.Normal,
                              #Ductualmixedwithothers.Stage3vsNorm =   Ductualmixedwithothers.stage3 - Normal.Normal,
                              #Ductualmixedwithothers.Stage1vsStage2 = Ductualmixedwithothers.stage1 - Ductualmixedwithothers.stage2,
                              #Ductualmixedwithothers.Stage1vsStage3 = Ductualmixedwithothers.stage1 - Ductualmixedwithothers.stage3,
                              #Ductualmixedwithothers.Stage2vsStage1 = Ductualmixedwithothers.stage2 - Ductualmixedwithothers.stage1,
                              #Ductualmixedwithothers.Stage2vsStage3 = Ductualmixedwithothers.stage2 - Ductualmixedwithothers.stage3,
                              #Ductualmixedwithothers.Stage3vsStage2 = Ductualmixedwithothers.stage3 - Ductualmixedwithothers.stage2,
                              #Ductualmixedwithothers.Stage3vsStage1 = Ductualmixedwithothers.stage3 - Ductualmixedwithothers.stage1,
                              
                              
                              Metaplasticcarcinoma.Stage1vsStage2 = Metaplasticcarcinoma.stage1 - Metaplasticcarcinoma.stage2,
                              Metaplasticcarcinoma.Stage1vsStage3 = Metaplasticcarcinoma.stage1 - Metaplasticcarcinoma.stage3,
                              Metaplasticcarcinoma.Stage2vsStage1 = Metaplasticcarcinoma.stage2 - Metaplasticcarcinoma.stage1,
                              Metaplasticcarcinoma.Stage2vsStage3 = Metaplasticcarcinoma.stage2 - Metaplasticcarcinoma.stage3,
                              Metaplasticcarcinoma.Stage3vsStage2 = Metaplasticcarcinoma.stage3 - Metaplasticcarcinoma.stage2,
                              Metaplasticcarcinoma.Stage3vsStage1 = Metaplasticcarcinoma.stage3 - Metaplasticcarcinoma.stage1,
                              
                             
                              Mucinousadenocarcinoma.Stage1vsStage2 = Mucinousadenocarcinoma.stage1 - Mucinousadenocarcinoma.stage2,
                              Mucinousadenocarcinoma.Stage1vsStage3 = Mucinousadenocarcinoma.stage1 - Mucinousadenocarcinoma.stage3,
                              Mucinousadenocarcinoma.Stage2vsStage1 = Mucinousadenocarcinoma.stage2 - Mucinousadenocarcinoma.stage1,
                              Mucinousadenocarcinoma.Stage2vsStage3 = Mucinousadenocarcinoma.stage2 - Mucinousadenocarcinoma.stage3,
                              Mucinousadenocarcinoma.Stage3vsStage2 = Mucinousadenocarcinoma.stage3 - Mucinousadenocarcinoma.stage2,
                              Mucinousadenocarcinoma.Stage3vsStage1 = Mucinousadenocarcinoma.stage3 - Mucinousadenocarcinoma.stage1,
                              
                              levels = colnames(design)) 

############


#group3
############
contr.matrix <- makeContrasts( Lobularcarcinoma.LumAvsLumB = LuminalA.Lobularcarcinoma - LuminalB.Lobularcarcinoma,     #unique comparisons 15
                               Lobularcarcinoma.LumAvsBasal = LuminalA.Lobularcarcinoma - BasalLike.Lobularcarcinoma,
                               Lobularcarcinoma.LumAvsHER2 = LuminalA.Lobularcarcinoma - HER2enriched.Lobularcarcinoma,
                               Lobularcarcinoma.LumAvsNormLike = LuminalA.Lobularcarcinoma - NormalLike.Lobularcarcinoma,

                               Lobularcarcinoma.LumBvsLumA = LuminalB.Lobularcarcinoma - LuminalA.Lobularcarcinoma,
                               Lobularcarcinoma.LumBvsBasal = LuminalB.Lobularcarcinoma - BasalLike.Lobularcarcinoma,
                               Lobularcarcinoma.LumBvsHER2 = LuminalB.Lobularcarcinoma - HER2enriched.Lobularcarcinoma,
                               Lobularcarcinoma.LumBvsNormLike = LuminalB.Lobularcarcinoma - NormalLike.Lobularcarcinoma,

                               Lobularcarcinoma.BasalvsLumA = BasalLike.Lobularcarcinoma- LuminalA.Lobularcarcinoma,
                               Lobularcarcinoma.BasalvsLumB = BasalLike.Lobularcarcinoma- LuminalB.Lobularcarcinoma,
                               Lobularcarcinoma.BasalvsHER2 = BasalLike.Lobularcarcinoma- HER2enriched.Lobularcarcinoma,
                               Lobularcarcinoma.BasalvsNormLike = BasalLike.Lobularcarcinoma - NormalLike.Lobularcarcinoma,

                               Lobularcarcinoma.HER2vsLumA = HER2enriched.Lobularcarcinoma- LuminalA.Lobularcarcinoma,
                               Lobularcarcinoma.HER2vsLumB = HER2enriched.Lobularcarcinoma- LuminalB.Lobularcarcinoma,
                               Lobularcarcinoma.HER2vsBasal = HER2enriched.Lobularcarcinoma- BasalLike.Lobularcarcinoma,
                               Lobularcarcinoma.HER2vsNormLike = HER2enriched.Lobularcarcinoma - NormalLike.Lobularcarcinoma,

                               Lobularcarcinoma.NormalLikevsLumA = NormalLike.Lobularcarcinoma- LuminalA.Lobularcarcinoma,
                               Lobularcarcinoma.NormalLikevsLumB = NormalLike.Lobularcarcinoma- LuminalB.Lobularcarcinoma,
                               Lobularcarcinoma.NormalLikevsBasal = NormalLike.Lobularcarcinoma- BasalLike.Lobularcarcinoma,
                               Lobularcarcinoma.NormalLikevsHER2 = NormalLike.Lobularcarcinoma- HER2enriched.Lobularcarcinoma,

                               
                               Ductalcarcinoma.LumAvsLumB = LuminalA.Ductalcarcinoma - LuminalB.Ductalcarcinoma,     #unique comparisons 15
                               Ductalcarcinoma.LumAvsBasal = LuminalA.Ductalcarcinoma - BasalLike.Ductalcarcinoma,
                               Ductalcarcinoma.LumAvsHER2 = LuminalA.Ductalcarcinoma - HER2enriched.Ductalcarcinoma,
                               Ductalcarcinoma.LumAvsNormLike = LuminalA.Ductalcarcinoma - NormalLike.Ductalcarcinoma,

                               Ductalcarcinoma.LumBvsLumA = LuminalB.Ductalcarcinoma - LuminalA.Ductalcarcinoma,
                               Ductalcarcinoma.LumBvsBasal = LuminalB.Ductalcarcinoma - BasalLike.Ductalcarcinoma,
                               Ductalcarcinoma.LumBvsHER2 = LuminalB.Ductalcarcinoma - HER2enriched.Ductalcarcinoma,
                               Ductalcarcinoma.LumBvsNormLike = LuminalB.Ductalcarcinoma - NormalLike.Ductalcarcinoma,

                               Ductalcarcinoma.BasalvsLumA = BasalLike.Ductalcarcinoma- LuminalA.Ductalcarcinoma,
                               Ductalcarcinoma.BasalvsLumB = BasalLike.Ductalcarcinoma- LuminalB.Ductalcarcinoma,
                               Ductalcarcinoma.BasalvsHER2 = BasalLike.Ductalcarcinoma- HER2enriched.Ductalcarcinoma,
                               Ductalcarcinoma.BasalvsNormLike = BasalLike.Ductalcarcinoma - NormalLike.Ductalcarcinoma,

                               Ductalcarcinoma.HER2vsLumA = HER2enriched.Ductalcarcinoma- LuminalA.Ductalcarcinoma,
                               Ductalcarcinoma.HER2vsLumB = HER2enriched.Ductalcarcinoma- LuminalB.Ductalcarcinoma,
                               Ductalcarcinoma.HER2vsBasal = HER2enriched.Ductalcarcinoma- BasalLike.Ductalcarcinoma,
                               Ductalcarcinoma.HER2vsNormLike = HER2enriched.Ductalcarcinoma - NormalLike.Ductalcarcinoma,

                               Ductalcarcinoma.NormalLikevsLumA = NormalLike.Ductalcarcinoma- LuminalA.Ductalcarcinoma,
                               Ductalcarcinoma.NormalLikevsLumB = NormalLike.Ductalcarcinoma- LuminalB.Ductalcarcinoma,
                               Ductalcarcinoma.NormalLikevsBasal = NormalLike.Ductalcarcinoma- BasalLike.Ductalcarcinoma,
                               Ductalcarcinoma.NormalLikevsHER2 = NormalLike.Ductalcarcinoma- HER2enriched.Ductalcarcinoma,

                               Ductalandlobularmixed.LumAvsLumB = LuminalA.Ductalandlobularmixed - LuminalB.Ductalandlobularmixed,     #unique comparisons 15
                               Ductalandlobularmixed.LumAvsBasal = LuminalA.Ductalandlobularmixed - BasalLike.Ductalandlobularmixed,

                               Ductalandlobularmixed.LumBvsLumA = LuminalB.Ductalandlobularmixed - LuminalA.Ductalandlobularmixed,
                               Ductalandlobularmixed.LumBvsBasal = LuminalB.Ductalandlobularmixed - BasalLike.Ductalandlobularmixed,

                               Ductalandlobularmixed.BasalvsLumA = BasalLike.Ductalandlobularmixed- LuminalA.Ductalandlobularmixed,
                               Ductalandlobularmixed.BasalvsLumB = BasalLike.Ductalandlobularmixed- LuminalB.Ductalandlobularmixed,

                               
                              #Ductualmixedwithothers.LumAvsLumB = LuminalA.Ductualmixedwithothers - LuminalB.Ductualmixedwithothers,     #unique comparisons 15
                              #Ductualmixedwithothers.LumAvsBasal = LuminalA.Ductualmixedwithothers - BasalLike.Ductualmixedwithothers,
                              #Ductualmixedwithothers.LumAvsHER2 = LuminalA.Ductualmixedwithothers - HER2enriched.Ductualmixedwithothers,
                              #Ductualmixedwithothers.LumAvsNormal = LuminalA.Ductualmixedwithothers - Normal.Normal,
                              #
                              #Ductualmixedwithothers.LumBvsLumA = LuminalB.Ductualmixedwithothers - LuminalA.Ductualmixedwithothers,
                              #Ductualmixedwithothers.LumBvsBasal = LuminalB.Ductualmixedwithothers - BasalLike.Ductualmixedwithothers,
                              #Ductualmixedwithothers.LumBvsHER2 = LuminalB.Ductualmixedwithothers - HER2enriched.Ductualmixedwithothers,
                              #Ductualmixedwithothers.LumBvsNormal = LuminalB.Ductualmixedwithothers - Normal.Normal,
                              #
                              #Ductualmixedwithothers.BasalvsLumA = BasalLike.Ductualmixedwithothers- LuminalA.Ductualmixedwithothers,
                              #Ductualmixedwithothers.BasalvsLumB = BasalLike.Ductualmixedwithothers- LuminalB.Ductualmixedwithothers,
                              #Ductualmixedwithothers.BasalvsHER2 = BasalLike.Ductualmixedwithothers- HER2enriched.Ductualmixedwithothers,
                              #Ductualmixedwithothers.BasalvsNormal = BasalLike.Ductualmixedwithothers - Normal.Normal,
                              #
                              #Ductualmixedwithothers.HER2vsLumA = HER2enriched.Ductualmixedwithothers- LuminalA.Ductualmixedwithothers,
                              #Ductualmixedwithothers.HER2vsLumB = HER2enriched.Ductualmixedwithothers- LuminalB.Ductualmixedwithothers,
                              #Ductualmixedwithothers.HER2vsBasal = HER2enriched.Ductualmixedwithothers- BasalLike.Ductualmixedwithothers,
                              #Ductualmixedwithothers.HER2vsNormal = HER2enriched.Ductualmixedwithothers - Normal.Normal,
                              #
                               
                               # metaplastic exists only in basal!

                               # mucinous only in LumA and LumB
                               
                               Mucinousadenocarcinoma.LumAvsLumB = LuminalA.Mucinousadenocarcinoma - LuminalB.Mucinousadenocarcinoma,     #unique comparisons 15

                               Mucinousadenocarcinoma.LumBvsLumA = LuminalB.Mucinousadenocarcinoma - LuminalA.Mucinousadenocarcinoma,

                               LuminalA.LobularvsNormal = LuminalA.Lobularcarcinoma - Normal.Normal,
                               LuminalA.DuctalvsNormal = LuminalA.Ductalcarcinoma - Normal.Normal,
                               LuminalA.DuctLobvsNormal = LuminalA.Ductalandlobularmixed - Normal.Normal,
                               #LuminalA.DuctOthersvsNormal = LuminalA.Ductualmixedwithothers - Normal.Normal,
                               LuminalA.MucinousvsNormal = LuminalA.Mucinousadenocarcinoma - Normal.Normal,
                               
                               LuminalA.LobularvsDuctal = LuminalA.Lobularcarcinoma - LuminalA.Ductalcarcinoma,
                               LuminalA.LobularvsDuctLob = LuminalA.Lobularcarcinoma - LuminalA.Ductalandlobularmixed,
                               #LuminalA.LobularvsDuctOther = LuminalA.Lobularcarcinoma - LuminalA.Ductualmixedwithothers,
                               LuminalA.LobularvsMucinous = LuminalA.Lobularcarcinoma - LuminalA.Mucinousadenocarcinoma,
                               
                               LuminalA.DuctalvsLobular = LuminalA.Ductalcarcinoma - LuminalA.Lobularcarcinoma,
                               LuminalA.DuctalvsDuctLob = LuminalA.Ductalcarcinoma - LuminalA.Ductalandlobularmixed,
                               #LuminalA.DuctalvsDuctOthers = LuminalA.Ductalcarcinoma - LuminalA.Ductualmixedwithothers,
                               LuminalA.DuctalvsMucinous = LuminalA.Ductalcarcinoma - LuminalA.Mucinousadenocarcinoma,
                               
                               LuminalA.DuctLobvsLobular = LuminalA.Ductalandlobularmixed - LuminalA.Lobularcarcinoma,
                               LuminalA.DuctLobvsDuctal  = LuminalA.Ductalandlobularmixed - LuminalA.Ductalcarcinoma,
                               #LuminalA.DuctLobvsDuctOther =LuminalA.Ductalandlobularmixed - LuminalA.Ductualmixedwithothers,
                               LuminalA.DuctLobvsMucinous  =LuminalA.Ductalandlobularmixed - LuminalA.Mucinousadenocarcinoma,
                               
                               #LuminalA.DuctOthersvsLobular = LuminalA.Ductualmixedwithothers - LuminalA.Lobularcarcinoma,
                               #LuminalA.DuctOthersvsDuctal = LuminalA.Ductualmixedwithothers - LuminalA.Ductalcarcinoma,
                               #LuminalA.DuctOthersvsDuctLob = LuminalA.Ductualmixedwithothers -  - LuminalA.Ductalandlobularmixed,
                               #LuminalA.DuctOthersvsMucinous  = LuminalA.Ductualmixedwithothers -  LuminalA.Mucinousadenocarcinoma,
                               
                               
                               LuminalA.MucinousvsLobular = LuminalA.Mucinousadenocarcinoma -  LuminalA.Lobularcarcinoma,
                               LuminalA.MucinousvsDuctal = LuminalA.Mucinousadenocarcinoma -LuminalA.Ductalcarcinoma,
                               LuminalA.MucinousvsDuctLob = LuminalA.Mucinousadenocarcinoma -  LuminalA.Ductalandlobularmixed,
                               #LuminalA.MucinousvsDuctOther = LuminalA.Mucinousadenocarcinoma - LuminalA.Ductualmixedwithothers,
                               
                               
                               LuminalB.LobularvsNormal = LuminalB.Lobularcarcinoma - Normal.Normal,
                               LuminalB.DuctalvsNormal = LuminalB.Ductalcarcinoma - Normal.Normal,
                               LuminalB.DuctLobvsNormal = LuminalB.Ductalandlobularmixed - Normal.Normal,
                               #LuminalB.DuctOthersvsNormal = LuminalB.Ductualmixedwithothers - Normal.Normal,
                               LuminalB.MucinousvsNormal = LuminalB.Mucinousadenocarcinoma - Normal.Normal,
                               
                               LuminalB.LobularvsDuctal = LuminalB.Lobularcarcinoma - LuminalB.Ductalcarcinoma,
                               LuminalB.LobularvsDuctLob = LuminalB.Lobularcarcinoma - LuminalB.Ductalandlobularmixed,
                               #LuminalB.LobularvsDuctOther = LuminalB.Lobularcarcinoma - LuminalB.Ductualmixedwithothers,
                               LuminalB.LobularvsMucinous = LuminalB.Lobularcarcinoma - LuminalB.Mucinousadenocarcinoma,
                               
                               LuminalB.DuctalvsLobular = LuminalB.Ductalcarcinoma - LuminalB.Lobularcarcinoma,
                               LuminalB.DuctalvsDuctLob = LuminalB.Ductalcarcinoma - LuminalB.Ductalandlobularmixed,
                               #LuminalB.DuctalvsDuctOthers = LuminalB.Ductalcarcinoma - LuminalB.Ductualmixedwithothers,
                               LuminalB.DuctalvsMucinous = LuminalB.Ductalcarcinoma - LuminalB.Mucinousadenocarcinoma,
                               
                               LuminalB.DuctLobvsLobular = LuminalB.Ductalandlobularmixed - LuminalB.Lobularcarcinoma,
                               LuminalB.DuctLobvsDuctal  = LuminalB.Ductalandlobularmixed - LuminalB.Ductalcarcinoma,
                               #LuminalB.DuctLobvsDuctOther =LuminalB.Ductalandlobularmixed - LuminalB.Ductualmixedwithothers,
                               LuminalB.DuctLobvsMucinous  =LuminalB.Ductalandlobularmixed - LuminalB.Mucinousadenocarcinoma,
                               
                               #LuminalB.DuctOthersvsLobular = LuminalB.Ductualmixedwithothers - LuminalB.Lobularcarcinoma,
                               #LuminalB.DuctOthersvsDuctal = LuminalB.Ductualmixedwithothers - LuminalB.Ductalcarcinoma,
                               #LuminalB.DuctOthersvsDuctLob = LuminalB.Ductualmixedwithothers -  - LuminalB.Ductalandlobularmixed,
                               #LuminalB.DuctOthersvsMucinous  = LuminalB.Ductualmixedwithothers -  LuminalB.Mucinousadenocarcinoma,
                               
                               
                               LuminalB.MucinousvsLobular = LuminalB.Mucinousadenocarcinoma -  LuminalB.Lobularcarcinoma,
                               LuminalB.MucinousvsDuctal = LuminalB.Mucinousadenocarcinoma -LuminalB.Ductalcarcinoma,
                               LuminalB.MucinousvsDuctLob = LuminalB.Mucinousadenocarcinoma -  LuminalB.Ductalandlobularmixed,
                               #LuminalB.MucinousvsDuctOther = LuminalB.Mucinousadenocarcinoma - LuminalB.Ductualmixedwithothers,
                               
                               
                               BasalLike.LobularvsNormal = BasalLike.Lobularcarcinoma - Normal.Normal,
                               BasalLike.DuctalvsNormal = BasalLike.Ductalcarcinoma - Normal.Normal,
                               BasalLike.DuctLobvsNormal = BasalLike.Ductalandlobularmixed - Normal.Normal,
                               #BasalLike.DuctOthersvsNormal = BasalLike.Ductualmixedwithothers - Normal.Normal,
                               BasalLike.MetaplastvsNormal = BasalLike.Metaplasticcarcinoma - Normal.Normal,
                               
                               BasalLike.LobularvsDuctal = BasalLike.Lobularcarcinoma - BasalLike.Ductalcarcinoma,
                               BasalLike.LobularvsDuctLob = BasalLike.Lobularcarcinoma - BasalLike.Ductalandlobularmixed,
                               #BasalLike.LobularvsDuctOther = BasalLike.Lobularcarcinoma - BasalLike.Ductualmixedwithothers,
                               BasalLike.LobularvsMetaplast = BasalLike.Lobularcarcinoma - BasalLike.Metaplasticcarcinoma,
                               
                               BasalLike.DuctalvsLobular = BasalLike.Ductalcarcinoma - BasalLike.Lobularcarcinoma,
                               BasalLike.DuctalvsDuctLob = BasalLike.Ductalcarcinoma - BasalLike.Ductalandlobularmixed,
                               #BasalLike.DuctalvsDuctOthers = BasalLike.Ductalcarcinoma - BasalLike.Ductualmixedwithothers,
                               BasalLike.DuctalvsMetaplast = BasalLike.Ductalcarcinoma - BasalLike.Metaplasticcarcinoma,
                               
                               BasalLike.DuctLobvsLobular = BasalLike.Ductalandlobularmixed - BasalLike.Lobularcarcinoma,
                               BasalLike.DuctLobvsDuctal  = BasalLike.Ductalandlobularmixed - BasalLike.Ductalcarcinoma,
                               #BasalLike.DuctLobvsDuctOther =BasalLike.Ductalandlobularmixed - BasalLike.Ductualmixedwithothers,
                               BasalLike.DuctLobvsMetaplast =BasalLike.Ductalandlobularmixed - BasalLike.Metaplasticcarcinoma,
                               
                               #BasalLike.DuctOthersvsLobular = BasalLike.Ductualmixedwithothers - BasalLike.Lobularcarcinoma,
                               #BasalLike.DuctOthersvsDuctal = BasalLike.Ductualmixedwithothers - BasalLike.Ductalcarcinoma,
                               #BasalLike.DuctOthersvsDuctLob = BasalLike.Ductualmixedwithothers -   BasalLike.Ductalandlobularmixed,
                               #BasalLike.DuctOthersvsMetaplast = BasalLike.Ductualmixedwithothers - BasalLike.Metaplasticcarcinoma,
                               
                               BasalLike.MetaplastvsLobular = BasalLike.Metaplasticcarcinoma -  BasalLike.Lobularcarcinoma,
                               BasalLike.MetaplastvsDuctal = BasalLike.Metaplasticcarcinoma - BasalLike.Ductalcarcinoma,
                               BasalLike.MetaplastvsDuctLob = BasalLike.Metaplasticcarcinoma -  BasalLike.Ductalandlobularmixed,
                               #BasalLike.MetaplastvsDuctOther = BasalLike.Metaplasticcarcinoma -  BasalLike.Ductualmixedwithothers,
                               
                               
                               HER2enriched.LobularvsNormal = HER2enriched.Lobularcarcinoma - Normal.Normal,
                               HER2enriched.DuctalvsNormal = HER2enriched.Ductalcarcinoma - Normal.Normal,
                               #HER2enriched.DuctOthersvsNormal = HER2enriched.Ductualmixedwithothers - Normal.Normal,
                               
                               HER2enriched.LobularvsDuctal = HER2enriched.Lobularcarcinoma - HER2enriched.Ductalcarcinoma,
                               #HER2enriched.LobularvsDuctOther = HER2enriched.Lobularcarcinoma - HER2enriched.Ductualmixedwithothers,
                               
                               HER2enriched.DuctalvsLobular = HER2enriched.Ductalcarcinoma - HER2enriched.Lobularcarcinoma,
                               #HER2enriched.DuctalvsDuctOthers = HER2enriched.Ductalcarcinoma - HER2enriched.Ductualmixedwithothers,
                               
                               #HER2enriched.DuctOthersvsLobular = HER2enriched.Ductualmixedwithothers - HER2enriched.Lobularcarcinoma,
                               #HER2enriched.DuctOthersvsDuctal = HER2enriched.Ductualmixedwithothers - HER2enriched.Ductalcarcinoma,
                               
                               NormalLike.LobularvsNormal = NormalLike.Lobularcarcinoma - Normal.Normal,
                               NormalLike.DuctalvsNormal = NormalLike.Ductalcarcinoma - Normal.Normal,
                               NormalLike.LobularvsDuctal = NormalLike.Lobularcarcinoma - NormalLike.Ductalcarcinoma,
                               NormalLike.DuctalvsLobular = NormalLike.Ductalcarcinoma - NormalLike.Lobularcarcinoma,
                               
            
                               
                               levels = colnames(design)) 


######

#group 4
###############
contr.matrix <- makeContrasts(x40andDown.LumAvsLumB = LuminalA.40andDown - LuminalB.40andDown,     #unique comparisons 15
                              x40andDown.LumAvsBasal = LuminalA.40andDown - BasalLike.40andDown,
                              x40andDown.LumAvsHER2 = LuminalA.40andDown - HER2enriched.40andDown,
                              x40andDown.LumAvsNormLike = LuminalA.40andDown - NormalLike.40andDown,
                              x40andDown.LumAvsNormal = LuminalA.40andDown - Normal.40andDown,
                              
                              x40andDown.LumBvsLumA = LuminalB.40andDown - LuminalA.40andDown,
                              x40andDown.LumBvsBasal = LuminalB.40andDown - BasalLike.40andDown,
                              x40andDown.LumBvsHER2 = LuminalB.40andDown - HER2enriched.40andDown,
                              x40andDown.LumBvsNormLike = LuminalB.40andDown - NormalLike.40andDown,
                              x40andDown.LumBvsNormal = LuminalB.40andDown - Normal.40andDown,
                              
                              x40andDown.BasalvsLumA = BasalLike.40andDown- LuminalA.40andDown,
                              x40andDown.BasalvsLumB = BasalLike.40andDown- LuminalB.40andDown,
                              x40andDown.BasalvsHER2 = BasalLike.40andDown- HER2enriched.40andDown,
                              x40andDown.BasalvsNormLike = BasalLike.40andDown - NormalLike.40andDown,
                              x40andDown.BasalvsNormal = BasalLike.40andDown - Normal.40andDown,
                              
                              x40andDown.HER2vsLumA = HER2enriched.40andDown- LuminalA.40andDown,
                              x40andDown.HER2vsLumB = HER2enriched.40andDown- LuminalB.40andDown,
                              x40andDown.HER2vsBasal = HER2enriched.40andDown- BasalLike.40andDown,
                              x40andDown.HER2vsNormLike = HER2enriched.40andDown - NormalLike.40andDown,
                              x40andDown.HER2vsNormal = HER2enriched.40andDown - Normal.40andDown,
                              
                              x40andDown.NormalLikevsLumA = NormalLike.40andDown- LuminalA.40andDown,
                              x40andDown.NormalLikevsLumB = NormalLike.40andDown- LuminalB.40andDown,
                              x40andDown.NormalLikevsBasal = NormalLike.40andDown- BasalLike.40andDown,
                              x40andDown.NormalLikevsHER2 = NormalLike.40andDown- HER2enriched.40andDown,
                              x40andDown.NormalLikevsNormal =NormalLike.40andDown- Normal.40andDown,
                              
                              
                              
                              x4049.LumAvsLumB = LuminalA.4049 - LuminalB.4049,     #unique compar
                              x4049.LumAvsBasal = LuminalA.4049 - BasalLike.4049,
                              x4049.LumAvsHER2 = LuminalA.4049 - HER2enriched.4049,
                              x4049.LumAvsNormLike = LuminalA.4049 - NormalLike.4049,
                              x4049.LumAvsNormal = LuminalA.4049 - Normal.4049,
                              x4049.LumBvsLumA = LuminalB.4049 - LuminalA.4049,
                              x4049.LumBvsBasal = LuminalB.4049 - BasalLike.4049,
                              x4049.LumBvsHER2 = LuminalB.4049 - HER2enriched.4049,
                              x4049.LumBvsNormLike = LuminalB.4049 - NormalLike.4049,
                              x4049.LumBvsNormal = LuminalB.4049 - Normal.4049,
                              x4049.BasalvsLumA = BasalLike.4049- LuminalA.4049,
                              x4049.BasalvsLumB = BasalLike.4049- LuminalB.4049,
                              x4049.BasalvsHER2 = BasalLike.4049- HER2enriched.4049,
                              x4049.BasalvsNormLike = BasalLike.4049 - NormalLike.4049,
                              x4049.BasalvsNormal = BasalLike.4049 - Normal.4049,
                              x4049.HER2vsLumA = HER2enriched.4049- LuminalA.4049,
                              x4049.HER2vsLumB = HER2enriched.4049- LuminalB.4049,
                              x4049.HER2vsBasal = HER2enriched.4049- BasalLike.4049,
                              x4049.HER2vsNormLike = HER2enriched.4049 - NormalLike.4049,
                              x4049.HER2vsNormal = HER2enriched.4049 - Normal.4049,
                              x4049.NormalLikevsLumA = NormalLike.4049- LuminalA.4049,
                              x4049.NormalLikevsLumB = NormalLike.4049- LuminalB.4049,
                              x4049.NormalLikevsBasal = NormalLike.4049- BasalLike.4049,
                              x4049.NormalLikevsHER2 = NormalLike.4049- HER2enriched.4049,
                              x4049.NormalLikevsNormal =NormalLike.4049- Normal.4049,
                              
                              
                              
                              
                              x5059.LumAvsLumB = LuminalA.5059 - LuminalB.5059,     #unique compar
                              x5059.LumAvsBasal = LuminalA.5059 - BasalLike.5059,
                              x5059.LumAvsHER2 = LuminalA.5059 - HER2enriched.5059,
                              x5059.LumAvsNormLike = LuminalA.5059 - NormalLike.5059,
                              x5059.LumAvsNormal = LuminalA.5059 - Normal.5059,
                              x5059.LumBvsLumA = LuminalB.5059 - LuminalA.5059,
                              x5059.LumBvsBasal = LuminalB.5059 - BasalLike.5059,
                              x5059.LumBvsHER2 = LuminalB.5059 - HER2enriched.5059,
                              x5059.LumBvsNormLike = LuminalB.5059 - NormalLike.5059,
                              x5059.LumBvsNormal = LuminalB.5059 - Normal.5059,
                              x5059.BasalvsLumA = BasalLike.5059- LuminalA.5059,
                              x5059.BasalvsLumB = BasalLike.5059- LuminalB.5059,
                              x5059.BasalvsHER2 = BasalLike.5059- HER2enriched.5059,
                              x5059.BasalvsNormLike = BasalLike.5059 - NormalLike.5059,
                              x5059.BasalvsNormal = BasalLike.5059 - Normal.5059,
                              x5059.HER2vsLumA = HER2enriched.5059- LuminalA.5059,
                              x5059.HER2vsLumB = HER2enriched.5059- LuminalB.5059,
                              x5059.HER2vsBasal = HER2enriched.5059- BasalLike.5059,
                              x5059.HER2vsNormLike = HER2enriched.5059 - NormalLike.5059,
                              x5059.HER2vsNormal = HER2enriched.5059 - Normal.5059,
                              x5059.NormalLikevsLumA = NormalLike.5059- LuminalA.5059,
                              x5059.NormalLikevsLumB = NormalLike.5059- LuminalB.5059,
                              x5059.NormalLikevsBasal = NormalLike.5059- BasalLike.5059,
                              x5059.NormalLikevsHER2 = NormalLike.5059- HER2enriched.5059,
                              x5059.NormalLikevsNormal =NormalLike.5059- Normal.5059,
                              
                              
                              
                              x6069.LumAvsLumB = LuminalA.6069 - LuminalB.6069,     #unique compar
                              x6069.LumAvsBasal = LuminalA.6069 - BasalLike.6069,
                              x6069.LumAvsHER2 = LuminalA.6069 - HER2enriched.6069,
                              x6069.LumAvsNormLike = LuminalA.6069 - NormalLike.6069,
                              x6069.LumAvsNormal = LuminalA.6069 - Normal.6069,
                              x6069.LumBvsLumA = LuminalB.6069 - LuminalA.6069,
                              x6069.LumBvsBasal = LuminalB.6069 - BasalLike.6069,
                              x6069.LumBvsHER2 = LuminalB.6069 - HER2enriched.6069,
                              x6069.LumBvsNormLike = LuminalB.6069 - NormalLike.6069,
                              x6069.LumBvsNormal = LuminalB.6069 - Normal.6069,
                              x6069.BasalvsLumA = BasalLike.6069- LuminalA.6069,
                              x6069.BasalvsLumB = BasalLike.6069- LuminalB.6069,
                              x6069.BasalvsHER2 = BasalLike.6069- HER2enriched.6069,
                              x6069.BasalvsNormLike = BasalLike.6069 - NormalLike.6069,
                              x6069.BasalvsNormal = BasalLike.6069 - Normal.6069,
                              x6069.HER2vsLumA = HER2enriched.6069- LuminalA.6069,
                              x6069.HER2vsLumB = HER2enriched.6069- LuminalB.6069,
                              x6069.HER2vsBasal = HER2enriched.6069- BasalLike.6069,
                              x6069.HER2vsNormLike = HER2enriched.6069 - NormalLike.6069,
                              x6069.HER2vsNormal = HER2enriched.6069 - Normal.6069,
                              x6069.NormalLikevsLumA = NormalLike.6069- LuminalA.6069,
                              x6069.NormalLikevsLumB = NormalLike.6069- LuminalB.6069,
                              x6069.NormalLikevsBasal = NormalLike.6069- BasalLike.6069,
                              x6069.NormalLikevsHER2 = NormalLike.6069- HER2enriched.6069,
                              x6069.NormalLikevsNormal =NormalLike.6069- Normal.6069,
                              
                              
                              
                              x70andUp.LumAvsLumB = LuminalA.70andUp - LuminalB.70andUp,     #unique compar
                              x70andUp.LumAvsBasal = LuminalA.70andUp - BasalLike.70andUp,
                              x70andUp.LumAvsHER2 = LuminalA.70andUp - HER2enriched.70andUp,
                              x70andUp.LumAvsNormLike = LuminalA.70andUp - NormalLike.70andUp,
                              x70andUp.LumAvsNormal = LuminalA.70andUp - Normal.70andUp,
                              x70andUp.LumBvsLumA = LuminalB.70andUp - LuminalA.70andUp,
                              x70andUp.LumBvsBasal = LuminalB.70andUp - BasalLike.70andUp,
                              x70andUp.LumBvsHER2 = LuminalB.70andUp - HER2enriched.70andUp,
                              x70andUp.LumBvsNormLike = LuminalB.70andUp - NormalLike.70andUp,
                              x70andUp.LumBvsNormal = LuminalB.70andUp - Normal.70andUp,
                              x70andUp.BasalvsLumA = BasalLike.70andUp- LuminalA.70andUp,
                              x70andUp.BasalvsLumB = BasalLike.70andUp- LuminalB.70andUp,
                              x70andUp.BasalvsHER2 = BasalLike.70andUp- HER2enriched.70andUp,
                              x70andUp.BasalvsNormLike = BasalLike.70andUp - NormalLike.70andUp,
                              x70andUp.BasalvsNormal = BasalLike.70andUp - Normal.70andUp,
                              x70andUp.HER2vsLumA = HER2enriched.70andUp- LuminalA.70andUp,
                              x70andUp.HER2vsLumB = HER2enriched.70andUp- LuminalB.70andUp,
                              x70andUp.HER2vsBasal = HER2enriched.70andUp- BasalLike.70andUp,
                              x70andUp.HER2vsNormLike = HER2enriched.70andUp - NormalLike.70andUp,
                              x70andUp.HER2vsNormal = HER2enriched.70andUp - Normal.70andUp,
                              x70andUp.NormalLikevsLumA = NormalLike.70andUp- LuminalA.70andUp,
                              x70andUp.NormalLikevsLumB = NormalLike.70andUp- LuminalB.70andUp,
                              x70andUp.NormalLikevsBasal = NormalLike.70andUp- BasalLike.70andUp,
                              x70andUp.NormalLikevsHER2 = NormalLike.70andUp- HER2enriched.70andUp,
                              x70andUp.NormalLikevsNormal =NormalLike.70andUp- Normal.70andUp,
                              
                              levels = colnames(design)) 
###############

##other contrasts
#######################
contr.matrix <- makeContrasts(LobularvsNormal = Lobularcarcinoma - Normal,
                              DuctalvsNormal = Ductalcarcinoma - Normal,
                              DuctLobvsNormal = Ductalandlobularmixed - Normal,
                              #DuctOthersvsNormal = Ductualmixedwithothers - Normal,
                              MetaplastvsNormal = Metaplasticcarcinoma - Normal,
                              MucinousvsNormal = Mucinousadenocarcinoma - Normal,
                              
                              LobularvsDuctal = Lobularcarcinoma - Ductalcarcinoma,
                              LobularvsDuctLob = Lobularcarcinoma - Ductalandlobularmixed,
                              #LobularvsDuctOther = Lobularcarcinoma - Ductualmixedwithothers,
                              LobularvsMetaplast = Lobularcarcinoma - Metaplasticcarcinoma,
                              LobularvsMucinous = Lobularcarcinoma - Mucinousadenocarcinoma,
                              
                              DuctalvsLobular = Ductalcarcinoma - Lobularcarcinoma,
                              DuctalvsDuctLob = Ductalcarcinoma - Ductalandlobularmixed,
                              #DuctalvsDuctOthers = Ductalcarcinoma - Ductualmixedwithothers,
                              DuctalvsMetaplast = Ductalcarcinoma - Metaplasticcarcinoma,
                              DuctalvsMucinous = Ductalcarcinoma - Mucinousadenocarcinoma,
                              
                              DuctLobvsLobular = Ductalandlobularmixed - Lobularcarcinoma,
                              DuctLobvsDuctal  = Ductalandlobularmixed - Ductalcarcinoma,
                              #DuctLobvsDuctOther =Ductalandlobularmixed - Ductualmixedwithothers,
                              DuctLobvsMetaplast =Ductalandlobularmixed - Metaplasticcarcinoma,
                              DuctLobvsMucinous  =Ductalandlobularmixed - Mucinousadenocarcinoma,
                              
                              #DuctOthersvsLobular = Ductualmixedwithothers - Lobularcarcinoma,
                             # DuctOthersvsDuctal = Ductualmixedwithothers - Ductalcarcinoma,
                              #DuctOthersvsDuctLob = Ductualmixedwithothers -  - Ductalandlobularmixed,
                              #DuctOthersvsMetaplast = Ductualmixedwithothers - - Metaplasticcarcinoma,
                              #DuctOthersvsMucinous  = Ductualmixedwithothers -  Mucinousadenocarcinoma,
                              
                              MetaplastvsLobular = Metaplasticcarcinoma -  Lobularcarcinoma,
                              MetaplastvsDuctal = Metaplasticcarcinoma - Ductalcarcinoma,
                              MetaplastvsDuctLob = Metaplasticcarcinoma -  Ductalandlobularmixed,
                              #MetaplastvsDuctOther = Metaplasticcarcinoma -  Ductualmixedwithothers,
                              MetaplastvsMucinous = Metaplasticcarcinoma  -   Mucinousadenocarcinoma,
                              
                              MucinousvsLobular = Mucinousadenocarcinoma -  Lobularcarcinoma,
                              MucinousvsDuctal = Mucinousadenocarcinoma -Ductalcarcinoma,
                              MucinousvsDuctLob = Mucinousadenocarcinoma -  Ductalandlobularmixed,
                              #MucinousvsDuctOther = Mucinousadenocarcinoma - Ductualmixedwithothers,
                              MucinousvsMetaplast = Mucinousadenocarcinoma - Metaplasticcarcinoma,
                              
                              levels = colnames(design)) 
  

contr.matrix <- makeContrasts( Stage1vsNorm = stage1 - Normal,  
                               Stage2vsNorm = stage2 - Normal,
                               Stage3vsNorm = stage3 - Normal,
                               Stage4vsNorm = stage4 - Normal,
                               
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
  
  setwd("~/Bioinformatics MSc UCPH/0_MasterThesis/TCGAbiolinks/CBL_scripts/data/DEA/Group1")
  #write.table(up, paste0(my_coef,"_up.txt"), sep = "\t",col.names = FALSE,row.names = FALSE, quote = FALSE)
  #write.table(down, paste0(my_coef,"_down.txt"),sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)
  
  return(list(tt=tt, dt=dt))
}
DEA_limmaVoom_out<-DEA_limmaVoom(y,design,contr.matrix)
fit<-DEA_limmaVoom_out$fit
#v<-DEA_limmaVoom_out$v



################### GO #################################
###### testing goana!  ######

####NB!!! entrez is has to be instead of gene symbol! 
# So : fix names like for rpkm, convert to entrez - rerun!


go <- goana(fit, coef="Stage4vsNorm",species = "Hs")

#The top set of most enriched GO terms can be viewed with the topGO function.
topGO(go, n=10)

topGO(go, ontology="BP")    # can specify GO domain of interest!


#The row names of the output are the universal identifiers of the GO terms, 
#with one term per row. The Term column gives the names of the GO terms.
#These terms cover three domains - biological process (BP), cellular component (CC) 
#and molecular function (MF), as shown in the Ont column. The N column represents 
#the total number of genes that are annotated with each GO term. The Up and Down 
#columns represent the number of differentially expressed genes that overlap with
#the genes in the GO term. The P.Up and P.Down columns contain the p-values for 
#over-representation of the GO term across the set of up- and down-regulated genes,
# respectively. The output table is sorted by the minimum of P.Up and P.Down by default.

#An additional refinement is to supply goana with the gene lengths using the
#covariate argument. (suppley as simple vector!)

library(BiasedUrn)
go_length <- goana(fit, coef="Stage4vsNorm",species = "Hs", covariate=dge$genes$Length)
topGO(go_length, ontology="BP") 

##  camera + MSIgDB
load("MSIgDB/human_c5_v5p2.rdata")
names(Hs.c5)

#The gene identifiers are Entrez Gene ID, which is the same as the rownames 
#of our voom object. We need to map the Entrez gene ids between the list of
#gene sets and our voom object. We can do this using the ids2indices function.

c5.ind <- ids2indices(Hs.c5, rownames(v))

#CAMERA takes as input the voom object v, the indexed list of gene sets c5.ind,
#the design matrix, the contrast being tested, as well as some other arguments.

#By default, CAMERA can estimate the correlation for each gene set separately. 
#However, in practise, it works well to set a small inter-gene correlation of 
#about 0.05 using the inter.gene.cor argument.


colnames(design[,colSums(design)!=0]) -> cols_to_keep
design_copy<- design[, colnames(design) %in% cols_to_keep]
dim(design_copy) 

dim(contr.matrix_copy)
contr.matrix_copy<-contr.matrix[rownames(contr.matrix) %in% cols_to_keep ,]
###### TRY remove zero cols before model estimation??? why are they even thre? is there data??
gst.camera <- camera(v,index=c5.ind,design=design,contrast = contr.matrix[,c("Stage4vsNorm")],inter.gene.cor=0.05)

#CAMERA outputs a dataframe of the resulting statistics, with each row denoting
#a different gene set. The output is ordered by p-value so that the most significant
#should be at the top. Let's look at the top 5 gene sets:
  
gst.camera[1:5,]

#The total number of significant gene sets at 5% FDR is

table(gst.camera$FDR < 0.05)
gst.camera[gst.camera$FDR < 0.05,]
#You can write out the camera results to a csv file to open in excel.

#write.csv(gst.camera,file="gst_BPregVsLac.csv")


stop()
###############


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

set_logFC=2
set_FDR=0.01

## decideTests        # here coef does not matter
dt <-DEA_MTC_save(fit, 6, set_logFC , set_FDR)$dt #coef can be anything, report dt for all anyway
summary(dt)

colnames(summary(dt))


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

setwd("~/Bioinformatics MSc UCPH/0_MasterThesis/TCGAbiolinks/CBL_scripts/data")



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


#dt_summaries<-list()
#dt_summaries<-get(load("groupsDEgenes.rda"))
#dt_summaries$group3FC2 <-dt
#names(dt_summaries)

save(dt_summaries, file="groupsDEgenes.rda")

dt_groups<-get(load("groupsDEgenes.rda"))
names(dt_groups)
summary(dt_groups$group1FC2)
               


save(all_genesDE, file="pam50_DE_genes.rda")

# with 3>12 : 7559, 396, 34, 51
# with 2>19 7622,399,34,50

# when group1 is in the model, for PAM50: 6818 / 412/ 52 /48

save(all_genesDE, file="group1_all_genesDE.rda")
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


#Group3
save(genesDEdf, file="Group3_DE_genes_numbers3.rda")
save(genesDE_AUTOdf, file="Group3_DE_AUTO_genes_numbers3.rda")
save(genesDE_AUTOCOREdf, file="Group3_DE_AUTOCORE_genes_numbers3.rda")
save(genesDE_AUTOTFdf, file="Group3_DE_AUTOTF_genes_numbers3.rda")


#GroupX 
save(genesDEdf, file="GroupX_DE_genes_numbers3.rda")
save(genesDE_AUTOdf, file="GroupX_DE_AUTO_genes_numbers3.rda")
save(genesDE_AUTOCOREdf, file="GroupX_DE_AUTOCORE_genes_numbers3.rda")
save(genesDE_AUTOTFdf, file="GroupX_DE_AUTOTF_genes_numbers3.rda")


#Group4
save(genesDEdf, file="Group4_DE_genes_numbers3.rda")
save(genesDE_AUTOdf, file="Group4_DE_AUTO_genes_numbers3.rda")
save(genesDE_AUTOCOREdf, file="Group4_DE_AUTOCORE_genes_numbers3.rda")
save(genesDE_AUTOTFdf, file="Group4_DE_AUTOTF_genes_numbers3.rda")



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
par(mfrow=c(1,4))


plotmyMD<- function (fit, dt, i){
  summa<-as.matrix(summary(dt))
  plotMD(fit, column=i, status=dt[,i],
       main=paste(colnames(fit)[i], "\n","up:", summa[1,i], "down:", summa[3,i]), 
       ylim=c(-10,10))
}

par(mfrow=c(1,4))
plotmyMD(fit, dt, 1)
plotmyMD(fit, dt, 2)
plotmyMD(fit, dt, 3)
plotmyMD(fit, dt, 4)

par(mfrow=c(1,3))
plotmyMD(fit, dt, 5)
plotmyMD(fit, dt, 6)
plotmyMD(fit, dt, 7)

plotmyMD(fit, dt, 8)
plotmyMD(fit, dt, 9)
plotmyMD(fit, dt, 10)

plotmyMD(fit, dt, 11)
plotmyMD(fit, dt, 12)
plotmyMD(fit, dt, 13)

par(mfrow=c(2,3))
plotmyMD(fit, dt, 1)
plotmyMD(fit, dt, 2)
plotmyMD(fit, dt, 3)
plotmyMD(fit, dt, 4)
plotmyMD(fit, dt, 5)
plotmyMD(fit, dt, 6)


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



