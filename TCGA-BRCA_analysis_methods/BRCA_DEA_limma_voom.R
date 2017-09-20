#! /usr/bin/env Rscript
rm(list=ls(all=TRUE))

library(SummarizedExperiment)
library(plyr)
library(edgeR)
library(limma)
library(ggplot2)
library(RColorBrewer)
library(dplyr)


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


####### SET UP FOR DEA ########### 


# edgeR object
dge <- DGEList(counts=dataSE )

# filtering low expression by cpm       #A CPM value of 1 is equivalent to a log-CPM value of 0

#how many genes have at least 19 samples with 2CPM  expression 
### these numbers are specific to this dataset and were guided by voom plot

cpm <- cpm(dge, log=FALSE)
keep.exprs <- rowSums(cpm > 2) >= 19
dge <- dge[keep.exprs,, keep.lib.sizes=FALSE]
dim(dge) # 

# in case you want to save the genes which passes low count filtering and will be used in DEA
#genes_for_DEA<-rownames(dge)  
#save(genes_for_DEA, file = "genes_for_DEA.rda")

#check autophagy in the current geneset
#autophagy_genes<- as.vector(read.table("autopahagy_genes.txt", as.is = T, header = FALSE))$V1
#shared <- intersect(autophagy_genes,rownames(dge))
#print(paste0("Total number of genes in analysis: ", length(rownames(dge))))
#print(paste0("Autophagy genes: ", length(shared)))


### add sample annotation data to dge object, 
    #relevel vectors to always have normal as referece (where makes sense),
    #replace NAs with character vectors (e.g. "unknown")
    #create Groups

dge<-addSampleData(dge,samples.matrix)



# Scale normalisation: correct for library size
#to match the between-sample distributions of gene counts
#in terms of parameters such as quantiles
y <- calcNormFactors(dge) #default is TMM




##### BUILDING DESIGN MATRICES #############

#pick which one to run

design <- model.matrix(~ 0 + PAM50 + age  + tss + year , data=y$samples) 
design <- model.matrix(~ 0 + stages + age  + tss + year, data=y$samples)
design <- model.matrix(~ 0 + morphology + age + tss + year , data=y$samples) 
design <- model.matrix(~ 0 + Group1 + age + tss + year , data=y$samples) 
design <- model.matrix(~ 0 + Group2 + age + tss + year , data=y$samples) 
design <- model.matrix(~ 0 + Group3 + age + tss + year , data=y$samples) 

#pick which one to run
colnames(design) <- gsub("PAM50", "", colnames(design))
colnames(design) <- gsub("stages", "", colnames(design))
colnames(design) <- gsub("morphology", "", colnames(design))
colnames(design) <- gsub("Group1", "", colnames(design))
colnames(design) <- gsub("Group2", "", colnames(design))
colnames(design) <- gsub("Group3", "", colnames(design))


# run all, always (adjusts colnames in the design matrix)
colnames(design) <- gsub(":age", ":", colnames(design))
colnames(design) <- gsub(" ", "", colnames(design))
colnames(design) <- gsub("-l", "L", colnames(design))
colnames(design) <- gsub("-", "", colnames(design))

colnames(design) 
#write.csv(design, file="testdesign.csv")


############# BUILDING CONTRAST MATRICES ############

# they are very big, so therefore are stored in collapsable ###,
# in the same order as options for design matrices above

### NB only run one!


## separate main effect models
#######################

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
                               ##
                               NormalLikevsLumA = NormalLike- LuminalA,
                               NormalLikevsLumB = NormalLike- LuminalB,
                               NormalLikevsBasal = NormalLike- BasalLike,
                               NormalLikevsHER2 = NormalLike- HER2enriched,
                               NormalLikevsNormal =NormalLike - Normal,
                               
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
                               
                               Stage4vsStage3 = stage4 - stage3,
                               Stage4vsStage2 = stage4 - stage2,
                               Stage4vsStage1 = stage4 - stage1,
                               
                               levels = colnames(design)) 


contr.matrix <- makeContrasts(LobularvsNormal = Lobularcarcinoma - Normal,
                              DuctalvsNormal = Ductalcarcinoma - Normal,
                              DuctLobvsNormal = Ductalandlobularmixed - Normal,
                              MetaplastvsNormal = Metaplasticcarcinoma - Normal,
                              MucinousvsNormal = Mucinousadenocarcinoma - Normal,
                              
                              LobularvsDuctal = Lobularcarcinoma - Ductalcarcinoma,
                              LobularvsDuctLob = Lobularcarcinoma - Ductalandlobularmixed,
                              LobularvsMetaplast = Lobularcarcinoma - Metaplasticcarcinoma,
                              LobularvsMucinous = Lobularcarcinoma - Mucinousadenocarcinoma,
                              
                              DuctalvsLobular = Ductalcarcinoma - Lobularcarcinoma,
                              DuctalvsDuctLob = Ductalcarcinoma - Ductalandlobularmixed,
                              DuctalvsMetaplast = Ductalcarcinoma - Metaplasticcarcinoma,
                              DuctalvsMucinous = Ductalcarcinoma - Mucinousadenocarcinoma,
                              
                              DuctLobvsLobular = Ductalandlobularmixed - Lobularcarcinoma,
                              DuctLobvsDuctal  = Ductalandlobularmixed - Ductalcarcinoma,
                              DuctLobvsMetaplast =Ductalandlobularmixed - Metaplasticcarcinoma,
                              DuctLobvsMucinous  =Ductalandlobularmixed - Mucinousadenocarcinoma,
                              
                              
                              MetaplastvsLobular = Metaplasticcarcinoma -  Lobularcarcinoma,
                              MetaplastvsDuctal = Metaplasticcarcinoma - Ductalcarcinoma,
                              MetaplastvsDuctLob = Metaplasticcarcinoma -  Ductalandlobularmixed,
                              MetaplastvsMucinous = Metaplasticcarcinoma  -   Mucinousadenocarcinoma,
                              
                              MucinousvsLobular = Mucinousadenocarcinoma -  Lobularcarcinoma,
                              MucinousvsDuctal = Mucinousadenocarcinoma -Ductalcarcinoma,
                              MucinousvsDuctLob = Mucinousadenocarcinoma -  Ductalandlobularmixed,
                              MucinousvsMetaplast = Mucinousadenocarcinoma - Metaplasticcarcinoma,
                              
                              levels = colnames(design)) 


################


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
                              stage1.MetaplastvsNormal = Metaplasticcarcinoma.stage1 - Normal.Normal,
                              stage1.MucinousvsNormal = Mucinousadenocarcinoma.stage1 - Normal.Normal,
                              
                              stage1.LobularvsDuctal = Lobularcarcinoma.stage1 - Ductalcarcinoma.stage1,
                              stage1.LobularvsDuctLob = Lobularcarcinoma.stage1 - Ductalandlobularmixed.stage1,
                              stage1.LobularvsMetaplast = Lobularcarcinoma.stage1 - Metaplasticcarcinoma.stage1,
                              stage1.LobularvsMucinous = Lobularcarcinoma.stage1 - Mucinousadenocarcinoma.stage1,
                              
                              stage1.DuctalvsLobular = Ductalcarcinoma.stage1 - Lobularcarcinoma.stage1,
                              stage1.DuctalvsDuctLob = Ductalcarcinoma.stage1 - Ductalandlobularmixed.stage1,
                              stage1.DuctalvsMetaplast = Ductalcarcinoma.stage1 - Metaplasticcarcinoma.stage1,
                              stage1.DuctalvsMucinous = Ductalcarcinoma.stage1 - Mucinousadenocarcinoma.stage1,
                              
                              stage1.DuctLobvsLobular = Ductalandlobularmixed.stage1 - Lobularcarcinoma.stage1,
                              stage1.DuctLobvsDuctal  = Ductalandlobularmixed.stage1 - Ductalcarcinoma.stage1,
                              stage1.DuctLobvsMetaplast =Ductalandlobularmixed.stage1 - Metaplasticcarcinoma.stage1,
                              stage1.DuctLobvsMucinous  =Ductalandlobularmixed.stage1 - Mucinousadenocarcinoma.stage1,
        
                              stage1.MetaplastvsLobular = Metaplasticcarcinoma.stage1 -  Lobularcarcinoma.stage1,
                              stage1.MetaplastvsDuctal = Metaplasticcarcinoma.stage1 - Ductalcarcinoma.stage1,
                              stage1.MetaplastvsDuctLob = Metaplasticcarcinoma.stage1 -  Ductalandlobularmixed.stage1,
                              stage1.MetaplastvsMucinous = Metaplasticcarcinoma.stage1  -   Mucinousadenocarcinoma.stage1,
                              
                              stage1.MucinousvsLobular = Mucinousadenocarcinoma.stage1 -  Lobularcarcinoma.stage1,
                              stage1.MucinousvsDuctal = Mucinousadenocarcinoma.stage1 -Ductalcarcinoma.stage1,
                              stage1.MucinousvsDuctLob = Mucinousadenocarcinoma.stage1 -  Ductalandlobularmixed.stage1,
                              stage1.MucinousvsMetaplast = Mucinousadenocarcinoma.stage1 - Metaplasticcarcinoma.stage1,
                              
                              stage2.LobularvsNormal = Lobularcarcinoma.stage2 - Normal.Normal,
                              stage2.DuctalvsNormal = Ductalcarcinoma.stage2 - Normal.Normal,
                              stage2.DuctLobvsNormal = Ductalandlobularmixed.stage2 - Normal.Normal,
                              stage2.MetaplastvsNormal = Metaplasticcarcinoma.stage2 - Normal.Normal,
                              stage2.MucinousvsNormal = Mucinousadenocarcinoma.stage2 - Normal.Normal,
                              
                              stage2.LobularvsDuctal = Lobularcarcinoma.stage2 - Ductalcarcinoma.stage2,
                              stage2.LobularvsDuctLob = Lobularcarcinoma.stage2 - Ductalandlobularmixed.stage2,
                              stage2.LobularvsMetaplast = Lobularcarcinoma.stage2 - Metaplasticcarcinoma.stage2,
                              stage2.LobularvsMucinous = Lobularcarcinoma.stage2 - Mucinousadenocarcinoma.stage2,
                              
                              stage2.DuctalvsLobular = Ductalcarcinoma.stage2 - Lobularcarcinoma.stage2,
                              stage2.DuctalvsDuctLob = Ductalcarcinoma.stage2 - Ductalandlobularmixed.stage2,
                              stage2.DuctalvsMetaplast = Ductalcarcinoma.stage2 - Metaplasticcarcinoma.stage2,
                              stage2.DuctalvsMucinous = Ductalcarcinoma.stage2 - Mucinousadenocarcinoma.stage2,
                              
                              stage2.DuctLobvsLobular = Ductalandlobularmixed.stage2 - Lobularcarcinoma.stage2,
                              stage2.DuctLobvsDuctal  = Ductalandlobularmixed.stage2 - Ductalcarcinoma.stage2,
                              stage2.DuctLobvsMetaplast =Ductalandlobularmixed.stage2 - Metaplasticcarcinoma.stage2,
                              stage2.DuctLobvsMucinous  =Ductalandlobularmixed.stage2 - Mucinousadenocarcinoma.stage2,
                              
                              stage2.MetaplastvsLobular = Metaplasticcarcinoma.stage2 -  Lobularcarcinoma.stage2,
                              stage2.MetaplastvsDuctal = Metaplasticcarcinoma.stage2 - Ductalcarcinoma.stage2,
                              stage2.MetaplastvsDuctLob = Metaplasticcarcinoma.stage2 -  Ductalandlobularmixed.stage2,
                              stage2.MetaplastvsMucinous = Metaplasticcarcinoma.stage2  -   Mucinousadenocarcinoma.stage2,
                              
                              stage2.MucinousvsLobular = Mucinousadenocarcinoma.stage2 -  Lobularcarcinoma.stage2,
                              stage2.MucinousvsDuctal = Mucinousadenocarcinoma.stage2 -Ductalcarcinoma.stage2,
                              stage2.MucinousvsDuctLob = Mucinousadenocarcinoma.stage2 -  Ductalandlobularmixed.stage2,
                              stage2.MucinousvsMetaplast = Mucinousadenocarcinoma.stage2 - Metaplasticcarcinoma.stage2,
                              
                              
                              stage3.LobularvsNormal = Lobularcarcinoma.stage3 - Normal.Normal,
                              stage3.DuctalvsNormal = Ductalcarcinoma.stage3 - Normal.Normal,
                              stage3.DuctLobvsNormal = Ductalandlobularmixed.stage3 - Normal.Normal,
                              stage3.MetaplastvsNormal = Metaplasticcarcinoma.stage3 - Normal.Normal,
                              stage3.MucinousvsNormal = Mucinousadenocarcinoma.stage3 - Normal.Normal,
                              
                              stage3.LobularvsDuctal = Lobularcarcinoma.stage3 - Ductalcarcinoma.stage3,
                              stage3.LobularvsDuctLob = Lobularcarcinoma.stage3 - Ductalandlobularmixed.stage3,
                              stage3.LobularvsMetaplast = Lobularcarcinoma.stage3 - Metaplasticcarcinoma.stage3,
                              stage3.LobularvsMucinous = Lobularcarcinoma.stage3 - Mucinousadenocarcinoma.stage3,
                              
                              stage3.DuctalvsLobular = Ductalcarcinoma.stage3 - Lobularcarcinoma.stage3,
                              stage3.DuctalvsDuctLob = Ductalcarcinoma.stage3 - Ductalandlobularmixed.stage3,
                              stage3.DuctalvsMetaplast = Ductalcarcinoma.stage3 - Metaplasticcarcinoma.stage3,
                              stage3.DuctalvsMucinous = Ductalcarcinoma.stage3 - Mucinousadenocarcinoma.stage3,
                              
                              stage3.DuctLobvsLobular = Ductalandlobularmixed.stage3 - Lobularcarcinoma.stage3,
                              stage3.DuctLobvsDuctal  = Ductalandlobularmixed.stage3 - Ductalcarcinoma.stage3,
                              stage3.DuctLobvsMetaplast =Ductalandlobularmixed.stage3 - Metaplasticcarcinoma.stage3,
                              stage3.DuctLobvsMucinous  =Ductalandlobularmixed.stage3 - Mucinousadenocarcinoma.stage3,
                              
                              
                              stage3.MetaplastvsLobular = Metaplasticcarcinoma.stage3 -  Lobularcarcinoma.stage3,
                              stage3.MetaplastvsDuctal = Metaplasticcarcinoma.stage3 - Ductalcarcinoma.stage3,
                              stage3.MetaplastvsDuctLob = Metaplasticcarcinoma.stage3 -  Ductalandlobularmixed.stage3,
                              stage3.MetaplastvsMucinous = Metaplasticcarcinoma.stage3  -   Mucinousadenocarcinoma.stage3,
                              
                              stage3.MucinousvsLobular = Mucinousadenocarcinoma.stage3 -  Lobularcarcinoma.stage3,
                              stage3.MucinousvsDuctal = Mucinousadenocarcinoma.stage3 -Ductalcarcinoma.stage3,
                              stage3.MucinousvsDuctLob = Mucinousadenocarcinoma.stage3 -  Ductalandlobularmixed.stage3,
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



######## Differential expression analysis (DEA) ##########

#functions
DEA_limmaVoom <- function(y, design, contr.matrix=NULL) {
  
  #the voom transformation is applied to the normalized and filtered DGEList object (y)
  v <- voom(y, design, plot=F); print ("done Voom")
  
  #After this, the usual limma pipelines for differential expression can be applied
  #fit a separate model to the expression values for each gene
  vfit <- lmFit(v, design); print ("done lmFit")
  
  #check if conrast matrix was specified and apply it
  if (!is.null(contr.matrix)){
    print ("Using contrast matrix... ")
    vfit <- contrasts.fit(vfit, contrasts=contr.matrix)  
  }
  
  #apply empirical Bayes to estimate  gene-wise variability
    efit <- eBayes(vfit); print ("done eBayes")
  
    #optional
  #The model's residual variances are plotted against average expression values
  # plotSA(efit,main="plotSA()")
  
  #testing FDR p-val distribution #always check p-val disrubution to make sure FDR worked
  #hist(as.vector(as.matrix(efit$p.value[,c(1:ncol(efit$p.value))])), col="red", main="P-values before MTC") #NBBBB not sure if this is the right p-vals
  #(mfrow=c(1,1))
  
  return(list(fit=efit, v=v)) 
}  

DEA_filter_save <-function(fit, my_coef, logFC, FDR){
  
  #For a quick look at differential expression levels, the number of significantly up-
  #and down-regulated genes can be summarised in a table.
  dt<-decideTests(fit, p.value=FDR, lfc= logFC, adjust.method="BH")
  
  print (my_coef)
  print (summary(dt)[,c(my_coef)])[2]
 
  
  tt <- topTable(fit, coef=my_coef, adjust='BH', p.value = FDR, lfc = logFC, n=Inf) # 
  
  #this catches situation when there are no DEGs
  if (dim(tt)[1]==0) {
    return(list(tt=c(1,2,3), dt=dt)) #dummy vector
  } else{

  ### below is only used for saving
    
  index.up <- which(tt$logFC >= logFC & tt$adj.P.Val < FDR)
  index.down <- which(tt$logFC <= -logFC & tt$adj.P.Val < FDR)
  direction <- c()
  direction[index.up] <- "up"
  direction[index.down] <- "down"
  tt <- cbind(tt,direction)
  #final <- list(up, down)
  
  up <- data.frame(rownames(tt[tt$direction == "up", ]))
  down <- data.frame(rownames(tt[tt$direction == "down", ]))
  
  ####   uncomment to save lists of up/doem genes for each contrast
  #setwd("DifferentialExpData")
  #write.table(up, paste0(my_coef,"_up.txt"), sep = "\t",col.names = FALSE,row.names = FALSE, quote = FALSE)
  #write.table(down, paste0(my_coef,"_down.txt"),sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)
  #setwd("../")
  return(list(tt=tt, dt=dt))
  }
}

#set parameters!
set_logFC=1
set_FDR=0.05


#using functions

##  1) get fit of design and contrast matrix

DEA_limmaVoom_out<-DEA_limmaVoom(y,design,contr.matrix)
fit<-DEA_limmaVoom_out$fit #fit object
v<-DEA_limmaVoom_out$v #voom object

## 2) first visualise DEGs for each contast (data not stored yet)

                       # here coef does not matter
dt <-DEA_filter_save(fit, my_coef=1, set_logFC , set_FDR)$dt #coef can be anything, the function reports dt for all anyway
#visualise DEGs (up and down, after fold change and FDR filtering)
summary(dt)

##  3) store genes visualised in previoiscommand in a list object
contrast_DE_list<-list()
for (i in colnames(dt)){          

  this_coef<-DEA_filter_save(fit, i, set_logFC , set_FDR)$tt
  if (typeof(this_coef)!="double") {
    contrast_DE_list[[i]]<-this_coef
  }else{
      #skip the ones with dummy c(1,2,3); those contrast will be ignore in further steps!
    }

}
# this is the list object 
names(contrast_DE_list)
#eg
contrast_DE_list$LumAvsLumB


### this is the end of DEA 





####### PREPARING DATA OBJECTS FOR enrichemnt analysis (EA) ####################



# function that extracts gene lists by up/down/both DE
getGeneNames<-function(comparison, direction){
  
  if (direction != 'both'){
    genes<-rownames(comparison[which(comparison$direction==direction),])
  } else{
    genes_up<-rownames(comparison[which(comparison$direction=='up'),])
    genes_down<-rownames(comparison[which(comparison$direction=='down'),])
    genes<-append(genes_up,genes_down)
  }    
  return(as.character(genes))
}


## 1) create data frames with gene number counts in each contast to use in EA 
# (for all genes, autophagy genes, autophagy core genes, and autophagy TF genes)

#create a df to hold all DE genes counts 
genesDEdf<-data.frame(matrix( NA,  nrow = length(names(contrast_DE_list)), ncol = 3))
rownames(genesDEdf)<-names(contrast_DE_list)
colnames(genesDEdf)<-c("up","down","both")
head(genesDEdf)

#create a df to hold DE AUTOPHAGY genes counts 
genesDE_AUTOdf<-data.frame(matrix( NA,  nrow = length(names(contrast_DE_list)), ncol = 3))
rownames(genesDE_AUTOdf)<-names(contrast_DE_list)
colnames(genesDE_AUTOdf)<-c("up","down","both")
head(genesDE_AUTOdf)

#create a df to hold DE AUTOPHAGY CORE genes counts 
genesDE_AUTOCOREdf<-data.frame(matrix( NA,  nrow = length(names(contrast_DE_list)), ncol = 3))
rownames(genesDE_AUTOCOREdf)<-names(contrast_DE_list)
colnames(genesDE_AUTOCOREdf)<-c("up","down","both")
head(genesDE_AUTOCOREdf)

#create a df to hold DE AUTOPHAGY TF genes counts 
genesDE_AUTOTFdf<-data.frame(matrix( NA,  nrow = length(names(contrast_DE_list)), ncol = 3))
rownames(genesDE_AUTOTFdf)<-names(contrast_DE_list)
colnames(genesDE_AUTOTFdf)<-c("up","down","both")
head(genesDE_AUTOTFdf)


## 2) fillin in the created data frames 

# place holders fro total count
all_genesDE <- c() #  count how many genes were DE in all comparisons
autophagy_genesDE <- c()
autophagyCORE_genesDE <- c()
autophagyTF_genesDE <- c()

# fill in df with numbers
for (i in 1:length(rownames(genesDEdf))){
  this_contrast<- rownames(genesDEdf)[i]
  
  # get numbers for all DE
  # for each contast will store number (length function) of up/down/both DEGs genes 
  genesDEdf[i,"up"]<- length(getGeneNames(contrast_DE_list[[this_contrast]],'up'))
  genesDEdf[i,"down"]<- length(getGeneNames(contrast_DE_list[[this_contrast]],'down'))
  genesDEdf[i,"both"]<- length(getGeneNames(contrast_DE_list[[this_contrast]],'both'))
  
  #get numbers for autophagy DE only 
  # same but aonly for autophagy genes (function sharedWithAuto() is defined in functions.R)
  genesDE_AUTOdf[i,"up"]<- length(sharedWithAuto(getGeneNames(contrast_DE_list[[this_contrast]],'up')))
  genesDE_AUTOdf[i,"down"]<- length(sharedWithAuto(getGeneNames(contrast_DE_list[[this_contrast]],'down')))
  genesDE_AUTOdf[i,"both"]<- length(sharedWithAuto(getGeneNames(contrast_DE_list[[this_contrast]],'both')))
  
  #get numbers for autophagy CORE DE only 
  genesDE_AUTOCOREdf[i,"up"]<- length(sharedWithAutoCORE(getGeneNames(contrast_DE_list[[this_contrast]],'up')))
  genesDE_AUTOCOREdf[i,"down"]<- length(sharedWithAutoCORE(getGeneNames(contrast_DE_list[[this_contrast]],'down')))
  genesDE_AUTOCOREdf[i,"both"]<- length(sharedWithAutoCORE(getGeneNames(contrast_DE_list[[this_contrast]],'both')))
  
  #get numbers for autophagy TF DE only 
  genesDE_AUTOTFdf[i,"up"]<- length(sharedWithAutoTF(getGeneNames(contrast_DE_list[[this_contrast]],'up')))
  genesDE_AUTOTFdf[i,"down"]<- length(sharedWithAutoTF(getGeneNames(contrast_DE_list[[this_contrast]],'down')))
  genesDE_AUTOTFdf[i,"both"]<- length(sharedWithAutoTF(getGeneNames(contrast_DE_list[[this_contrast]],'both')))
  
  all_genesDE <- unique(sort(append(all_genesDE, getGeneNames(contrast_DE_list[[this_contrast]],'both'))))
  autophagy_genesDE <- unique(sort(append(autophagy_genesDE, sharedWithAuto(getGeneNames(contrast_DE_list[[this_contrast]],'both')))))
  autophagyCORE_genesDE <- unique(sort(append(autophagyCORE_genesDE, sharedWithAutoCORE(getGeneNames(contrast_DE_list[[this_contrast]],'both')))))
  autophagyTF_genesDE <- unique(sort(append(autophagyTF_genesDE, sharedWithAutoTF(getGeneNames(contrast_DE_list[[this_contrast]],'both')))))

  
}
print("Total number of DEGs:");length(all_genesDE)
print("Autophagy DEGs:");length(autophagy_genesDE)
print("Autophagy core DEGs:");length(autophagyCORE_genesDE)
print("Autophagy TF DEGs:");length(autophagyTF_genesDE)




### 3)  saving data for Enrichment Analysis


setwd("EnrichmentAnalysisData")


# run the block of commands for thr model that you tried

#PAM50
save(genesDEdf, file="PAM50_DE_genes_numbers.rda")
save(genesDE_AUTOdf, file="PAM50_DE_AUTO_genes_numbers.rda")
save(genesDE_AUTOCOREdf, file="PAM50_DE_AUTOCORE_genes_numbers.rda")
save(genesDE_AUTOTFdf, file="PAM50_DE_AUTOTF_genes_numbers.rda")

### STAGES
save(genesDEdf, file="stages_DE_genes_numbers.rda")
save(genesDE_AUTOdf, file="stages_DE_AUTO_genes_numbers.rda")
save(genesDE_AUTOCOREdf, file="stages_DE_AUTOCORE_genes_numbers.rda")
save(genesDE_AUTOTFdf, file="stages_DE_AUTOTF_genes_numbers.rda")

### morphology
save(genesDEdf, file="morph_DE_genes_numbers.rda")
save(genesDE_AUTOdf, file="morph_DE_AUTO_genes_numbers.rda")
save(genesDE_AUTOCOREdf, file="morph_DE_AUTOCORE_genes_numbers.rda")
save(genesDE_AUTOTFdf, file="morph_DE_AUTOTF_genes_numbers.rda")


#Group1
save(genesDEdf, file="Group1_DE_genes_numbers.rda")
save(genesDE_AUTOdf, file="Group1_DE_AUTO_genes_numbers.rda")
save(genesDE_AUTOCOREdf, file="Group1_DE_AUTOCORE_genes_numbers.rda")
save(genesDE_AUTOTFdf, file="Group1_DE_AUTOTF_genes_numbers.rda")


#Group2
save(genesDEdf, file="Group2_DE_genes_numbers.rda")
save(genesDE_AUTOdf, file="Group2_DE_AUTO_genes_numbers.rda")
save(genesDE_AUTOCOREdf, file="Group2_DE_AUTOCORE_genes_numbers.rda")
save(genesDE_AUTOTFdf, file="Group2_DE_AUTOTF_genes_numbers.rda")


#Group3
save(genesDEdf, file="Group3_DE_genes_numbers.rda")
save(genesDE_AUTOdf, file="Group3_DE_AUTO_genes_numbers.rda")
save(genesDE_AUTOCOREdf, file="Group3_DE_AUTOCORE_genes_numbers.rda")
save(genesDE_AUTOTFdf, file="Group3_DE_AUTOTF_genes_numbers.rda")


