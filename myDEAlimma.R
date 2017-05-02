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
dim(samples.matrix)
sampleTokeepSE <- PAMNPout$samplesToKeep
dataSE<- dataSE[,colnames(dataSE) %in% sampleTokeepSE]
dim(dataSE) 

## 

addSampleData<-function(y, samples.matrix){
  
  # adding samples information
  y$samples$condition <- as.factor(samples.matrix$condition)
  y$samples$PAM50 <- as.factor(samples.matrix$PAM50)
  y$samples$morphology <- as.factor(samples.matrix$tumourTypes) # from sample lists!
  y$samples$stages <- as.factor(samples.matrix$tumourStages) # from sample lists!
  y$samples$year <- as.factor(samples.matrix$year_diagnosed)
  y$samples$tss <- as.factor(samples.matrix$tss)
  y$samples$age <- as.factor(samples.matrix$ageGroups)
  
  y$samples$Group1 <- factor(paste( gsub(" ", "", samples.matrix$PAM50) , samples.matrix$tumourStages,sep="."))
  
  
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

  
  return (y)
}


#add genes annotation ----- FIGURE THIS OUT LATER
#library("BSgenome.Hsapiens.UCSC.hg19")
#geneid <- rownames(y) 
#genes <- select(BSgenome.Hsapiens.UCSC.hg19, keys=geneid, columns=c("SYMBOL", "TXCHROM"), keytype="ENTREZID")
#dim(genes)
#genes <- genes[!duplicated(genes$ENTREZID),]


#editing genes names in dataSE: these 3 genes have entrzid instaead of geneame
#2    155060 LOC155060
#3      8225    GTPBP6
#4     90288   EFCAB12

head(rownames(dataSE))
rownames(dataSE)[2]<-"LOC155060"
rownames(dataSE)[3]<-"GTPBP6"
rownames(dataSE)[4]<-"EFCAB12"

# edgeR object
dge <- DGEList(counts=dataSE) #here will be adding genes: genes = Ann)

#### adding genes ######

dge$genes <- data.frame(SYMBOL=rownames(dge))

#library(mygene)
#q <-queryMany(dge$genes$SYMBOL, scopes="symbol", fields=c("entrezgene"), species="human")

#che k if got all the genes
#'%nin%' <- Negate('%in%')
#dge$genes$SYMBOL[dge$genes$SYMBOL%nin% q$query]

#qtest<-queryMany(c("CSDAP1", "MGC2752","NCRNA00171"), scopes="symbol", fields=c("entrezgene","ensembl.gene"), species="human")
# 3 genes taht dont have entrez
#1 ENSG00000261614  117.6303      CSDAP1 ENSG00000261614
#2 ENSG00000268784  117.6677     MGC2752 ENSG00000268784
#3 ENSG00000229653  117.6849  NCRNA00171 ENSG00000229653
#q[q$query=='CSDAP1',]$`_id`<-NA
#q[q$query=='MGC2752',]$`_id`<-NA
#q[q$query=='NCRNA00171',]$`_id`<-NA

#remove duplicated with ENSembl id
#q<-q[!(startsWith(q$`_id`, "ENSG0000") & !is.na(q$`_id`)),]
#remove just duplicates
#q<- subset(q, !duplicated(query))
#update
#dge$genes$ENTREZID <- q$entrezgene

head(dge$genes)
##################


dge<-addSampleData(dge,samples.matrix)

# counts per million
log.cpm <- cpm(dge, log=TRUE)
cpm <- cpm(dge, log=FALSE)

# filtering low expression

#how many genes have at least 50 samples with 0 expression
table(rowSums(dge$counts==0)==25) 
#Keep genes with total counts more than 50.
keep.exprs <- rowSums(cpm > 1) >= 50        #A CPM value of 1 is equivalent to a log-CPM value of 0.
dge <- dge[keep.exprs,, keep.lib.sizes=FALSE]
dim(dge) #removed 1647

#check autophagy in the current geneset
#autophagy_genes<- as.vector(read.table("autopahagy_genes.txt", as.is = T, header = FALSE))$V1
#shared <- intersect(autophagy_genes,rownames(dge))
#print(paste0("Total number of genes in analysis: ", length(rownames(dge))))
#print(paste0("Autophagy genes: ", length(shared)))


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

#try later
#library(Glimma)
#glMDSPlot(log.cpm, labels=paste(PAM50, condition, sep="_"), groups=y$samples[,c("PAM50", "condition")], launch=TRUE)





design <- model.matrix(~0+condition, data=y$samples) #setting up model contrasts is more straight forward in the absence of an intercept for
design <- model.matrix(~0+PAM50+condition, data=y$samples) #setting up model contrasts is more straight forward in the absence of an intercept for
design <- model.matrix(~PAM50 + PAM50:age  +morphology , data=y$samples) #nested interaction # makes all possible pairings (alt to manual contarsts)

design <- model.matrix(~0 + PAM50 +condition + age + stages  + morphology + tss + year , data=y$samples) #nested interaction # makes all possible pairings (alt to manual contarsts)
design <- model.matrix(~0 + stages  +condition + PAM50 , data=y$samples) #nested interaction # makes all possible pairings (alt to manual contarsts)


design <- model.matrix(~PAM50+condition, data=y$samples) 
design <- model.matrix(~morphology+condition, data=y$samples) #setting up model contrasts is more straight forward in the absence of an intercept for
design <- model.matrix(~0+stages+condition, data=y$samples) #setting up model contrasts is more straight forward in the absence of an intercept for
design <- model.matrix(~0+Group1, data=y$samples) #setting up model contrasts is more straight forward in the absence of an intercept for

colnames(design) <- gsub("Group1", "", colnames(design))
colnames(design) <- gsub("PAM50", "", colnames(design))
colnames(design) <- gsub("morphology", "", colnames(design))
colnames(design) <- gsub("stages", "", colnames(design))
colnames(design) <- gsub(":age", ":", colnames(design))


colnames(design) <- gsub(" ", "", colnames(design))
colnames(design) <- gsub("-l", "L", colnames(design))
colnames(design) <- gsub("-", "", colnames(design))


colnames(design)


# contr matrices
#######       

contr.matrix <- makeContrasts(Basal_st1vs4= BasalLike.stage1-BasalLike.stage4,
                              LumA_st1vs4= LuminalA.stage1-LuminalA.stage4,
                              allBasal = (BasalLike.stage1+BasalLike.stage2+BasalLike.stage3+BasalLike.stage4+BasalLike.unknown)/5 - Normal.Normal,
                              allLumA = (LuminalA.stage1 + LuminalA.stage2 + LuminalA.stage3 +LuminalA.stage4 +LuminalA.unknown)/5 - Normal.Normal,
                              Basal_stage1 = BasalLike.stage1 - Normal.Normal,
                              Basal_stage4 = BasalLike.stage4 - Normal.Normal,
                              LumA_stage1 = LuminalA.stage1 - Normal.Normal,
                              LumA_stage4 = LuminalA.stage4 - Normal.Normal,
                              BasalVSLumA_stage1 = (BasalLike.stage1 - Normal.Normal) - (LuminalA.stage1 - Normal.Normal),
                              BasalVSLumA_stage4 = (BasalLike.stage4 - Normal.Normal) - (LuminalA.stage4 - Normal.Normal),
                              average_stage4= (BasalLike.stage4+LuminalA.stage4)/2 -(BasalLike.stage1+LuminalA.stage1)/2, #average stage4 effect regardless of subtype
                              LumA_Basal_diff=(LuminalA.stage1+LuminalA.stage4)/2 +(BasalLike.stage1+BasalLike.stage4)/2,
                              
                              levels = colnames(design)) 
                              ####HERE!
                              
contr.matrix <- makeContrasts( Stage1vsNorm = stage1 - Normal,  #unique comparisons 15
                               Stage2vsNorm = stage2 - Normal,
                               Stage3vsNorm = stage3 - Normal,
                               Stage4vsNorm = stage4 - Normal,
                               
                               Stage1vsStage2 = stage1 - stage2,
                               Stage1vsStage3 = stage1 - stage3,
                               Stage1vsStage4 = stage1 - stage4,
                               Stage2vsStage3 = stage2 - stage3,
                               Stage2vsStage4 = stage2 - stage4,
                               Stage3vsStage4 = stage3 - stage4,
                               
                               levels = colnames(design)) 


contr.matrix <- makeContrasts( LumAvsLumB = LuminalA - LuminalB,     #unique comparisons 15
                               LumAvsBasal = LuminalA - BasalLike,
                               LumAvsHER2 = LuminalA - HER2enriched,
                               LumAvsNormLike = LuminalA - NormalLike,
                               LumAvsNormal = LuminalA - Normal,
                               #
                               LumBvsBasal = LuminalB - BasalLike,
                               LumBvsHER2 = LuminalB - HER2enriched,
                               LumBvsNormLike = LuminalB - NormalLike,
                               LumBvsNormal = LuminalB - Normal,
                               #
                               BasalvsHER2 = BasalLike- HER2enriched,
                               BasalvsNormLike = BasalLike - NormalLike,
                               BasalvsNormal = BasalLike - Normal,
                               #
                               HER2vsNormLike = HER2enriched - NormalLike,
                               HER2vsNormal = HER2enriched - Normal,
                               #
                               NormalLikevsNormal = NormalLike - Normal,

                               levels = colnames(design)) 
#####



DEA_limmaVoom <- function(y, design, contr.matrix=NULL) {
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
  efit <- eBayes(vfit); print ("done eBayes")

  #The model's residual variances are plotted against average expression values
 # plotSA(efit,main="plotSA()")

  #testing FDR p-val distribution #always check p-val disrubution to make sure FDR worked
  #hist(as.vector(as.matrix(efit$p.value[,c(1:ncol(efit$p.value))])), col="red", main="P-values before MTC") #NBBBB not sure if this is the right p-vals
  par(mfrow=c(1,1))
 
  return(list(fit=efit, v=v)) #change if choose stricter option
}  

DEA_limmaVoom_out<-DEA_limmaVoom(y,design,contr.matrix)
fit<-DEA_limmaVoom_out$fit



DEA_MTC_save <-function(fit, my_coef, logFC, FDR){
  
  #For a quick look at differential expression levels, the number of significantly up-
  #and down-regulated genes can be summarised in a table.
  dt<-decideTests(fit, p.value=FDR, lfc= logFC, adjust.method="BH")
  
  # for stricter filtering on significance
  
      # The treat method can be used to calculate p-values from
      # empirical Bayes moderated t-statistics with a minimum log-FC requirement
      #when testing requires genes to have a log-FC that is significantly 
      #greater than 1 (equivalent to a 2-fold difference between cell types on the original scale)
  #tfit <- treat(vfit, lfc=1) 
  #dt <- decideTests(tfit)
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
  #up <- data.frame(rownames(tt[tt$direction == "up", ]))
  #down <- data.frame(rownames(tt[tt$direction == "down", ]))
  #write.table(up, "up_filename.txt", sep = "\t",col.names = FALSE,row.names = FALSE, quote = FALSE)
  #write.table(down, "down_filename.txt", sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)
  
  return(list(tt=tt, dt=dt))
}



#set parameters!
set_logFC=1
set_FDR=0.05

colnames(design)
#it thinks that <40 is the base factor for age group --- is that OK?

LumA_age4049 <-DEA_MTC_save(fit, "LuminalA:4049", set_logFC , set_FDR)$tt
LumA_age5059 <-DEA_MTC_save(fit, "LuminalA:5059", set_logFC , set_FDR)$tt
LumA_age6069 <-DEA_MTC_save(fit, "LuminalA:6069", set_logFC , set_FDR)$tt
LumA_age70_up <-DEA_MTC_save(fit, "LuminalA:70+", set_logFC , set_FDR)$tt

dt<-decideTests(fit, p.value=set_FDR, lfc= set_logFC, adjust.method="BH")
my_coef1 = c("LuminalA:4049", "LuminalA:5059", "LuminalA:6069", "LuminalA:70+")
my_coef2 = c("LuminalB:4049", "LuminalB:5059", "LuminalB:6069", "LuminalB:70+")
my_coef3 = c("BasalLike:4049", "BasalLike:5059", "BasalLike:6069", "BasalLike:70+")
my_coef4 = c("HER2enriched:4049", "HER2enriched:5059", "HER2enriched:6069", "HER2enriched:70+")
my_coef5 = c("NormalLike:4049", "NormalLike:5059", "NormalLike:6069", "NormalLike:70+")

print (summary(dt)[,c(my_coef1)])
print (summary(dt)[,c(my_coef2)])
print (summary(dt)[,c(my_coef3)])
print (summary(dt)[,c(my_coef4)])
print (summary(dt)[,c(my_coef5)])



#set parameters!
set_logFC=2
set_FDR=0.001

## topTable PAM50 model
LumAvsLumB_tt <-DEA_MTC_save(fit, "LumAvsLumB", set_logFC , set_FDR)$tt
LumAvsBasal_tt <-DEA_MTC_save(fit, "LumAvsBasal", set_logFC , set_FDR)$tt
LumAvsHER2_tt <-DEA_MTC_save(fit, "LumAvsHER2", set_logFC , set_FDR)$tt
LumAvsNormLike_tt <-DEA_MTC_save(fit, "LumAvsNormLike", set_logFC , set_FDR)$tt
LumAvsNormal_tt <-DEA_MTC_save(fit, "LumAvsNormal", set_logFC , set_FDR)$tt
#
LumBvsBasal_tt <-DEA_MTC_save(fit, "LumBvsBasal", set_logFC , set_FDR)$tt
LumBvsHER2_tt <-DEA_MTC_save(fit, "LumBvsHER2", set_logFC , set_FDR)$tt
LumBvsNormLike_tt <-DEA_MTC_save(fit, "LumBvsNormLike", set_logFC , set_FDR)$tt
LumBvsNormal_tt <-DEA_MTC_save(fit, "LumBvsNormal", set_logFC , set_FDR)$tt
#
BasalvsHER2_tt <-DEA_MTC_save(fit, "BasalvsHER2", set_logFC , set_FDR)$tt
BasalvsNormLike_tt <-DEA_MTC_save(fit, "BasalvsNormLike", set_logFC , set_FDR)$tt
BasalvsNormal_tt <-DEA_MTC_save(fit, "BasalvsNormal", set_logFC , set_FDR)$tt
#
HER2vsNormLike_tt <-DEA_MTC_save(fit, "HER2vsNormLike", set_logFC , set_FDR)$tt
HER2vsNormal_tt <-DEA_MTC_save(fit, "HER2vsNormal", set_logFC , set_FDR)$tt
#
NormalLikevsNormal_tt <-DEA_MTC_save(fit, "NormalLikevsNormal", set_logFC , set_FDR)$tt

## decideTests
dt <-DEA_MTC_save(fit, "LumAvsLumB", set_logFC , set_FDR)$dt #coef can be anything, report dt fot all anyway
summary(dt)

## topTable stages model
set_logFC=1
set_FDR=0.05
Stage1vsNorm_tt <- DEA_MTC_save(fit, "Stage1vsNorm", set_logFC , set_FDR)$tt
Stage2vsNorm_tt <- DEA_MTC_save(fit, "Stage2vsNorm", set_logFC , set_FDR)$tt
Stage3vsNorm_tt <- DEA_MTC_save(fit, "Stage3vsNorm", set_logFC , set_FDR)$tt
Stage4vsNorm_tt <- DEA_MTC_save(fit, "Stage4vsNorm", set_logFC , set_FDR)$tt
Stage1vsStage2_tt <- DEA_MTC_save(fit, "Stage1vsStage2", set_logFC , set_FDR)$tt
Stage1vsStage3_tt <- DEA_MTC_save(fit, "Stage1vsStage3", set_logFC , set_FDR)$tt
Stage1vsStage4_tt <- DEA_MTC_save(fit, "Stage1vsStage4", set_logFC , set_FDR)$tt
Stage2vsStage3_tt <- DEA_MTC_save(fit, "Stage2vsStage3", set_logFC , set_FDR)$tt
Stage2vsStage4_tt <- DEA_MTC_save(fit, "Stage2vsStage4", set_logFC , set_FDR)$tt
Stage3vsStage4_tt <- DEA_MTC_save(fit, "Stage3vsStage4", set_logFC , set_FDR)$tt

dt <-DEA_MTC_save(fit, "Stage1vsNorm", set_logFC , set_FDR)$dt #coef can be anything, report dt fot all anyway
summary(dt)

## toprable group 1
set_logFC=1
set_FDR=0.05

Basal_st1vs4_tt <- DEA_MTC_save(fit, "Basal_st1vs4",set_logFC , set_FDR)$tt
LumA_st1vs4_tt <- DEA_MTC_save(fit, "LumA_st1vs4",set_logFC , set_FDR)$tt
allBasal_tt <- DEA_MTC_save(fit, "allBasal",set_logFC , set_FDR)$tt
allLumA_tt <- DEA_MTC_save(fit, "allLumA",set_logFC , set_FDR)$tt
Basal_stage1_tt <- DEA_MTC_save(fit,"Basal_stage1", set_logFC , set_FDR)$tt
Basal_stage4_tt <- DEA_MTC_save(fit,"Basal_stage4", set_logFC , set_FDR)$tt
LumA_stage1_tt <- DEA_MTC_save(fit,"LumA_stage1", set_logFC , set_FDR)$tt
LumA_stage4_tt <- DEA_MTC_save(fit,"LumA_stage4", set_logFC , set_FDR)$tt
BasalVSLumA_stage1_tt <- DEA_MTC_save(fit, "BasalVSLumA_stage1", set_logFC , set_FDR)$tt
BasalVSLumA_stage4_tt <- DEA_MTC_save(fit, "BasalVSLumA_stage4", set_logFC , set_FDR)$tt
average_stage4_tt <- DEA_MTC_save(fit, "average_stage4", set_logFC , set_FDR)$tt
LumA_Basal_diff_tt <- DEA_MTC_save(fit, "LumA_Basal_diff", set_logFC , set_FDR)$tt

dt <-DEA_MTC_save(fit, "allBasal", set_logFC , set_FDR)$dt #coef can be anything, report dt fot all anyway
summary(dt)



#### visualisation of limma results####

## manual volcano #### NEEDS fixing

#ggplot(LumAvsBasal_tt, aes(x=logFC, y=-log10(P.Value),col=p.adjust(P.Value, method="BH")<0.001)) +
#  geom_point(alpha=0.25) + 
#  geom_vline(aes(xintercept=2), col="blue") + 
#  geom_vline(aes(xintercept=-2), col="blue") + 
#  scale_color_manual(values=c("black", "red")) + 
#  theme(legend.position="none")



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




#volcano plot ########## RUN on server#####

#dataDEG <- read.csv("both_PAM_DEG_23-04.csv", sep=",", quote = "\n", row.names = 1)

TCGAVisualize_volcano(dataDEG$logFC,dataDEG$adj.P.Val,
                      filename = "volcanoplot_both_PAM_DEG_23-04.png",
                      x.cut = 5,y.cut = 10^-7,
                      names = rownames(dataDEG),
                      color = c("black","red","darkgreen"),
                      names.size = 2,
                      xlab = " Gene expression fold change (Log2)",
                      legend = "State",
                      title = "Volcano plot",width = 6,height = 4)#



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


#group1
all_genesDE <- unique(sort(c(getGeneNames(Basal_st1vs4_tt, 'both'),
                             getGeneNames(LumA_st1vs4_tt,'both'),
                             getGeneNames(allBasal_tt, 'both'),
                             getGeneNames(allLumA_tt,'both'),
                             getGeneNames(Basal_stage1_tt,'both'),
                             getGeneNames(Basal_stage4_tt,'both'),
                             getGeneNames(LumA_stage1_tt,'both'),
                             getGeneNames(LumA_stage4_tt,'both'),
                             getGeneNames(BasalVSLumA_stage1_tt,'both'),
                             getGeneNames(BasalVSLumA_stage4_tt,'both'),
                             getGeneNames(average_stage4_tt,'both'),
                             getGeneNames(LumA_Basal_diff_tt,'both')
                             )))



#stages
all_genesDE <- unique(sort(c(getGeneNames(Stage1vsNorm_tt,'up'),getGeneNames(Stage1vsNorm_tt,'down'),
                             getGeneNames(Stage2vsNorm_tt,'up'),getGeneNames(Stage2vsNorm_tt,'down'),
                             getGeneNames(Stage3vsNorm_tt,'up'),getGeneNames(Stage3vsNorm_tt,'down'),
                             getGeneNames(Stage4vsNorm_tt,'up'),getGeneNames(Stage4vsNorm_tt,'down'),
                             getGeneNames(Stage1vsStage2_tt,'up'),getGeneNames(Stage1vsStage2_tt,'down'),
                             getGeneNames(Stage1vsStage3_tt,'up'),getGeneNames(Stage1vsStage3_tt,'down'),
                             getGeneNames(Stage1vsStage4_tt,'up'),getGeneNames(Stage1vsStage4_tt,'down'),
                             getGeneNames(Stage2vsStage3_tt,'up'),getGeneNames(Stage2vsStage3_tt,'down'),
                             getGeneNames(Stage2vsStage4_tt,'up'),getGeneNames(Stage2vsStage4_tt,'down'),
                             getGeneNames(Stage3vsStage4_tt,'up'),getGeneNames(Stage3vsStage4_tt,'down'))))
  
                             #
#pam50
all_genesDE <- unique(sort(c(getGeneNames(LumAvsLumB_tt,'up'),getGeneNames(LumAvsLumB_tt,'down'),
                             getGeneNames(LumAvsBasal_tt,'up'),getGeneNames(LumAvsBasal_tt,'down'),
                             getGeneNames(LumAvsHER2_tt,'up'),getGeneNames(LumAvsHER2_tt,'down'),
                             getGeneNames(LumAvsNormLike_tt,'up'),getGeneNames(LumAvsNormLike_tt,'down'),
                             getGeneNames(LumAvsNormal_tt,'up'),getGeneNames(LumAvsNormal_tt,'down'),
                             #
                             getGeneNames(LumBvsBasal_tt,'up'),getGeneNames(LumBvsBasal_tt,'down'),
                             getGeneNames(LumBvsHER2_tt,'up'),getGeneNames(LumBvsHER2_tt,'down'),
                             getGeneNames(LumBvsNormLike_tt,'up'),getGeneNames(LumBvsNormLike_tt,'down'),
                             getGeneNames(LumBvsNormal_tt,'up'),getGeneNames(LumBvsNormal_tt,'down'),
                             #
                             getGeneNames(BasalvsHER2_tt,'up'),getGeneNames(BasalvsHER2_tt,'down'),
                             getGeneNames(BasalvsNormLike_tt,'up'),getGeneNames(BasalvsNormLike_tt,'down'),
                             getGeneNames(BasalvsNormal_tt,'up'),getGeneNames(BasalvsNormal_tt,'down'),
                             #
                             getGeneNames(HER2vsNormLike_tt,'up'),getGeneNames(HER2vsNormLike_tt,'down'),
                             getGeneNames(HER2vsNormal_tt,'up'),getGeneNames(HER2vsNormal_tt,'down'),
                             #
                             getGeneNames(NormalLikevsNormal_tt,'up'),getGeneNames(NormalLikevsNormal_tt,'down'))))
                             # )))
length(all_genesDE)

### quick compate with autophagy
autophagy_genes<- as.vector(read.table("autopahagy_genes.txt", as.is = T, header = FALSE))$V1
shared <- intersect(autophagy_genes,all_genesDE)
print(paste0("Total number of DE genes: ", length(all_genesDE)))
print(paste0("Autophagy genes: ", length(shared)))



# Get expression values
v<-DEA_limmaVoom_out$v
all_geneExp <- v$E[rownames(v$E) %in% all_genesDE , ]
dim(all_geneExp)
#View(all_geneExp)

#write.table(all_L, paste0(my.dir, my.pattern, "_DE_Limma.txt"), sp = "\t", quote = FALSE)
library(gplots)
mycol <- colorpanel(1000,"blue","white","red") 
heatmap.2(all_geneExp , scale="none", labRow=all_genesDE, 
          labCol=y$samples$PAM50, col=mycol, trace="none", density.info="none", 
          margin=c(8,6), lhei=c(2,10), dendrogram="column")



###### patient classes setup #####

#get only interesting cols
names(y$samples)
hm.design<-subset(y$samples[,c('morphology','stages','PAM50')])#, 'tss'

#recorder annotations
hm.design$PAM50 = factor(hm.design$PAM50, levels = c("Basal-like", "HER2-enriched", "Luminal A", "Luminal B", "Normal-like", "Normal"))
hm.design$stages = factor(hm.design$stages, levels = c("stage1", "stage2", "stage3", "stage4", "unknown")) 

hm.design$morphology = factor(hm.design$morphology, levels = c("Normal",  "Mucinous adenocarcinoma", 
                                       "Ductal carcinoma","Lobular carcinoma", "Ductual mixed with others ",
                                       "Ductal and lobular mixed", "Metaplastic carcinoma" ,"Other")) #other is a mix of few samples of different ones

# colours
PAM50Col <- c( "#E69F00", "#0072B2", "#F0E442","red3", "#CC79A7","#009E73")##D55E00->altred
tumourTypCol = c("#009E73",  "#0072B2", "#D55E00","#56B4E9","#E69F00","#F0E442","#CC79A7", "white")
tumourStgCol = c( "red","green","blue", "yellow2", "white")

names(PAM50Col)<-levels(hm.design$PAM50)
names(tumourTypCol)<-levels(hm.design$morphology)
names(tumourStgCol)<-levels(hm.design$stages)

annColour <-list(
  PAM50=PAM50Col,
  morphology=tumourTypCol,
  stages=tumourStgCol
  )
######


pheatmap::pheatmap(mat = as.matrix(all_geneExp), color = brewer.pal(name = "YlGnBu", n = 9),
                   clustering_distance_rows = 'manhattan', 
                   clustering_distance_cols = 'manhattan', 
                   #scale="row",
                   annotation_col=hm.design,
                   annotation_colors = annColour,
                   #annotation_row=ok_g_functions,
                   cluster_cols = T, cluster_rows = T, 
                   show_rownames = F,show_colnames = F,
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
plotMD(efit, column=1, status=dt[,5], main=colnames(efit)[5])#, xlim=c(-8,13))
plotMD(efit, column=1, status=dt[,10], main=colnames(efit)[10])#, xlim=c(-8,13))
plotMD(efit, column=1, status=dt[,15], main=colnames(efit)[15])#, xlim=c(-8,13))




