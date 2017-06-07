library(Mfuzz)
library(edgeR)

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
getGeneLenght_out<-getGeneLenght(dataSE)
gene_lengths<-getGeneLenght_out$gene_lengths
dataSE<-getGeneLenght_out$dataSE
dim(dataSE)
dim(gene_lengths)

# edgeR object
dge <- DGEList(counts=dataSE, genes = data.frame(Length= as.numeric(gene_lengths$gene_length)) ) #here will be adding genes: genes = Ann)
#head(dge$genes)

#### RPKM
yr <- calcNormFactors(dge)
# rpkm will use the normalized effective library sizes to compute rpkm instead of the raw library sizes. 
RPKM <- rpkm(yr)
dim(RPKM) #17368   984 matrix

# get genes that have exp >1 rpkm in at least 12 samples
keep.exprs<-c()
filt.exprs <- rowSums(RPKM > 1) >=12 
for (i in 1:length(filt.exprs )){
  if (as.vector(filt.exprs [i]) == TRUE) {
    keep.exprs<-append(keep.exprs, names(filt.exprs [i]))
  }
}
length(keep.exprs)

RPKM <- RPKM[rownames(RPKM) %in% keep.exprs, ]
dim(RPKM)

##### if want to look at specific subtype


keep_subtype<-samples.matrix[samples.matrix$PAM50=="Normal-like" | samples.matrix$PAM50=="Normal",]$barcode
RPKM<- RPKM[,colnames(RPKM) %in% keep_subtype]
dim(RPKM)


# remove samples whose stage is unknowm
known_stage<-samples.matrix[samples.matrix$tumourStages!='unknown',]$barcode
RPKM<- RPKM[,colnames(RPKM) %in% known_stage]
dim(RPKM)

# create a rpkm matrix that has a column per 'time point', i.e. stage
getAverageExpSampleSet<- function(dataSE, set){
  
  #get list of samples for this set
  set_samples<-samples.matrix[samples.matrix$tumourStages==set,]$barcode
  
  #get their exp vals
  dataSE<- dataSE[,colnames(dataSE) %in% set_samples]
  print(paste0(set, " has ", dim(dataSE)[2], " samples."))
  
  #caclulate average expression 
  rpkm_set_mean <- as.matrix(rowMeans(dataSE))
  colnames(rpkm_set_mean) <- set
  
  return(rpkm_set_mean)
}


#calculate mean exp for each stage
stage_averages<- matrix(nrow=dim(RPKM)[1], ncol=5)
colnames(stage_averages)<-c('Normal', 'stage1','stage2','stage3','stage4')
rownames(stage_averages)<- rownames(RPKM)

for (i in colnames(stage_averages)){
  this_average<-getAverageExpSampleSet(RPKM, i)
  stage_averages[,i]<- this_average}
head(stage_averages)



#### extract data only for autophagy genes
getAutophagyGenes <- function(dataSE){
  
  autophagy_genes<- as.vector(read.table
    ("~/Bioinformatics MSc UCPH/0_MasterThesis/TCGAbiolinks/CBL_scripts/data/transcription_factors.txt", as.is = T, header = FALSE))$V1
  all_genes<-rownames(dataSE)
  shared <- intersect(autophagy_genes,all_genes)
  newdataSE<- dataSE[c(shared),]
  
  print(paste0("Total number of genes: ", length(all_genes)))
  print(paste0("Autophagy genes: ", length(shared)))
  
  return(newdataSE)
} 
dim (stage_averages)
stage_averages <- getAutophagyGenes(stage_averages)
dim (stage_averages)



#we create the Mfuzz object (ExpressionSet)
exprSet=ExpressionSet(assayData=stage_averages)



#As a first step,  we exclude genes with more than 25% of the measurements missing
#-> genes with 0 RPKM in 25% of the conditions
exprSet.r <- filter.NA(exprSet, thres=0.25) 

#Fuzzy  c-means  like  many  other  cluster  algorithms,  does  not  allow  for  missing  values.
#Thus,  we  timelace  remaining  missing  values  by  the  average  values  expression  value
#of  the corresponding gene.
#Methods for replacement of missing values. Missing values should be indicated by NA in #the expression matrix.
#Mode method for replacement of missing values:
#mean- missing values will be replaced by the mean expression value of the gene,
#median- missing values will be replaced by the median expression value of the gene,
#knn- missing values will be replaced by the averging over the corresponding expression    #values of the k-nearest neighbours,
#knnw-same replacement method as knn, but the expression values averaged are #weighted by the distance to the corresponding neighbour

#replace missing values with means
exprSet.f <- fill.NA(exprSet.r,mode="mean") 


#As soft clustering is noise robust, pre-filtering can usually be avoided.
#However, if the number of genes with small expression changes is large, such pre-filtering #may be necessary to reduce noise.
#This function can be used to exclude genes with low standard deviation.
#min.std : threshold for minimum standard deviation.
#If the standard deviation of a gene's expression is smaller than min.std the corresponding #gene will be excluded.

#Calculation the standard deviation shows, however,
#that the transition between low and high values for variation in 
#gene expression is smooth and no particular cut-off point is indicated 
tmp <- filter.std(exprSet.f,min.std=0.3, visu=FALSE) #genes with lower std will be excluded


#setdiff(rownames(exprSet.f), rownames(tmp))





#Since the clustering is performed in Euclidian space, the expression values of 
#genes were standardised to have a mean value of zero and a standard deviation of one. 
#Importantly, Mfuzz assumes that the given expression data are fully preprocessed 
#including  any  data  normalisation.
#The  function standardise does  not  replace  the  normalisation step (eg RPKN #normalization).
#This step ensures that vectors of genes with similar changes in expression are close
#in Euclidean space:

exprSet.s <- standardise(tmp)  #NB this step makes expression between genes COMPARABLE
dim(exprSet.s)
# c is the number of clusters 
# m is the fuzziness parameter ->

# we would like to choose a value which prevents clustering of random data. 
#Note, that fuzzy clustering can be tuned in such manner, that random data is not 
#clustered. This is a clear advantage to hard clustering (such as k-means)


m1=mestimate(exprSet.s) 
m1



current_test <- "all_genes_all_samples"
current_test <- "autophagy_genes_all_samples"
current_test <- "autocore_genes_all_samples"
current_test <- "transcrfactors_all_samples"


#c1 =cselection(exprSet.s, m=m1, crange=seq(4,32,2),repeats=5,visu=TRUE)

setwd("~/Bioinformatics MSc UCPH/0_MasterThesis/TCGAbiolinks/CBL_scripts/data/mfuzz_with_genedata/NormalLike_samples/")

# for preview
cl <- mfuzz(exprSet.s,c=3,m=m1) 
mfuzz.plot(exprSet.s,cl=cl,mfrow=c(3,3),min.mem=0, time.labels=colnames(stage_averages), new.window = T)

#for saving
pdf(file=paste0(current_test,".pdf"))
mfuzz.plot(exprSet.s,cl=cl,mfrow=c(3,3),min.mem=0, time.labels=colnames(stage_averages), new.window = F)
dev.off()

pdf(file=paste0(current_test,"_members_above_0.6_membership.pdf"))
mfuzz.plot(exprSet.s,cl=cl,mfrow=c(3,3),min.mem=0.6, time.labels=colnames(stage_averages), new.window = F)
dev.off()
#par(mfrow=c(1,1))
#O <- overlap(cl) 
#ov.tmp <- overlap.plot(cl,over=O,thres=0.05)


#To extract list of genes belonging to the cluster cores, the acore function can be used.
cores<-acore(exprSet.s,cl=cl,min.acore=0.6)

#make a list object with gene names for each cluster
cluster_data<-list()
# get total number of genes after filetering : std=0.3, membership  >0.6
gene_sum=0

for ( i in 1:dim(summary(cores))[1]){
  print (paste( "Cluster", i, ":" , dim(cores[[i]])[1] ))
  this_cluster<-as.character(paste0("cluster",i))
  cluster_data[[this_cluster]]<-rownames(cores[[i]])
  gene_sum=gene_sum+dim(cores[[i]])[1]
}

cluster_data
gene_sum

#clusters don't overlap!

save(cluster_data, file= paste0(current_test, "_cluster_data.rda"))



###############################################

#testing loading the data:
data<-get(load("autophagy_genes_all_samples_cluster_data.rda"))

#explore data list (run this loop)

gene_sum=0
for ( i in names(data)){
  print (paste( i, ":" , length(data[[i]]) ))
  gene_sum=gene_sum+dim(cores[[i]])[1]
}


# total number of genes after filtering and clustring in this data group (std=0.3, membership  >0.6)
gene_sum

#e.g.
cluster1<-data[[1]]

