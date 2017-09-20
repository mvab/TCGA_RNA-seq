
library(edgeR)
setwd("~/Bioinformatics MSc UCPH/0_MasterThesis/TCGAbiolinks/CBL_scripts/data")
source("../functions.R") 


## most recent all types and stages cancer and normal
dataSE<-get(load("BRCA_Illumina_HiSeqnew_updatedSE_allFemale_PreprocessedData_wo_batch_GC_010_updatedSE_allTypes_allStages_Female.rda"))
samples.matrix<-get(load("BRCA_Illumina_HiSeqnew_updatedSE_allFemale_PreprocessedData_wo_batch_GC_010_sampleMatrix_allTypes_allStages_Female.rda"))
dim(samples.matrix)
dim(dataSE)


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


##### SELECTION OPTIONS

#### load and extract genes from DEA testing
all_genesDE<- get(load(file="all_genesDE.rda"))
autophagy_genesDE<- get(load(file="autophagy_genesDE.rda"))

forDEA_genes<-get(load("genes_for_DEA.rda"))

dataSE<- dataSE[rownames(dataSE) %in% forDEA_genes,] #15779
dim(dataSE) 


#### extract data only for autophagy genes
getAutophagyGenes <- function(dataSE){
  
  autophagy_genes<- as.vector(read.table("autopahagy_genes.txt", as.is = T, header = FALSE))$V1
  all_genes<-rownames(dataSE)
  shared <- intersect(autophagy_genes,all_genes)
  newdataSE<- dataSE[c(shared),]
  
  print(paste0("Total number of genes: ", length(all_genes)))
  print(paste0("Autophagy genes: ", length(shared)))
  
  return(newdataSE)
} 
dim (dataSE)
dataSE <- getAutophagyGenes(dataSE)



# normalise with log
dge <- DGEList(dataSE)
dge <- calcNormFactors(object=dge, method="TMM")
dim(dge)
EM <- cpm(x=dge, log=TRUE)

# remove genes that are not use in DEG: 
#cpm <- cpm(dge, log=FALSE)
#keep.exprs <- rowSums(cpm > 1) >= 50        #A CPM value of 1 is equivalent to a log-CPM value of 0.
#keep.exprs <- rowSums(cpm > 2) >= 19
#dge <- dge[keep.exprs,, keep.lib.sizes=FALSE]
#dim(dge)
#EM <- cpm(x=dge, log=TRUE) #new

### calculate stage average for each PAM50 subtype
getAverageExpSampleSet<- function(dataSE, samples.matrix,  my_stage, my_subtype) {
  
  #get list of samples for this set
  stage_samples<-samples.matrix[(samples.matrix$tumourStages==my_stage & samples.matrix$PAM50==my_subtype),]$barcode
  
  #get their exp vals
  dataSE<- dataSE[,colnames(dataSE) %in% stage_samples]
  
  print(paste(my_stage, " ", my_subtype ," has ", dim(dataSE)[2], " samples."))
  
  #caclulate average expression 
  set_mean <- as.matrix(rowMeans(dataSE))
  colnames(set_mean) <- as.vector(paste(my_stage, my_subtype))
  
  return(set_mean)
}


#calculate mean exp for each stage
averages<- matrix(nrow=dim(EM)[1], ncol=30)
rownames(averages)<- rownames(EM)
colnames(averages)<-seq(1,30,1)

stage_list<-c('Normal', 'stage1','stage2','stage3','stage4')
subtype_list<-c("Luminal A", "Luminal B","Basal-like"  ,"HER2-enriched" ,"Normal-like","Normal")

count=0
for (i in stage_list){
  for (j in subtype_list){
    count=count+1
    colnames(averages)[count]<-as.vector(paste(i,j))
    this_average<-getAverageExpSampleSet(EM, samples.matrix, i, j)
    print (colnames(averages)[count])
    averages[,(paste(i,j))]<- this_average
    
  }
}  
dim(averages)
head(averages)  
head(averages[,c(1:5,12,18,24,29, 30)])
averages<-averages[,-c(1:5,12,18,24,29,30)] #removes cols with nan
colnames(averages)
colnames(averages)[1]<-"Normal"

my_palette2<-c("#009E73",
        "#F0E442", "red3","#E69F00", "#0072B2", "#CC79A7",
        "#F0E442", "red3","#E69F00", "#0072B2", "#CC79A7",
        "#F0E442", "red3","#E69F00", "#0072B2", "#CC79A7",
        "#F0E442", "red3","#E69F00", "#0072B2")
#library(scales)
#show_col(my_palette2)

library(rafalib)
mypar()

windowsFonts(A = windowsFont("Times New Roman"))

d <- dist( t(averages) )

hc <- hclust(d)
plot(hc,labels=as.vector(colnames(averages)),cex=0.5)
dev.new()

myplclust(hc, labels=as.vector(colnames(averages)), 
          lab.col=my_palette2,
          cex=0.8, main="Clustering of PAM50 subtype and stages groups", 
          family="serif")

