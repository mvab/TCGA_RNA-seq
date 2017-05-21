#! /usr/bin/env Rscript
rm(list=ls(all=TRUE))

library(TCGAbiolinks)
library(SummarizedExperiment)
library(plyr)
require("getopt", quietly=TRUE)

library(limma)
library(sva)


# User input (old):
#1) -i --projectID: the project of interest (e.g. TCGA-LUAD).
#2) -d --data: Count data set to be used (.rda-file containing SummarizedExperiment object consisting of: genes (rows) & samples (columns))
#3) -1 --condition1: a .txt containing a list of samples in question with newline as delimiter (\n).
#4) -2 --condition2: a .txt containing a list of samples in question with newline as delimiter (\n).
# IF conditions are not specified, Tumour/Normal conditions will be used. 

########## NB the current DEA works on no-voom and no-batch effect data!


# flag specification # ref: https://tgmstat.wordpress.com/2014/05/21/r-scripts/
spec = matrix(c(
  "projectID", "i", 1, "character",
  "data", "d", 1, "character",
  "condition1", "1", 2, "character",
  "condition2", "2", 2, "character"
), byrow=TRUE, ncol=4)


opt = getopt(spec)

argsProjectID <- opt$projectID
argsData <- opt$data
dataSE<-get(load(argsData))


if (is.null(opt$condition1)) {
  argsCond1 <- FALSE # default: will be assigned to TP
} else {
  argsCond1 <- opt$condition1
}

if (is.null(opt$condition2)) {
  argsCond2 <- FALSE # default: will be assigned to NT
} else {
  argsCond2 <- opt$condition2
}



getPairedData <- function(argsProjectID){
  query.exp <- GDCquery(project = argsProjectID,
                        legacy = TRUE,
                        data.category = "Gene expression",
                        data.type = "Gene expression quantification",
                        platform = "Illumina HiSeq", 
                        file.type = "results",
                        experimental.strategy = "RNA-Seq",
                        sample.type = c("Primary solid Tumor","Solid Tissue Normal"))
  
  # select samples that are primary solid tumor
  dataSmTP <- TCGAquery_SampleTypes(query.exp$results[[1]]$cases,"TP")
  
  # select samples that are solid tissue normal
  dataSmNT <- TCGAquery_SampleTypes(query.exp$results[[1]]$cases,"NT")
  
  return(list(TP=dataSmTP, NT=dataSmNT))
}


# load sample list if condition flags (-1 and -2) are used and save them as vectors.
if( (argsCond1 != FALSE ) & (argsCond2 != FALSE) ) {
  argsCond1Samples<-read.table(argsCond1, as.is = T, header = FALSE) 
  argsCond1Samples<-as.vector(argsCond1Samples)
  argsCond1Samples<-argsCond1Samples$V1
  
  argsCond2Samples<-read.table(argsCond2, as.is = T, header = FALSE) 
  argsCond2Samples<-as.vector(argsCond2Samples)
  argsCond2Samples<-argsCond2Samples$V1
}else{
  pairedSamples<-getPairedData(argsProjectID)
  argsCond1Samples<-pairedSamples$TP
  argsCond2Samples<-pairedSamples$NT
  argsCond1 <- "Tumour"
  argsCond2 <- "Normal"
}


print ("Extraction of selected groups is finished.")
cat("\n")
print ("################### Data summary: ####################")
print (paste0("Preprocessed data contains ",dim(dataSE)[1]," genes and ", dim(dataSE)[2], " samples. "))
print (paste0("There are ", length(argsCond1Samples), " samples in condition 1 (", argsCond1, ")"))
print (paste0("There are ", length(argsCond2Samples), " samples in condition 2 (", argsCond2, ")"))
print ("######################################################")
cat("\n")
print("Starting DEA ...")


# DEA anaylsis from TCGAbiolinks ---- to be replaced!
dataDEGs <- TCGAanalyze_DEA(mat1 = dataSE[,argsCond1Samples],
                            mat2 = dataSE[,argsCond2Samples],
                            Cond1type = argsCond1,
                            Cond2type = argsCond2,
                            fdr.cut = 0.01 ,
                            logFC.cut = 1,
                            method = "glmLRT")  

#print (dataDEGs)
