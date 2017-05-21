#! /usr/bin/env Rscript
#
rm(list=ls(all=TRUE))

library(TCGAbiolinks)
library(SummarizedExperiment)
library(plyr)
require("getopt", quietly=TRUE)

library(limma)
library(sva)


# User input:
#1) -d --data: Count data set to be used (.rda-file containing SummarizedExperiment object consisting of: genes (rows) & samples (columns))
#2) -c --cutoff: cor.cutoff value
#3) -n --normMethod: Normalization method ("gcContent", "geneLength")
#4) -q --qntCut: qnt.cut value
#5) -f --filtMethod: Filtering method ("quantile", "varFilter", "filter1", "filter2")
#6) -b --batchCorrect: Correcting for batch effect, if argument is set to TRUE
#7) -r --removePairedSamples: removing paired samples if argument is set to TRUE
#8) -t --transformation: transform count data to log2-counts per million (logCPM) by using voom.

# flag specification # ref: https://tgmstat.wordpress.com/2014/05/21/r-scripts/
spec = matrix(c(
  "data", "d", 1, "character",
  "corCut", "c", 2, "double",
  "normMethod", "n", 2, "character",
  "qntCut", "q", 2, "double",
  "filtMethod", "f", 2, "character",
  "batchCorrect", "b", 2, "logical",
  "removePairedSamples", "r", 2, "logical",
  "transformation", "t", 2, "logical"
), byrow=TRUE, ncol=4)

opt = getopt(spec)

argsData <- opt$data


if (is.null(opt$corCut)) {
  argsCorCut <- 0.6 # default cor.cut value.
} else {
  argsCorCut <- opt$corCut
}

if (is.null(opt$normMethod)) {
  argsNormMethod <- "both" # default normalization method selection.
} else {
  argsNormMethod <- opt$normMethod
}

if (is.null(opt$qntCut)) {
  argsQntCut <- 0.25 # default qnt.cut value.
} else {
  argsQntCut <- opt$qntCut
}

if (is.null(opt$filtMethod)) {
  argsFiltMethod <- "quantile" # default filtering method selection.
} else {
  argsFiltMethod <- opt$filtMethod
}

if (is.null(opt$removePairedSamples)) {
  argsRemovePairedSamples <- TRUE # default: paired samples are removed.
} else {
  argsRemovePairedSamples <- opt$removePairedSamples
}

if (is.null(opt$batchCorrect)) {
  argsBatchCorrect <- TRUE # default: batch effect is corrected.
} else {
  argsBatchCorrect <- opt$batchCorrect
}

if (is.null(opt$transformation)) {
  argsTransformation <- TRUE # default: voom transformation.
} else {
  argsTransformation <- opt$transformation
}

# Thildes function to get batch number
get_IDs <- function(data) {
  IDs <- strsplit(c(colnames(data)), "-")
  IDs <- ldply(IDs, rbind)
  colnames(IDs) <- c('project', 'tss','participant', 'sample', "portion", "plate", "center")
  cols <- c("project", "tss", "participant")
  IDs$patient <- apply(IDs[,cols],1,paste,collapse = "-" )
  barcode <- colnames(data)
  IDs <- cbind(IDs, barcode)
  condition <- gsub("11+[[:alpha:]]", "normal", as.character(IDs$sample))
  condition  <- gsub("01+[[:alpha:]]", "cancer", condition)
  IDs$condition <- condition
  IDs$myorder  <- 1:nrow(IDs)
  return(IDs)
}




# -----------------------------------------------------------------------------------------------------------------------------
# Load and preprocess data
# -----------------------------------------------------------------------------------------------------------------------------
# load data
dataAssy<-get(load(argsData))

# preprocessing
dataPrep <- TCGAanalyze_Preprocessing(object = dataAssy,
                                      cor.cut = argsCorCut,
                                      datatype = "raw_count")

# extract metadata from dataAssy
sampleMetaData <- colData(dataAssy)




# -----------------------------------------------------------------------------------------------------------------------------
# Remove paired samples if set to TRUE by the user (This step only removes the TP part).
# -----------------------------------------------------------------------------------------------------------------------------
if(argsRemovePairedSamples){
  print("Removing paired samples from count data...")
  paired <- TCGAquery_MatchedCoupledSampleTypes(sampleMetaData$barcode,c("NT","TP")) # find paired samples
  paired_tumor_index <- lapply(paired, FUN = function(x) substr(unlist(strsplit(x, split = "-"))[[4]],1,2) == "01") # find index for TP part of paired samples
  paired_tumor <- paired[unlist(paired_tumor_index)] # find TP samples of paired samples
  sampleMetaData <- sampleMetaData[!sampleMetaData$barcode %in% paired_tumor, ] # update sampleMetaData with only unpaired samples
  
  
  
  # save sampleMetaData with only unpaired samples
  write.table(subset(sampleMetaData, select = c("barcode", "shortLetterCode")),"sampleMetaData_pairedSamplesRemoved.txt", quote=F,sep = ",", row.names = F)
  
  # update dataPrep with only unpaired samples
  dataPrep <- dataPrep[,!colnames(dataPrep) %in% paired_tumor]
}





# -----------------------------------------------------------------------------------------------------------------------------
# Normalization
# -----------------------------------------------------------------------------------------------------------------------------
if(argsNormMethod != "non"){
  if(argsNormMethod == "both"){
    dataNorm_ <- TCGAanalyze_Normalization(tabDF = dataPrep,
                                           geneInfo = geneInfo,
                                           method = "gcContent")
    
    dataNorm <- TCGAanalyze_Normalization(tabDF = dataNorm_,
                                          geneInfo = geneInfo,
                                          method = "geneLength")
    
    print("Normalization done")
    
  }else{
    dataNorm <- TCGAanalyze_Normalization(tabDF = dataPrep,
                                          geneInfo = geneInfo,
                                          method = argsNormMethod)
    
    print("Normalization done")
  }
}else{
  dataNorm <- dataPrep
}

# convert row names from entrez-gene into gene names in the data frame, if no normalization method is used
#(since this happens automatically in the normalization method).
if(argsNormMethod == "non"){
  rownames(dataNorm)<-unlist(lapply(rownames(dataNorm), FUN = function(names) unlist(strsplit(names, "\\|"))[1]))
}



# -----------------------------------------------------------------------------------------------------------------------------
# Filtering
# -----------------------------------------------------------------------------------------------------------------------------
if(argsFiltMethod != "non"){
  print("Performing filtering...")
  dataFilt <- TCGAanalyze_Filtering(tabDF = dataNorm,
                                    method = argsFiltMethod,
                                    qnt.cut =  argsQntCut)
  
  print("Filtering done")
}else{
  dataFilt <- dataNorm
}





# get metadata from Thildes function
meta<-get_IDs(dataFilt)

# create mod_design object for voom transformation & batch correction
mod_design <- model.matrix(~as.factor(meta$condition))
colnames(mod_design) <- c("cancer", "normal")
#   cancer normal
#1      1      0
#2      1      1
#3      1      0


# -----------------------------------------------------------------------------------------------------------------------------
# voom transformation
# -----------------------------------------------------------------------------------------------------------------------------


# voom transformation: transform count data to log2-counts per million (logCPM)
if(argsTransformation){
  print("Performing voom transformation: Transform count data to log2-counts per million (logCPM)...")
  
  v <- voom(dataFilt, mod_design, plot=TRUE)
  dataFilt <- v$E
  
  print("voom transformation done")
}


# save preprocessed data set without batch correction (_Preprocessed_wo_batch.rda), to be used in the DEA script.
save(dataFilt,file=paste0(dirname(argsData),"/",unlist(strsplit(basename(argsData),".", fixed = T))[1],"_PreprocessedData_wo_batch.rda"))
print("Preprocessed data without batch correction saved.")




# -----------------------------------------------------------------------------------------------------------------------------
# Batch correction
# -----------------------------------------------------------------------------------------------------------------------------
if(argsBatchCorrect){
  print("Batch correction...")
  
  
  # Batch correction
  batch_corr <- ComBat(dataFilt, as.factor(meta$plate), mod_design, par.prior=TRUE,prior.plots=FALSE)
  dataFilt <- batch_corr
}






# -----------------------------------------------------------------------------------------------------------------------------
# Save preprocessed data to .rda-file
# -----------------------------------------------------------------------------------------------------------------------------
save(dataFilt,file=paste0(dirname(argsData),"/",unlist(strsplit(basename(argsData),".", fixed = T))[1],"_PreprocessedData.rda"))

print("Data saved.")