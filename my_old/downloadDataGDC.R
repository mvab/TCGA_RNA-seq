#! /usr/bin/env Rscript


rm(list=ls(all=TRUE))

library(TCGAbiolinks)
library(SummarizedExperiment)
require("getopt", quietly=TRUE)




# User input (old):
#1) -i --projectID: the project of interest (e.g. TCGA-LUAD).
#2) -l --legacy: specify whether legacy data of newest data is used.
#3) -c --dataCategory: specify the data category of interest.
#4) -t --dataType: the data type of interest.
#5) -p --platform: the platform type of interest.
#6) -f --fileType: the fileType of interest.
#7) -w --workflow: the workflow of interest.
#8) -e --expStrategy: the experimental strategy of interest.
#9) -s --sampleList: a .txt containing a list of samples in question with newline as delimiter (\n).
#10) -o --pairedOnly: TRUE) Only barcodes which contain paired samples (TP and NT) will be returned. FALSE) All barcodes will be returned.



# flag specification # ref: https://tgmstat.wordpress.com/2014/05/21/r-scripts/
spec = matrix(c(
  "prosonectID", "i", 1, "character",
  "legacy", "l", 2, "logical",
  "dataCategory", "c", 2, "character",
  "dataType", "t", 2, "character",
  "platform", "p", 2, "character",
  "fileType", "f", 2, "character",
  "workflow", "w", 2, "character",
  "expStrategy", "e", 2, "character",
  "sampleList", "s", 2, "logical",
  "pairedOnly", "o", 2, "logical"
), byrow=TRUE, ncol=4)



opt = getopt(spec)

argsProjectID <- opt$projectID

if (is.null(opt$legacy)) {
  argsLegacy <- TRUE # default data category.
} else {
  argsLegacy <- opt$legacy
}

if (is.null(opt$dataCategory)) {
  argsDataCategory <- FALSE # default data category.
} else {
  argsDataCategory <- opt$dataCategory
}

if (is.null(opt$dataType)) {
  argsDataType <- FALSE # default data type.
} else {
  argsDataType <- opt$dataType
}

if (is.null(opt$platform)) {
  argsPlatform <- FALSE # default platform type.
} else {
  argsPlatform <- opt$platform
}

if (is.null(opt$fileType)) {
  argsFileType <- FALSE # default file type.
} else {
  argsFileType <- opt$fileType
}

if (is.null(opt$workFlow)) {
  argsWorkFlow <- FALSE # default work flow.
} else {
  argsWorkFlow <- opt$workFlow
}

if (is.null(opt$expStrategy)) { # NEW NEW NEW
  argsExpStrategy <- FALSE # default paired only.
} else {
  argsExpStrategy <- opt$expStrategy
}

if (is.null(opt$sampleList)) {
  argsSampleList <- FALSE # default sample list.
} else {
  argsSampleList <- opt$sampleList
}


if (is.null(opt$pairedOnly)) {
  argsPairedOnly <- FALSE # default paired only.
} else {
  argsPairedOnly <- opt$pairedOnly
}




# -----------------------------------------------------------------------------------------------------------------------------
# Download query specified by user input flags
# -----------------------------------------------------------------------------------------------------------------------------

# retrieve data info specified by the query
query.exp <- GDCquery(project = argsProjectID,
                      legacy = argsLegacy,
                      data.category = argsDataCategory,
                      data.type = argsDataType,
                      platform = argsPlatform,
                      file.type = argsFileType,
                      workflow.type = argsWorkFlow,
                      experimental.strategy = argsExpStrategy,
                      sample.type = c("Primary solid Tumor","Solid Tissue Normal"))


# get clinical data
dataClin <- GDCquery_clinic(project = argsProjectID, "clinical")

# select samples that are primary solid tumor
dataSmTP <- TCGAquery_SampleTypes(query.exp$results[[1]]$cases,"TP")

# select samples that are solid tissue normal
dataSmNT <- TCGAquery_SampleTypes(query.exp$results[[1]]$cases,"NT")




# load sample list if the -s or --sampleList flag is used and save it as a vector.
if(argsSampleList == TRUE){
  argsSamples<-read.table("sampleList.txt", as.is = T, header = FALSE)
  argsSamples<-as.vector(argsSamples$V1)
}else{
  argsSamples<-c(dataSmTP,dataSmNT)
}

# only use paired samples (NT and TP) if the -o --pairedOnly flag is set to TRUE (or T).
if(argsPairedOnly == TRUE){
  argsSamples <- TCGAquery_MatchedCoupledSampleTypes(barcode = argsSamples, c("NT","TP"))
}



# update query
query.exp <- GDCquery(project = argsProjectID,
                      legacy = argsLegacy,
                      data.category = argsDataCategory,
                      data.type = argsDataType,
                      platform = argsPlatform,
                      file.type = argsFileType,
                      workflow.type = argsWorkFlow,
                      experimental.strategy = argsExpStrategy,
                      barcode = argsSamples)


# download data specified by the query
GDCdownload(query.exp)


# -----------------------------------------------------------------------------------------------------------------------------
# save data as a data.frame into categorized subfolders
# -----------------------------------------------------------------------------------------------------------------------------
if(argsLegacy == TRUE){
  data.exp <- GDCprepare(query = query.exp, save = TRUE, save.filename = paste0("GDCdata/",argsProjectID,"/","legacy","/",gsub("[ ]","_",argsDataCategory),"/", gsub("[ ]","_",argsDataType),"/",unlist(strsplit(argsProjectID,"-"))[2],"_",gsub("[ ]","_",argsPlatform),".rda"))
}else{
  data.exp <- GDCprepare(query = query.exp, save = TRUE, save.filename = paste0("GDCdata/",argsProjectID,"/","harmonized","/",gsub("[ ]","_",argsDataCategory),"/", gsub("[ ]","_",argsDataType),"/",unlist(strsplit(argsProjectID,"-"))[2],"_",gsub("[ ]","_",argsPlatform),".rda"))
}
