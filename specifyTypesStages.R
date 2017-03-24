#! /usr/bin/env Rscript
rm(list=ls(all=TRUE))

library(TCGAbiolinks)
library(SummarizedExperiment)
library(plyr)
require("getopt", quietly=TRUE)

library(limma)
library(sva)



#################################

#Rscript thisscript.R -i TCGA-BRCA -d brcaExp_PreprocessedData_wo_batch.rda -f female_onlyTandN.txt -t 6 -s 4  -1 female_85203.txt -2 female_84803.txt -3 female_85003.txt -4 female_85223.txt -5 #female_85233.txt -6 female_85753.txt -7 stage1.txt -8 stage2.txt -9 stage3.txt -0 stage4.txt

##################################


# User input:
#1) -i --projectID: the project of interest (e.g. TCGA-LUAD).
#2) -d --data: Count data set to be used (.rda-file containing SummarizedExperiment object consisting of: genes (rows) & samples (columns))
#3) -f --tumournormal: a .txt containing a list of samples in question (here: all female T adn N) with newline as delimiter (\n).
#4) -t --tumoutTypeCount: specify how many tumour types you will input as lists
#5) -s --stageCount: specify how many tumour stages you will input as lists
#6) -1 --condition1: a .txt containing a list of samples in question with newline as delimiter (\n).
#7) -2 --condition2: a .txt containing a list of samples in question with newline as delimiter (\n).
#8,9,10,11) -3,-4,-5,-6 more conditions (subtypes)
# 12,13,14,15) -7,-8,-9,-0 more conditons (stages)


## NB subtypes must always come before stages, unless only stages are used. 


# flag specification # ref: https://tgmstat.wordpress.com/2014/05/21/r-scripts/
spec = matrix(c(
  "projectID", "i", 1, "character",
  "data", "d", 1, "character",
  "tumournormal", "f", 1, "character", # extracts all FEMALE t/n samples. #specify this to get all normal samples
  "tumourTypeCount", "t", 2, "double", #how many subtypes of tumours youre going to specify
  "stageCount", "s", 2, "double",   # how many stages youre going to specigy
  "condition1", "1", 2, "character", #takes up to 10 conditions
  "condition2", "2", 2, "character",
  "condition3", "3", 2, "character",
  "condition4", "4", 2, "character",
  "condition5", "5", 2, "character",
  "condition6", "6", 2, "character",
  "condition7", "7", 2, "character",
  "condition8", "8", 2, "character",
  "condition9", "9", 2, "character",
  "condition0", "0", 2, "character"  
), byrow=TRUE, ncol=4)


opt = getopt(spec)

argsProjectID <- opt$projectID
argsData <- opt$data
dataSE<-get(load(argsData))

conditions_specified <-vector(mode="character", length=0)


if (is.null(opt$tumourTypeCount)) {
  argsTTcount <- 0
} else {
  argsTTcount <- opt$tumourTypeCount 
}
if (is.null(opt$stageCount)) {
  argsScount <- 0
} else {
  argsScount <- opt$stageCount 
}


if (is.null(opt$condition1)) {
  argsCond1 <- FALSE 
} else {
  argsCond1 <- opt$condition1
  conditions_specified<- append(conditions_specified, argsCond1)
}

if (is.null(opt$condition2)) {
  argsCond2 <- FALSE 
} else {
  argsCond2 <- opt$condition2
  conditions_specified<- append(conditions_specified, argsCond2)
}
if (is.null(opt$condition3)) {
  argsCond3 <- FALSE 
} else {
  argsCond3 <- opt$condition3
  conditions_specified<- append(conditions_specified, argsCond3)
}

if (is.null(opt$condition4)) {
  argsCond4 <- FALSE 
} else {
  argsCond4 <- opt$condition4
  conditions_specified<- append(conditions_specified, argsCond4)
}
if (is.null(opt$condition5)) {
  argsCond5 <- FALSE 
} else {
  argsCond5 <- opt$condition5
  conditions_specified<- append(conditions_specified, argsCond5)
}

if (is.null(opt$condition6)) {
  argsCond6 <- FALSE 
} else {
  argsCond6 <- opt$condition6
  conditions_specified<- append(conditions_specified, argsCond6)
}
if (is.null(opt$condition7)) {
  argsCond7 <- FALSE 
} else {
  argsCond7 <- opt$condition7
  conditions_specified<- append(conditions_specified, argsCond7)
}

if (is.null(opt$condition8)) {
  argsCond8 <- FALSE 
} else {
  argsCond8 <- opt$condition8
  conditions_specified<- append(conditions_specified, argsCond8)
}
if (is.null(opt$condition9)) {
  argsCond9 <- FALSE 
} else {
  argsCond9 <- opt$condition9
  conditions_specified<- append(conditions_specified, argsCond9)
}
if (is.null(opt$condition0)) {
  argsCond0 <- FALSE 
} else {
  argsCond0 <- opt$condition0
  conditions_specified<- append(conditions_specified, argsCond0)
}


if( sum(argsTTcount,argsScount) != length(conditions_specified) )  {
  # throw an error
  stop("The number of specified tumour types ans stages does not match the specified condition files! \n OR you did not spesity the number of files using -t and -s flags. ")
} else {
  print(paste0(length(conditions_specified), " conditions were specified:"))
  print (conditions_specified)
  cat("\n")
}


if (is.null(opt$tumournormal)) {
  argsTN <- FALSE 
} else {
  argsTN <- opt$tumournormal #if this is the only argument supplied, then the provided list will be split into TP/NT (not full dataset)
}



getDataBarcodes <- function(argsProjectID, barcodeList, paired=FALSE){
  # NB even though we are searching by short barcodes, GDCquery will return long ones 
  query.exp <- GDCquery(project = argsProjectID,
                        legacy = TRUE,
                        data.category = "Gene expression",
                        data.type = "Gene expression quantification",
                        platform = "Illumina HiSeq", 
                        file.type = "results",
                        experimental.strategy = "RNA-Seq",
                        barcode = barcodeList)        
  
  if (paired == TRUE){     # tumour/normal case
    
    # select samples that are primary solid tumor
    dataSmTP <- TCGAquery_SampleTypes(query.exp$results[[1]]$cases,"TP")
    # select samples that are solid tissue normal
    dataSmNT <- TCGAquery_SampleTypes(query.exp$results[[1]]$cases,"NT")  
    output = (list(TP=dataSmTP, NT=dataSmNT))
    
  } else {     # it is not a paired case; return just tumour samples for the queries barcodes (subtype)
    
    # select samples that are primary solid tumor
    dataSmTP <- TCGAquery_SampleTypes(query.exp$results[[1]]$cases,"TP")
    output = dataSmTP
  }
  return (output)  
}    



# extract  tumour/normal

argsTN<-as.vector(read.table(argsTN, as.is = T, header = FALSE) )$V1
print (paste0("Querying long sample barcodes for Tumour/Normal."))
pairedSamples<-getDataBarcodes(argsProjectID, argsTN, paired=TRUE)
allTumourSamples <- pairedSamples$TP
allNormalSamples <- pairedSamples$NT
print (paste0("Normal samples count: ", length(allNormalSamples)))
tumourName <- "Tumour"
normalName <- "Normal"


# extract samples in each conditiona nd store thim in R list with condX as handles (same indexing can be used to extact condi actual name)
if ( length(conditions_specified) != 0 ) {
  
  samplesToExtract<-list()
  for (i in 1:length(conditions_specified)) { 
    print (paste0("Querying long sample barcodes for Condition ", i, ": ", conditions_specified[i]))
    this_condition_Samples<-as.vector(read.table(conditions_specified[i], as.is = T, header = FALSE))$V1
    condSamples<-getDataBarcodes(argsProjectID, this_condition_Samples)
    samplesToExtract[[ substr(conditions_specified[i],1,nchar(conditions_specified[i])-4) ]] <- condSamples
  }
  # add also normal samples to conditions list
  #samplesToExtract[["normal"]] <- allNormalSamples
}

print (names(samplesToExtract))



print ("Extraction of selected groups is finished.")
cat("\n")
print ("################### Data summary: ####################")
print (paste0("Preprocessed data contains ",dim(dataSE)[1]," genes and ", dim(dataSE)[2], " samples. "))

for (i in 1:length(samplesToExtract)) { 
  print (paste0("There are ", length(samplesToExtract[[i]]), " samples in condition ", names(samplesToExtract)[i] ," (", conditions_specified[i], ")"))
}
print ("######################################################")
cat("\n")



get_IDs <- function(data) {
  IDs <- strsplit(c(colnames(data)), "-") #split by hyphen into a list
  IDs <- ldply(IDs, rbind) #make a matrix samples VS barcode bits
  colnames(IDs) <- c('project', 'tss','participant', 'sample', "portion", "plate", "center")
  cols <- c("project", "tss", "participant")
  IDs$patient <- apply(IDs[,cols],1,paste,collapse = "-" ) #take 'cols' columns and make a new one
  barcode <- colnames(data) #get original sample names from input data
  IDs <- cbind(IDs, barcode) #add them to matrix
  condition <- gsub("11+[[:alpha:]]", "normal", as.character(IDs$sample)) #replace barcode nomenclature 11 for tumour
  condition  <- gsub("01+[[:alpha:]]", "cancer", condition) # 01 for normal # [[:alpha]] matches any letter (not important)
  IDs$condition <- condition #add condition column in the matrix
  IDs$myorder  <- 1:nrow(IDs) 
  ##test<-IDs[sort.list(IDs[,3]), ] #sort by participant (to see pairs)
  return(IDs)
}



# keep only samples of the subtypes we're investigation in the data frame 


newdataSE<-dataSE[, c(unique(unlist(samplesToExtract, use.names=FALSE)), allNormalSamples)]
samplesMatrix <- get_IDs(newdataSE)
print (paste0("Currently looking at ",dim(samplesMatrix)[1], " samples."))


addCondition <- function(samplesMatrix, conditionsList, TTcount=argsTTcount, TScount=argsScount){
  
  all_samples <- samplesMatrix$barcode
  
  # this will work if both types and stages been specified (but this block only deals with typess)
  if (TTcount != 0) {
    tumourTypes <- vector(mode="character", length=0)
    TT_list <- conditionsList[1:TTcount]
    print(paste0("Found " , TTcount, " tumour subtypes:"))
    print(names(TT_list))
    for (i in 1:length(all_samples)){
      barcode<-all_samples[i]
      #now iterating over tumour types
      for (j in 1:length(TT_list)){
        current_type<-unlist(TT_list[j], use.names = FALSE) # get barcodes of the curr tumout type
        if (barcode %in% current_type){
          tumourTypes[i]<- names(TT_list[j])
          break
        }else{
          tumourTypes[i]<- "unknown"
        }
        
        
      }
    }
    print (paste0("Labelled this many tumour types not unknown: ",length(tumourTypes[tumourTypes != "unknown"])) )
  }
  
  cat("\n")
  
  # this will work if both types and stages been specified (but this block only deals with stages)
  if ( (TTcount != 0) & (TScount != 0) ){
    tumourStages <- vector(mode="character", length=0)
    TS_list <- tail(conditionsList, TScount) #assuming stages come after types
    print(paste0("Found " , TScount, " tumour stages:"))
    print(names(TS_list))
    for (i in 1:length(all_samples)){
      barcode<-all_samples[i]
      #now iterating over tumour types
      for (j in 1:length(TS_list)){
        current_stage<-unlist(TS_list[j], use.names = FALSE) # get barcodes of the curr tumout type
        if (barcode %in% current_stage){
          tumourStages[i]<- names(TS_list[j])
          break
        }else{
          tumourStages[i]<- "unknown"
        }
      }
    }  
    print (paste0("Labelled this many tumour not unknown: ",length(tumourStages[tumourStages != "unknown"])) )
  }
  
  # this will work if ONLY stages been specified
  if ((TTcount == 0) & (TScount != 0)){
    tumourStages <- vector(mode="character", length=0)
    TS_list <- conditionsList[1:TScount] #assuming stages come after types
    print(paste0("Found " , TScount, " tumour stages:"))
    print(names(TS_list))
    for (i in 1:length(all_samples)){
      barcode<-all_samples[i]
      #now iterating over tumour types
      for (j in 1:length(TS_list)){
        current_stage<-unlist(TS_list[j], use.names = FALSE) # get barcodes of the curr tumout type
        if (barcode %in% current_stage){
          tumourStages[i]<- names(TS_list[j])
          break
        }else{
          tumourStages[i]<- "unknown"
        }
      }
    }  
    print (paste0("Labelled this many tumour not unknown: ",length(tumourStages[tumourStages != "unknown"])) )
  }
  
  if (TTcount != 0) {
    samplesMatrix$tumourTypes <- tumourTypes }
  
  if (TScount != 0) {
    samplesMatrix$tumourStages <- tumourStages}
  
  return (samplesMatrix)   
}



samplesMatrix<- addCondition(samplesMatrix, samplesToExtract, TTcount=argsTTcount, TScount=argsScount)

cat("\n")
print ("Assigned the following tumour types")
print (unique(samplesMatrix$tumourTypes))
print ("Assigned the following tumour stages")
print (unique(samplesMatrix$tumourStages))
print ("Assigned the following conditions")
print (unique(samplesMatrix$condition))


#replace 'unknowm' in normal samples with NA
samplesMatrix<- within(samplesMatrix, tumourTypes[condition == 'normal'] <- NA)
samplesMatrix <- within(samplesMatrix, tumourStages[condition == 'normal'] <- NA)


cat("\n")
cat("\n")
print (paste0("The new Summarized Experiment has dimentions: ", dim(newdataSE)[1]," " ,dim(newdataSE)[2]))
print (paste0("The new samples matrix has dimentions: ", dim(samplesMatrix)[1]," " ,dim(samplesMatrix)[2] ))
cat("\n")

# saving a new SE (only samples in question)
save(newdataSE,file=paste0(dirname(argsData),"/",unlist(strsplit(basename(argsData),".", fixed = T))[1],"_updatedSE_allTypes_allStages.rda"))

# saving the 'my_IDs' equivant but with types 
save(samplesMatrix, file=paste0(dirname(argsData),"/",unlist(strsplit(basename(argsData),".", fixed = T))[1],"_sampleMatrix_allTypes_allStages.rda"))
print("Data saved.")


#reading in the data! 
# testing<-get(load(paste0(dirname(argsData),"/",unlist(strsplit(basename(argsData),".", fixed = T))[1],"_allTypes_allStages.rda")))
# print (dim(testing))











#print("Starting DEA ...")


###### DEA anaylsis from TCGAbiolinks ---- to be replaced!
#dataDEGs <- TCGAanalyze_DEA(mat1 = dataSE[,cond1Samples],
#                            mat2 = dataSE[,cond2Samples],
#                            Cond1type = argsCond1,
#                            Cond2type = argsCond2,
#                            fdr.cut = 0.01 ,
#                            logFC.cut = 1,
#                            method = "glmLRT")  

#print (dataDEGs)
