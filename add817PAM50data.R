rm(list=ls(all=TRUE))
library(dplyr)

setwd("~/Bioinformatics MSc UCPH/0_MasterThesis/TCGAbiolinks/CBL_scripts/data")

input<-read.csv("817BRCAfrompaper.csv",header=T)
dim(input)
names(input)

extradata <- select(input, Case.ID, PAM50)#, mRNA)
head(extradata)


# I got
#Basal-like HER2-enriched     Luminal A     Luminal B        normal   Normal-like 
#97            58           230           121           112             8 

table(extradata$PAM50)
#Basal   Her2   LumA   LumB Normal 
#136     65    415    176     25 

#rename for merging
colnames(extradata)<-c("patient", "PAM50new")

#cancer and normal samples
samplesMatrix<-get(load("BRCA_Illumina_HiSeqnew_updatedSE_allFemale_PreprocessedData_wo_batch_GC_010_sampleMatrix_allTypes_allStages_Female.rda"))
dim(samplesMatrix)

dataSubt<-get(load("dataSubt.rda"))

#putting in to a copy for temp convevience
samples.matrix<-samplesMatrix
#add extra patient information 
subdiagnosis <- merge(samples.matrix, dataSubt, by="patient", all.x=T) 
subdiagnosis <- subdiagnosis[order(subdiagnosis$myorder), ]
dim(subdiagnosis)
# make normal samples from pairs "normal" while retaining info on subtype.
subdiagnosis$PAM50 <- ifelse(subdiagnosis$condition == "normal", "normal", as.character(subdiagnosis$PAM50.mRNA))
samplesMatrixPAM<-subdiagnosis 
dim(samplesMatrixPAM)
table(samplesMatrixPAM$PAM50) #this contains all known labels but also PAM50


length(intersect(samplesMatrix$patient, extradata$patient))#807 in 1081

#removing 10 unidentified samples
x<-extradata$patient
y<-samplesMatrix$patient
odd10 <- x[which(!x %in% y)]
extradata <- extradata[!(extradata$patient %in% odd10), ]
dim(extradata)

#updating PAM50new to matching names
PAM50_updname <- vector(mode="character", length=0)
for (i in 1:length(extradata$PAM50new)) {
  
    if (extradata$PAM50new[i] =="LumA"){
      PAM50_updname[i] <- "Luminal A"
    }else if  (extradata$PAM50new[i] =="LumB"){
      PAM50_updname[i] <- "Luminal B" 
    }else if  (extradata$PAM50new[i] =="Her2"){
      PAM50_updname[i] <- "HER2-enriched"
    }else if (extradata$PAM50new[i] =="Basal"){
      PAM50_updname[i] <- "Basal-like"
    }else if(extradata$PAM50new[i] =="Normal"){
      PAM50_updname[i] <- "Normal-like"
    }else{
      print ("ERROR")
      print(extradata$PAM50new[i])
    }
}
extradata$PAM50new <- PAM50_updname
head(extradata)
dim(extradata)
table(extradata$PAM50new)


save(extradata, file ="extradata_var.rda")
#ready for merge!



#need a matrix with 'my' PAM50, but also all unassigned samples in it!
#merge PAM50new with original samples matrix
samplesMatrixUpd<- merge(extradata, samplesMatrixPAM, by="patient", all.x=T, all.y=T)
dim(samplesMatrixUpd)


# ad hoc comparison of new and existing PAM
comparePAM<-subset(samplesMatrixUpd[, c("patient", "PAM50","PAM50new")])
dim(comparePAM)
comparePAM_cancerOnly<-(comparePAM[which(comparePAM$PAM50 != "normal" | is.na(comparePAM$PAM50)),])
dim(comparePAM_cancerOnly)

#remove samples that are NA in both

comparePAM_cancerOnly_noNA<-data.frame()
for (i in 1:length(comparePAM_cancerOnly$patient)){
  
  if (is.na(comparePAM_cancerOnly$PAM50[i]) && is.na(comparePAM_cancerOnly$PAM50new[i])){
     print(paste0("This sample ", comparePAM_cancerOnly$patient[i], " is NA"))
     
  }else{
     comparePAM_cancerOnly_noNA<-rbind(comparePAM_cancerOnly_noNA, comparePAM_cancerOnly[i,])
  }
}
names(comparePAM_cancerOnly_noNA)<-names(comparePAM_cancerOnly)
dim(comparePAM_cancerOnly_noNA)
head(comparePAM_cancerOnly_noNA) ##### no double NA, no normal




#print cases when PAM 50 is mismatched
wrong_count=0
right_count=0
extraPAMnew=0
extraPAMold=0

PAM50_upd <- data.frame()
PAM50old<-data.frame()

for (i in 1:length(comparePAM_cancerOnly_noNA$patient)){
  
    if (is.na(comparePAM_cancerOnly_noNA$PAM50[i])){
      extraPAMnew=extraPAMnew+1
    }else if (is.na(comparePAM_cancerOnly_noNA$PAM50new[i])){  
      extraPAMold=extraPAMold+1
      PAM50old<-rbind(PAM50old, comparePAM_cancerOnly_noNA[i,])
      
    }else{
        if (comparePAM_cancerOnly_noNA$PAM50[i] != comparePAM_cancerOnly_noNA$PAM50new[i]){
            print (paste0(comparePAM_cancerOnly_noNA$PAM50[i], "---", comparePAM_cancerOnly_noNA$PAM50new[i]))
            wrong_count=wrong_count+1
            
            
            # for the wrong ones, save patient and newname
            PAM50_upd<- rbind(PAM50_upd, comparePAM_cancerOnly_noNA[i,c("patient", "PAM50new")])
            
        }else if(comparePAM_cancerOnly_noNA$PAM50[i] == comparePAM_cancerOnly_noNA$PAM50new[i]){
            #PAM50_upd[i]<- comparePAM_cancerOnly_noNA$PAM50new[i]
            right_count=right_count+1
        }
    }
}
right_count#371
wrong_count#39
extraPAMnew#397
extraPAMold#104

length(comparePAM_cancerOnly_noNA$PAM50new)#911

table(comparePAM_cancerOnly_noNA$PAM50new)


### new data + only in old data

dim(PAM50old)
PAM50old$PAM50new <-NULL
colnames(PAM50old)<- c("patient", "PAM50new") #this is actually old, but called new for rbind
names(PAM50old)

extradata_w_old<- rbind(extradata, PAM50old)
dim(extradata_w_old)
save(extradata_w_old, file ="extradata_var_w_old.rda")


####

dim(PAM50_upd)
save(PAM50_upd, file= "potentially_wrong_PAM.rda")
