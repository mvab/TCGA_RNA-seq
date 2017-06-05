setwd("~/Bioinformatics MSc UCPH/0_MasterThesis/TCGAbiolinks/CBL_scripts/data")

# clinical data subset
clindata<-get(load("ClinicalData_needsfiltering.rda"))

dim(clindata) #1097
colnames(clindata)

#convert age to years
clindata$age_at_diag_years <- floor(clindata$age_at_diagnosis/365)

#add year when diagnosed
clindata$diag_year <- clindata$year_of_birth + clindata$age_at_diag_years

head(clindata)

clindata<-clindata[,c(3,2,5,6)]

to_replace <-c("not reported")
clindata$race[clindata$race %in% to_replace] <- NA

names(clindata)[3]<- "age_at_diagnosis"
names(clindata)[4]<- "year_diagnosed"

hist(clindata$age_at_diagnosis)
hist(clindata$year_diagnosed)
table(clindata$race)

barplot(height=table(clindata$race))

#replce NAs with 0 for converience
clindata[c("age_at_diagnosis", "year_diagnosed")][is.na(clindata[c("age_at_diagnosis", "year_diagnosed")])] <- 0



ageGroups <- vector(mode="character", length=0)
for (i in 1:length(clindata$age_at_diagnosis)){
      print(clindata$age_at_diagnosis[i])

      if (clindata$age_at_diagnosis[i] == 0){
        print ("error")
        ageGroups[i]<-  NA 
      } else if (clindata$age_at_diagnosis[i] < 40) {
          ageGroups[i]<- "below 40"
      } else if ((clindata$age_at_diagnosis[i] >= 40) & (clindata$age_at_diagnosis[i] < 50)){
          ageGroups[i]<- "40-49"
      } else if  ((clindata$age_at_diagnosis[i] >= 50) & (clindata$age_at_diagnosis[i] < 60)){
          ageGroups[i]<- "50-59"
      } else if ((clindata$age_at_diagnosis[i] >= 60) & (clindata$age_at_diagnosis[i] < 70)){
        ageGroups[i] <- "60-69"
      } else if ((clindata$age_at_diagnosis[i] >= 70) & (clindata$age_at_diagnosis[i] < 80)){
          ageGroups[i]<- "70-79"  
      } else if (clindata$age_at_diagnosis[i] >= 80) {
        ageGroups[i]<- "above 80"    
      }  
}

clindata$ageGroups<-ageGroups
length(ageGroups)


yearGroups <-  vector(mode="character", length=0)
for (i in 1:length(clindata$year_diagnosed)){
  #print(clindata$year_diagnosed[i])
  
  if (clindata$year_diagnosed[i] == 0){
    #print ("error")
    yearGroups[i] <-  NA 
  } else if (clindata$year_diagnosed[i] < 2000) {
    yearGroups[i]<- "before 2000"
  } else if ((clindata$year_diagnosed[i] >= 2000) & (clindata$year_diagnosed[i] < 2006)){
    yearGroups[i]<- "2000-2005"
  } else if  ((clindata$year_diagnosed[i] >= 2006) & (clindata$year_diagnosed[i] < 2011)){
    yearGroups[i]<- "2006-2010"
  } else if (clindata$year_diagnosed[i] >= 2011) {
    yearGroups[i]<- "2011-2016"    
  }
}
clindata$yearGroups<-yearGroups
length(yearGroups)

head(clindata)

save (clindata, file = "clinData_prepared.rda")
