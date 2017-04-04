setwd("~/Bioinformatics MSc UCPH/0_MasterThesis/TCGAbiolinks/CBL_scripts/data")

input<-read.csv("817BRCAfrompaper.csv",header=T)
dim(input)

library(dplyr)

names(input)

data <- select(input, Case.ID, PAM50, mRNA)
head(data)


# I got
#Basal-like HER2-enriched     Luminal A     Luminal B        normal   Normal-like 
#97            58           230           121           112             8 

table(data$PAM50)
#Basal   Her2   LumA   LumB Normal-like 
#136     65    415    176     25 