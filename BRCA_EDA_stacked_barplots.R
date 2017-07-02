#! /usr/bin/env Rscript
rm(list=ls(all=TRUE))


library(SummarizedExperiment)

library(ggplot2)
library(RColorBrewer)
library(dplyr)



setwd("~/Bioinformatics MSc UCPH/0_MasterThesis/TCGAbiolinks/CBL_scripts/data")
source("../functions.R") 

samples.matrix<-get(load("BRCA_Illumina_HiSeqnew_updatedSE_allFemale_PreprocessedData_wo_batch_GC_010_sampleMatrix_allTypes_allStages_Female.rda"))
dim(samples.matrix)

#### rename morphology
samples.matrix<-renameMorph(samples.matrix)

####  add clinical data
samples.matrix<-addClinData(samples.matrix)


# only PAM50 samples 
PAMNPout<-addXtraPAMandNormal(samples.matrix)# (807 +104 -39) + 112 = 984
samples.matrix<-PAMNPout$samples.matrix 
dim(samples.matrix)

head(samples.matrix)

# exclude mrphology samples Ductal mixed with Others
samples.matrix<-samples.matrix[samples.matrix$tumourTypes!="Ductual mixed with others",]



#### exploring data with stacked barplots 
cbPalette <- c( "#E69F00", "#0072B2", "#F0E442","#D55E00",  "#009E73", "#CC79A7",   "#56B4E9", "black")


sum(DF2[DF2$year_diagnosed=='2010',]$n)

#####   TSS

set1<-subset(samples.matrix[ c("tss","PAM50")])
head(set1)

# Convert character vector to list of symbols
dots <- lapply(names(set1), as.symbol)
# Perform frequency counts
DF1 <- set1 %>%
  group_by_(.dots=dots) %>%
  summarise(n = n())

head(DF1)

ggplot(DF1, aes(x =tss, y = n, fill = PAM50)) + 
  geom_bar(stat = "identity")+
  scale_fill_manual(values=cbPalette)+
  ylab("Count") + xlab("Sample Source Site")+
  ggtitle("Stacked barplot showing count of each PAM50 subtype \n collected at each sample source site")+
  scale_y_continuous(expand = c(0, 0),breaks = seq(from = 0, by =25, to =250))+
  theme_classic()+
  theme(legend.position="bottom",
      plot.title = element_text(hjust = 0.5))


#### YEAR

set2<-subset(samples.matrix[ c("year_diagnosed","PAM50")])
head(set2)
# Convert character vector to list of symbols
dots <- lapply(names(set2), as.symbol)
# Perform frequency counts
DF2 <- set2 %>%
  group_by_(.dots=dots) %>%
  summarise(n = n())

head(DF2)

ggplot(DF2, aes(x =year_diagnosed, y = n, fill = PAM50)) + 
  geom_bar(stat = "identity")+
  scale_fill_manual(values=cbPalette)+
  theme_classic()+
  ylab("Count") + xlab("Year sample collected")+
  ggtitle("Stacked barplot showing count of each PAM50 subtype \n collected in different years")+
  scale_y_continuous(expand = c(0, 0),breaks = seq(from = 0, by =25, to =250))+
  theme_classic()+
  theme(legend.position="bottom",
        axis.text.x = element_text(angle = 90, hjust = 1),
        plot.title = element_text(hjust = 0.5))



### AGE

set3<-subset(samples.matrix[ c("ageGroups","PAM50")])
head(set3)
# Convert character vector to list of symbols
dots <- lapply(names(set3), as.symbol)
# Perform frequency counts
DF3 <- set3 %>%
  group_by_(.dots=dots) %>%
  summarise(n = n())

head(DF3)

ggplot(DF3, aes(x =ageGroups, y = n, fill = PAM50)) + 
  geom_bar(stat = "identity")+
  scale_fill_manual(values=cbPalette)+
  theme_classic()+
  ylab("Count") + xlab("Diagnosed at age")+
  ggtitle("Stacked barplot showing count of each \n PAM50 subtype  collected in patient age groups")+
  scale_y_continuous(expand = c(0, 0),breaks = seq(from = 0, by =25, to =300))+
  theme_classic()+
  theme(legend.position="bottom",
        axis.text.x = element_text(angle = 90, hjust = 1),
        plot.title = element_text(hjust = 0.5))



### STAGE


set4<-subset(samples.matrix[ c("tumourStages","PAM50")])
head(set4)
# Convert character vector to list of symbols
dots <- lapply(names(set4), as.symbol)
# Perform frequency counts
DF4 <- set4%>%
  group_by_(.dots=dots) %>%
  summarise(n = n())

head(DF4)

ggplot(DF4, aes(x =tumourStages, y = n, fill = PAM50)) + 
  geom_bar(stat = "identity")+
  scale_fill_manual(values=cbPalette)+
  theme_classic()+
  ylab("Count") + xlab("Stage")+
  ggtitle("Stacked barplot showing count of each \n PAM50 subtype in different cancer stages")+
  scale_y_continuous(expand = c(0, 0),breaks = seq(from = 0, by =50, to =1000))+
  theme(legend.position="bottom",
        axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(hjust = 0.5))

#or
c <- ggplot(DF4, aes(x = tumourStages, y = variable))
c + facet_wrap(~ value) + geom_bar(stat = "identity")



### MORPHOLOGY


set5<-subset(samples.matrix[ c("tumourTypes","PAM50")])
head(set5)
#### exploring data with stacked barplots 
# Convert character vector to list of symbols
dots <- lapply(names(set5), as.symbol)
# Perform frequency counts
DF5 <- set5%>%
  group_by_(.dots=dots) %>%
  summarise(n = n())

head(DF5)

ggplot(DF5, aes(x =tumourTypes, y = n, fill = PAM50)) + 
  geom_bar(stat = "identity")+
  scale_fill_manual(values=cbPalette)+
  theme_classic()+
  ylab("Count") + xlab("Morphology")+
  ggtitle("Stacked barplot showing count of each \n PAM50 subtype in different morphology groups")+
  scale_y_continuous(expand = c(0, 0),breaks = seq(from = 0, by =50, to =1000))+
  theme(legend.position="right",
        axis.text.x = element_text(angle = 35, hjust = 1),
        plot.title = element_text(hjust = 0.9))

#or
c <- ggplot(DF5, aes(x = tumourTypes, y = variable))
c + facet_wrap(~ value) + geom_bar(stat = "identity")


### testing year + tss

set5<-subset(samples.matrix[ c("year_diagnosed","tss")])
head(set5)
#### exploring data with stacked barplots 
# Convert character vector to list of symbols
dots <- lapply(names(set5), as.symbol)
# Perform frequency counts
DF5 <- set5%>%
  group_by_(.dots=dots) %>%
  summarise(n = n())

head(DF5)

fortyColours<-c("#3579b4","#c8c049","#8996da","#ee5345","#c84297","#43ea3b","#50a376","#281340","#6e5c41","#94f5d2","#fd0d32","#f19832","#b1f555","#d727b1","#f27456","#4bfe9f","#61789b","#2896be","#db1453","#c7a233","#d9a5c8","#1e785f","#3183e5","#82117f","#e5cbb0","#2dc194","#8f2ccf","#4e8fec","#e7ad8a","#234220","#4cee30","#d7b51c","#c96629","#472134","#36d1c8","#9f6f63","#ac8d3c","#a63dbd","#1db9d9","#10c399")


ggplot(DF5, aes(x =year_diagnosed, y = n, fill = tss)) + 
  geom_bar(stat = "identity")+
  scale_fill_manual(values=fortyColours)+
  theme_classic()+
  ylab("Count") + xlab("Year")+
  scale_y_continuous(expand = c(0, 0),breaks = seq(from = 0, by =25, to =1000))+
  theme(legend.position="right",
        axis.text.x = element_text(angle = 90, hjust = 1))

#or
c <- ggplot(DF5, aes(x = tumourTypes, y = variable))
c + facet_wrap(~ value) + geom_bar(stat = "identity")













### other MELT
names(samples.matrix)

set6<-subset(samples.matrix[ c("tumourStages","metastasis")])
head(set6)
#### exploring data with stacked barplots 
cbPalette <- c( "#E69F00", "#0072B2", "#F0E442","#D55E00",  "#009E73", "#CC79A7",   "#56B4E9", "black")


DF6 <- melt(set6, id.var="tumourStages")
head(DF6)

ggplot(DF6, aes(x =tumourStages, y = variable, fill = value)) + 
  geom_bar(stat = "identity")+
  scale_fill_manual(values=cbPalette)+
  theme_classic()+
  theme(legend.position="right",
        legend.title=element_blank(),
        axis.text.x = element_text(angle = 35, hjust = 1))

#
c + facet_wrap(~ value) + geom_bar(stat = "identity")
