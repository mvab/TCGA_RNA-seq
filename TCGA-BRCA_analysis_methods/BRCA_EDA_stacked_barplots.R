#! /usr/bin/env Rscript
rm(list=ls(all=TRUE))


library(SummarizedExperiment)
library(ggplot2)
library(RColorBrewer)
library(dplyr)



setwd("~/Bioinformatics MSc UCPH/0_MasterThesis/TCGAbiolinks/CBL_scripts/TCGA-BRCA_analysis_methods/data/")
source("../BRCA_functions.R") 


## most recent all types and stages cancer and normal
samples.matrix<-get(load("BRCA_Illumina_HiSeqnew_updatedSE_allFemale_PreprocessedData_wo_batch_GC_010_withSampleInfo_SM.rda"))
dim(samples.matrix)



########## PRE-ANALYIS STEPS ################

#### rename morphology
samples.matrix<-renameMorph(samples.matrix)


####  add clinical data
samples.matrix<-addClinData(samples.matrix)


#### add PAM50 annotation to samples.matrix and update data SE
PAMNPout<-addPAM50annotation(samples.matrix)

samples.matrix<-PAMNPout$samples.matrix #add annotation
dim(samples.matrix) 
head(samples.matrix)


# exclude mrphology samples Ductal mixed with Others
samples.matrix<-samples.matrix[samples.matrix$tumourTypes!="Ductual mixed with others",]


#### colur patette 
cbPalette <- c( "#E69F00", "#0072B2", "#F0E442","#D55E00",  "#009E73", "#CC79A7",   "#56B4E9", "black")



##### BARLOTS ######

# creating count barplots fo various sample groups, with proportions colourd by PAM50, always



#####   TSS (tumour source site)
set1<-subset(samples.matrix[ c("tss","PAM50")])
head(set1)

# Convert character vector to list of symbols
dots <- lapply(names(set1), as.symbol)
# Perform frequency counts
DF1 <- set1 %>%
  group_by_(.dots=dots) %>%
  summarise(n = n())

head(DF1)
dev.new()
ggplot(DF1, aes(x =tss, y = n, fill = PAM50)) + 
  geom_bar(stat = "identity")+
  scale_fill_manual(values=cbPalette)+
  ylab("Count") + xlab("Sample Source Site")+
  ggtitle("Stacked barplots of sample counts per source site")+
  scale_y_continuous(expand = c(0, 0),breaks = seq(from = 0, by =25, to =250))+
  theme_classic()+
  theme(legend.position="bottom",
      plot.title = element_text(hjust = 0.5),
      text=element_text(size=14,  family="serif"))


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
dev.new()
ggplot(DF2, aes(x =year_diagnosed, y = n, fill = PAM50)) + 
  geom_bar(stat = "identity")+
  scale_fill_manual(values=cbPalette)+
  theme_classic()+
  ylab("Count") + xlab("Year sample collected")+
  ggtitle("Stacked barplots of sample counts collected in different years")+
  scale_y_continuous(expand = c(0, 0),breaks = seq(from = 0, by =25, to =250))+
  theme_classic()+
  theme(legend.position="bottom",
        axis.text.x = element_text(angle = 90, hjust = 1),
        plot.title = element_text(hjust = 0.5),
        text=element_text(size=14,  family="serif"))



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
dev.new()
ggplot(DF4, aes(x =tumourStages, y = n, fill = PAM50)) + 
  geom_bar(stat = "identity")+
  scale_fill_manual(values=cbPalette)+
  theme_classic()+
  ylab("Count") + xlab("")+
  ggtitle("Stacked barplots of sample count per cancer stage")+
  scale_y_continuous(expand = c(0, 0),breaks = seq(from = 0, by =50, to =1000))+
  theme(legend.position="none",
        axis.text.x = element_text(angle = 45, hjust = 1, size=16),
        axis.text.y = element_text(size=16),
        plot.title = element_text(hjust = 0.5),
        text=element_text(size=16,  family="serif"))




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
dev.new()
ggplot(DF5, aes(x =tumourTypes, y = n, fill = PAM50)) + 
  geom_bar(stat = "identity")+
  scale_fill_manual(values=cbPalette)+
  theme_classic()+
  ylab("Count") + xlab("")+
  ggtitle("Stacked barplots of sample count per morphology types")+
  scale_y_continuous(expand = c(0, 0),breaks = seq(from = 0, by =50, to =1000))+
  theme(legend.position="right",
        axis.text.x = element_text(angle = 35, hjust = 1,size=16),
        axis.text.y = element_text(size=16),
        plot.title = element_text(hjust = 0.3),
        text=element_text(size=16,  family="serif"))




#### EXTRA PLOT


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






