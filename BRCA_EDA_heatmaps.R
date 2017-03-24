#! /usr/bin/env Rscript
rm(list=ls(all=TRUE))

library(TCGAbiolinks)
library(SummarizedExperiment)
library(plyr)
library(edgeR)

library(ggplot2)
library(RColorBrewer)

library(limma)
#library(sva)

setwd("~/Bioinformatics MSc UCPH/0_MasterThesis/TCGAbiolinks/CBL_scripts/data")


#femaly only T/N data upload
#dataSE<-get(load("brcaExp_PreprocessedData_wo_batch_updatedSE_allFemale.rda"))
#samplesMatrix<-get(load("brcaExp_PreprocessedData_wo_batch_sampleMatrix_allFemale.rda"))

## most recent all types and stages data upload
dataSE<-get(load("brcaExp_PreprocessedData_wo_batch_updatedSE_allTypes_allStages_Female.rda"))
samplesMatrix<-get(load("brcaExp_PreprocessedData_wo_batch_sampleMatrix_allTypes_allStages_Female.rda"))



addPAM50removePairedCancer<-function(samples.matrix){
  
  # get information on subtype/pation details
  dataSubt <- TCGAquery_subtype(tumor = "BRCA") 
  
  #add extra patient information, but some have NAs; if FALSE- no NAs but lose ~200 patients #~735
  subdiagnosis <- merge(samples.matrix, dataSubt, by="patient", all.x=TRUE) 
  subdiagnosis <- subdiagnosis[order(subdiagnosis$myorder), ]
  
  # get barcodes of samples with no subtype label (Luminal etc.)
  no_diag <- subdiagnosis[is.na(subdiagnosis$PAM50.mRNA), ] 
  no_diag <- no_diag$barcode
  print (length(no_diag))
  
  # get barcodes of patients with paired samples
  paired <- samples.matrix[ duplicated(samples.matrix$participant, fromLast=FALSE) | duplicated(samples.matrix$participant, fromLast=TRUE),] #get patients with 2 entries
  paired <- as.character(paired$participant)  #get patients with 2 entries
  
  #get cancer samples from paired data
  paired <- samples.matrix[samples.matrix$participant %in% paired & samples.matrix$condition == "cancer", ]  
  paired <- as.character(paired$barcode)
  
  # concatnate barcodes to remove 
  
  #join NA and paired cancer samples
  to_remove<- unique(sort(c(as.character(no_diag), as.character(paired)))) 
  # OTHER OPTION dont remove paired
  #to_remove <- unique(sort(as.character(no_diag)))
  print(length(to_remove))
  
  # remove samples with no subdtype label and paired tumor sample from dataframe, subtype vector and ID-information.
  samples.matrix <- samples.matrix[!samples.matrix$barcode %in% to_remove, ]    
  subdiagnosis <- subdiagnosis[!subdiagnosis$barcode %in% to_remove, ]
  
  # make normal samples from pairs "normal" while retaining info on subtype.
  subdiagnosis$PAM50 <- ifelse(subdiagnosis$condition == "normal", "normal", as.character(subdiagnosis$PAM50.mRNA))
  
  return(list(samples.matrix=subdiagnosis, samplesToremove=to_remove))
}  

PAMNPout<-addPAM50removePairedCancer(samplesMatrix)
samplesMatrix<-PAMNPout$samples.matrix #578/1191
sampleToremoveSE <- PAMNPout$samplesToremove

# if used addPAM50removePairedCancer DO THIS:
newdataSE<- dataSE[,!colnames(dataSE) %in% sampleToremoveSE]
dim(newdataSE) #7449x515 or 579 if include paired

table(samplesMatrix$PAM50)


## if without PAM50, then
#newdataSE<-dataSE

###
autophagy_genes<- as.vector(read.table("autopahagy_genes.txt", as.is = T, header = FALSE))$V1
dataSE_autophagy <-subset(newdataSE, row.names(newdataSE) %in% autophagy_genes)
dim(dataSE_autophagy)
####







##### heatmaps #####

###### classes setup #####

#get only interesting cols
names(samplesMatrix)
design<-subset(samplesMatrix[,c('tumourTypes','tumourStages','PAM50')])#, 'tss'
row.names(design)<-samplesMatrix$barcode

#recorder annotations
design$PAM50 = factor(design$PAM50, levels = c("Basal-like", "HER2-enriched", "Luminal A", "Luminal B", "Normal-like", "normal"))
design$tumourTypes = factor(design$tumourTypes, levels = c("female_84803", "female_85003", "female_85203", "female_85223","female_85233", "female_85753 ", "unknown")) #here unknown is mixed
design$tumourStages = factor(design$tumourStages, levels = c("stage1", "stage2", "stage3", "stage4", "unknown")) 


PAM50Col= c( "turquoise4", "darkorange2", "slateblue2", "palevioletred2", "gold2", "olivedrab3")
tumourTypCol = c( "#bebada", "#8dd3c7","#b3de69",  "#fb8072", "#80b1d3", "#fdb462","white")
#tumourTypCol = c(  "khaki2","palegreen1",  "tomato1","royalblue4","yellowgreen","deeppink3", "aquamarine1")
tumourStgCol = c( "coral1","coral4","salmon3", "sandybrown", "white")
names(PAM50Col)<-levels(design$PAM50)
names(tumourTypCol)<-levels(design$tumourTypes)
names(tumourStgCol)<-levels(design$tumourStages)

annColour <-list(
  PAM50=PAM50Col,
  tumourTypes=tumourTypCol,
  tumourStages=tumourStgCol
)
######



dge <- DGEList(newdataSE)
dge <- calcNormFactors(object=dge, method="TMM")
EM <- cpm(x=dge, log=TRUE)

# Normalize without log
dgeH <- DGEList(newdataSE)
dgeH <- calcNormFactors(object=dgeH, method="TMM")
EM_nolog <- cpm(x=dgeH, log=FALSE)

# Calculate distance with manhattan
dm <- dist(t(scale(EM_nolog)), method="manhattan")

#patients
pheatmap::pheatmap(mat = as.matrix(dm), color = brewer.pal(name = "RdPu", n = 9),
                   clustering_distance_rows = dm, clustering_distance_cols = dm,
                   annotation_col=design,
                   annotation_colors = annColour,
                   cluster_cols = TRUE, cluster_rows = FALSE, 
                   show_rownames = FALSE,show_colnames = FALSE)
#genes

# Normalize without log
dgeH_auto <- DGEList(dataSE_autophagy)
dgeH_auto <- calcNormFactors(object=dgeH_auto, method="TMM")
EM_nolog_auto <- cpm(x=dgeH_auto, log=FALSE)

# Calculate distance with manhattan
dm_auto <- dist(t(scale(EM_nolog_auto)), method="manhattan")
dim(dm_auto)

pheatmap::pheatmap(mat = as.matrix(dm_auto), color = brewer.pal(name = "RdPu", n = 9),
                   clustering_distance_rows = dm_auto, clustering_distance_cols = dm_auto,
                   annotation_col=design,
                   annotation_colors = annColour,
                   cluster_cols = T, cluster_rows = T, 
                   show_rownames = F,show_colnames = F)


### why it uses patient on both sides?



# Calculate eucliedan distances
distmat <- dist(t(EM))


#heatmap of distance matrix
pheatmap::pheatmap(mat=as.matrix(distmat), color=brewer.pal(name="YlOrRd", n=9), 
                   clustering_distance_rows=distmat, clustering_distance_cols=distmat,
                   annotation_row=design,
                   cluster_cols = FALSE, cluster_rows = TRUE, 
                   show_rownames = FALSE,show_colnames = FALSE)

#DO NOT RUN
cormat <- cor(EM)
#corrplot::corrplot(corr=cormat, method="ellipse", diag=FALSE, order="AOE")

#now plot using correlation instead of distance, and then annotation by both
#condition and PC-scores of the first 3 PCs. 
#Here we could also have highlighted the result of a clustering analysis.

#Let's first quickly rescale correlations to distances:
cor_as_dist <- as.dist((1-cormat) / 2) # rescale correlations to 01 range distances 


# heatmap of the correlation matrix by PATIENTS
pheatmap::pheatmap(mat = t(cormat), color = brewer.pal(name = "YlGnBu", n = 9), 
                   clustering_distance_rows = cor_as_dist, clustering_distance_cols = cor_as_dist, 
                   annotation_col=design,
                   annotation_colors = annColour,
                   cluster_cols = TRUE, cluster_rows = FALSE,
                   show_rownames = FALSE, show_colnames = FALSE)



#heatmap of correlation matrix by GENES
pheatmap::pheatmap(mat = t(cormat), color = brewer.pal(name = "YlGnBu", n = 9), 
                   clustering_distance_rows = cor_as_dist, clustering_distance_cols = cor_as_dist, 
                   annotation_col=design,
                   annotation_colors = annColour,
                   cluster_cols = FALSE, cluster_rows = TRUE,
                   show_rownames = FALSE, show_colnames = FALSE)



#Heatmap of large datasets

#To avoid overplotting, pheatmap has a feature for aggregating rows prior to 
#plotting using k-means, and then clustering the k centers instead.

#Using this feature we can create heatmap of even very big datasets.
#However, the resulting plot will depend heavily of the choice of k.
#When plotting actual expression values, it is often useful to scale the dendrogram, 
#so we do not primarily see differences in expression level, but rather differences
#in expression patterns.



EM.t <- cpm(x=dge, log=TRUE)

design2<-design
rownames(design2) <- colnames(EM)

pheatmap::pheatmap(mat = EM, kmeans_k = 50, annotation_col= design, annotation_colors = annColour,
                   scale = "column", gp = gpar(fill = "red"),
                   cluster_cols = TRUE, cluster_rows = TRUE, 
                   show_rownames = FALSE,show_colnames = FALSE)

EM.t<-t(EM)
dim(EM.t)
rownames()


pheatmap::pheatmap(mat = EM.t, kmeans_k = 100, #annotation_col= design2, 
                   scale = "row", gp = gpar(fill = "red"),
                   cluster_cols = TRUE, cluster_rows = TRUE, 
                   show_rownames = FALSE,show_colnames = TRUE)
