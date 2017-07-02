#! /usr/bin/env Rscript
rm(list=ls(all=TRUE))

#library(TCGAbiolinks)
library(SummarizedExperiment)
library(plyr)
library(edgeR)

library(ggplot2)
library(RColorBrewer)

library(limma)
#library(sva)

setwd("~/Bioinformatics MSc UCPH/0_MasterThesis/TCGAbiolinks/CBL_scripts/data")
source("../functions.R")


## most recent all types and stages cancer and normal
dataSE<-get(load("BRCA_Illumina_HiSeqnew_updatedSE_allFemale_PreprocessedData_wo_batch_GC_010_updatedSE_allTypes_allStages_Female.rda"))
samples.matrix<-get(load("BRCA_Illumina_HiSeqnew_updatedSE_allFemale_PreprocessedData_wo_batch_GC_010_sampleMatrix_allTypes_allStages_Female.rda"))
dim(samples.matrix)
dim(dataSE)


#trying cancer only with subtypes (1081)
#dataSE<-get(load("BRCA_Illumina_HiSeqnew_PreprocessedData_wo_batch_GC_010_updatedSE_allTypes_allStages_Female_cancerOnly_2017-04-02_16-22.rda"))
#samples.matrix<-get(load("BRCA_Illumina_HiSeqnew_PreprocessedData_wo_batch_GC_010_sampleMatrix_allTypes_allStages_Female_cancerOnly_2017-04-02_16-23.rda"))
#dim(samples.matrix)
#dim(dataSE)



#### rename morphology
samples.matrix<-renameMorph(samples.matrix)

samples.matrix<-addClinData(samples.matrix)

#### extract data only for autophagy genes
dataSE <- getAutophagyGenes(dataSE)


##########################


PAMNPout<-addXtraPAMandNormal(samples.matrix)# (807 +104 -39) + 112 = 984
samples.matrix<-PAMNPout$samples.matrix 
dim(samples.matrix)

sampleTokeepSE <- PAMNPout$samplesToKeep
dataSE<- dataSE[,colnames(dataSE) %in% sampleTokeepSE]
dim(dataSE) 

##### testing fro male samples D;


# exclude mrphology samples Ductal mixed with Others
removeUnknownOther <- function(dataSE, samples.matrix){
  samples.matrix[samples.matrix$tumourTypes=="Ductual mixed with others",]$barcode -> to_remove
  to_remove<- as.character(to_remove)
  
  samples.matrix<-samples.matrix[!samples.matrix$barcode %in% to_remove,]
  dim(samples.matrix)
  dataSE<-dataSE[, !colnames(dataSE) %in% to_remove ]
  dim(dataSE)
  
  return(list(dataSE=dataSE, samples.matrix=samples.matrix))
}
removedsamples<-removeUnknownOther(dataSE, samples.matrix)
dataSE<-removedsamples$dataSE
samples.matrix<-removedsamples$samples.matrix
dim(dataSE)
dim(samples.matrix)



#for toy example
keep_subtype<-samples.matrix[samples.matrix$PAM50=="Basal-like" | samples.matrix$PAM50=="Normal",]$barcode
dataSE<- dataSE[,colnames(dataSE) %in% keep_subtype]
dim(dataSE)
samples.matrix<-samples.matrix[samples.matrix$PAM50=="Basal-like" | samples.matrix$PAM50=="Normal",]
dim(samples.matrix)
normal <- PAMNPout$normal_barcodes



##### heatmaps #####

#### genes classes setup ###

g_functions <- get(load("autophagy_functions.rda"))       #### add comments!
rownames(g_functions)<-g_functions$genes 

ok_autop_genes<-rownames(dataSE)
ok_g_functions<-g_functions[c(ok_autop_genes),]

thirteen_cols<- c("#3579b4","#c8c049","#8996da","#ee5345",
                  "#c84297","#43ea3b","#50a376","#281340",
                  "#6e5c41","#fd0d32","#f19832","#94f5d2","#b1f555")
ok_g_functions$gene_functions=factor(ok_g_functions$gene_functions,
              levels=c("multifunction","lipid","phosphatidyl","endo_exosomes",    
              "transport","rabs","docking_fusion","mito","autoph_core",
              "transcr_factors","receptors_ligands", "mTOR_induction","lysosome"))
names(thirteen_cols)<-levels(ok_g_functions$gene_functions)
ok_g_functions$genes <-NULL

###### patient classes setup #####

#get only interesting cols
names(samples.matrix)
design<-subset(samples.matrix[,c('tumourTypes','tumourStages','PAM50', 'tss')])
row.names(design)<-samples.matrix$barcode
colnames(design)<-c("Morphology", "Stages", "PAM50", "tss")

#recorder annotations
design$PAM50 = factor(design$PAM50, levels = c("Basal-like", "HER2-enriched", "Luminal A", "Luminal B", "Normal-like", "Normal"))
design$Stages = factor(design$Stages, levels = c("stage1", "stage2", "stage3", "stage4", "unknown")) 
design$Morphology = factor(design$Morphology,  
              levels = c("Ductal and lobular mixed","Ductal carcinoma","Lobular carcinoma", 
                         "Metaplastic carcinoma" ,"Mucinous adenocarcinoma", "Normal",  "Other")) #other is a mix of few samples of different ones

design$tss = factor(design$tss)

# colours
PAM50Col <- c( "#E69F00", "#0072B2", "#F0E442","red3", "#CC79A7","#009E73")##D55E00->altred
tumourStgCol = c( "red","green","blue", "yellow2", "white")
tumourTypCol = c( "#F0E442", "#56B4E9", "#CC79A7", "#D55E00", "#0072B2", "#009E73", "#E69F00")
tssCol=c("#fc8785", "#3cb777", "#575a6c", "#c390ec", "#c5eebd",
                      "#4dc9d5", "#d1f6f2", "#e4f90d", "#5c9134", "#4159e2",
                      "#f3132d", "#d0abef", "#80eb6e", "#05dda1", "#af7fae",
                      "#349ad6", "#f2e5a9", "#43a405", "#f05017", "#f655f5",   
                      "#01e034", "#3f7f56", "#016ecd", "#e2523f", "#b270e4","#d7b51c")


names(PAM50Col)<-levels(design$PAM50)
names(tumourTypCol)<-levels(design$Morphology)
names(tumourStgCol)<-levels(design$Stages)
names(tssCol)<-levels(design$tss)

annColour <-list(
  PAM50=PAM50Col,
  Morphology=tumourTypCol,
  Stages=tumourStgCol,
  tss=tssCol
  #gene_functions=thirteen_cols
)
######




##############   for method!!  ####
names(samples.matrix)
design<-as.data.frame(samples.matrix$condition)
row.names(design)<-samples.matrix$barcode
colnames(design)<-c("condition")

#recorder annotations
design$condition = factor(design$condition, levels = c("cancer", "normal"))
# colours
tnPalette=c("deeppink2", "dodgerblue2")
names(tnPalette)<-levels(design$condition)
annColour <-list(
  condition=tnPalette
)
#######


# normalise with log
dge <- DGEList(dataSE)
dge <- calcNormFactors(object=dge, method="TMM")
EM <- cpm(x=dge, log=TRUE)



## selecting top 1000 genes:

# calc variance for each gene (log2 data), take top 1000 highest ->
#look at the stuff that is actually changing
#use only cancer samples for varince calculation
EMc <- EM[, !(colnames(EM) %in% normal)]

EM_var<- apply(EMc, MARGIN=1, var) #calc var across the row
length(EM_var)
top_genes<- sort(EM_var, decreasing = TRUE)[1:1000]
#barplot(top_genes)

# cancer and normal
EM2<- subset(EM[c(names(top_genes)),])
dim(EM2)

#only cancer
EM2c<- subset(EMc[c(names(top_genes)),])
dim(EM2c)

save (EM2, file="PAM50_EM_for_cluster_testing_1K_cancer.rda")



# check home any autophagy genes in top 1000 cancer_only and all
data_only_auto<- getAutophagyGenes(dataSE)
auto_genes<-rownames(data_only_auto)
length(auto_genes)
length(intersect(auto_genes, names(top_genes))) #30 (canOnly), 26 (all); shared 25
          
                                      #362 in 50% top genes
# pam50 autophagy test
all_auto<-as.vector(read.table("autopahagy_genes.txt", as.is = T, header = FALSE))$V1
pam50<-as.vector(read.table("pam50_list.txt", as.is = T, header = FALSE))$V1
length(intersect(pam50,all_auto)) #9
length(intersect(pam50,auto_genes))#9


##################


#library(scales)
#show_col(colorspace::heat_hcl(n=9))
#show_col(colorRampPalette(colorspace::heat_hcl(n=9))(20))
#colorRampPalette(brewer.pal(9,"Blues"))(20)


pheatmap::pheatmap(mat = as.matrix(EM2), color = colorRampPalette(brewer.pal(9,"YlGnBu"))(20),
#pheatmap::pheatmap(mat = as.matrix(EM2), color = rev(colorRampPalette(colorspace::terrain_hcl(n=9))(20)),
                   clustering_distance_rows = 'euclidean', 
                   clustering_distance_cols = 'euclidean',
                   #clustering_distance_rows = 'correlation',
                   #clustering_distance_cols = 'correlation',
                   #clustering_distance_rows = 'manhattan', 
                   #clustering_distance_cols = 'manhattan', 
                   #scale="row",
                   annotation_col=design,
                   annotation_colors = annColour,
                   #annotation_row=ok_g_functions,
                   cluster_cols = T, cluster_rows = T, 
                   show_rownames = F,show_colnames = F,
                   fontsize = 5)







 #use this not for pheatmap!
#dm<-dist(t(EM))
#rownames(dm)



#genes

stop()

###### CLUSTERING BY GENES


# Calculate eucliedan distances
dim(EM)
distmat <- dist((EM))


mycluster<- hclust(row_dist)
plot(mycluster)






dim(EM)

average.lin <- function(x) hclust(x, method="average") #average linkage

heatmap.2(EM, hclustfun=average.lin, #RowSideColors = c(rep("blue", 32),rep("yellow", 32), rep("magenta", 32)), 
          col=redgreen, trace="none",
          dendrogram = c("both"))#, main="Heatmap of gene expression")













### why it uses patient on both sides?






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
