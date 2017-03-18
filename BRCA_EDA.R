#! /usr/bin/env Rscript
rm(list=ls(all=TRUE))

library(TCGAbiolinks)
library(SummarizedExperiment)
library(plyr)
library(edgeR)

library(ggplot2)
library(RColorBrewer)

library(limma)
library(sva)

setwd("~/Bioinformatics MSc UCPH/0_MasterThesis/TCGAbiolinks/CBL_scripts")

## subtype data upload
dataSE<-get(load("brcaExp_PreprocessedData_wo_batch_updatedSE_subtypes.rda"))
samplesMatrix<-get(load("brcaExp_PreprocessedData_wo_batch_sampleMatrix_subtypes.rda"))


#femaly only T/N data upload
dataSE<-get(load("brcaExp_PreprocessedData_wo_batch_updatedSE_allFemale.rda"))
samplesMatrix<-get(load("brcaExp_PreprocessedData_wo_batch_sampleMatrix_allFemale.rda"))

########## exploration of different matrices : DO NOT RUN for standard case ########

addOnlyTopTSS <- function(samples.matrix){
  
  #filtering out patient with neither Positive nor Negative
  #tss_above50samples <-c("BH","A2","E2","A8","D8","E9", "AR", "B6","AC" ) #660 patients
  tss_2 <-c("A2","E9") 
  
  #samples.matrix <- samples.matrix[samples.matrix$tss %in% tss_above50samples,] 
  samples.matrix <- samples.matrix[samples.matrix$tss %in% tss_2,] 
  
  return (samples.matrix)
}   

  samples.matrixTSS<-addOnlyTopTSS(samplesMatrix)
  dim(samples.matrixTSS)
  # if used addOnlyTopTSS DO THIS:
  newdataSE<- dataSE[,c(samples.matrixTSS$barcode)]
  dim(newdataSE)


addBRCAReceptorStatus <- function(samples.matrix){

    # get information on subtype/pation details
    dataSubt <- TCGAquery_subtype(tumor = "BRCA") 
    
    #add extra patient information, but some have NAs; if FALSE- no NAs but lose ~200 patients #~735
    subdiagnosis <- merge(samples.matrix, dataSubt, by="patient", all.x=TRUE) 
    subdiagnosis <- subdiagnosis[order(subdiagnosis$myorder), ]
    
    # we only have female samples, so NA in gender can be replaced with FEMALE
    subdiagnosis[c("Gender")][is.na(subdiagnosis[c("Gender")])] <- "FEMALE"
    
    # note:  the are many thing not relevant - remove them and then do complete cases
    
    # testing how many pations how no NAs for ER PR HER cols
    complete_cases_receptors<- subdiagnosis[complete.cases(subdiagnosis[, c(1:12, 15:17)]),]# 704 patients 
    
    #filtering out patient with neither Positive nor Negative
    other_recp_status <-c("Indeterminate", "Not Performed", "Performed but Not Available", "Equivocal", "Not Available")
    complete_cases_receptors$ER.Status[complete_cases_receptors$ER.Status %in% other_recp_status] <- NA
    complete_cases_receptors$PR.Status[complete_cases_receptors$PR.Status %in% other_recp_status] <- NA
    complete_cases_receptors$HER2.Final.Status[complete_cases_receptors$HER2.Final.Status %in% other_recp_status] <- NA
    
    complete_cases_receptors<- complete_cases_receptors[complete.cases(complete_cases_receptors[, c(1:12, 15:17)]),]
    complete_cases_receptors<- subset(complete_cases_receptors[, c(1:12, 15:17)]) # 646 patients
    dim(complete_cases_receptors)
    
    
    # replacing Positive and Negative in receptor column with 0 and 1 
    
    complete_cases_receptors$ER.Status <- gsub('Positive', 1, complete_cases_receptors$ER.Status)
    complete_cases_receptors$ER.Status <- gsub('Negative', 0, complete_cases_receptors$ER.Status)
    complete_cases_receptors$PR.Status <- gsub('Positive', 1, complete_cases_receptors$PR.Status)
    complete_cases_receptors$PR.Status <- gsub('Negative', 0, complete_cases_receptors$PR.Status)
    complete_cases_receptors$HER2.Final.Status <- gsub('Positive', 1, complete_cases_receptors$HER2.Final.Status)
    complete_cases_receptors$HER2.Final.Status <- gsub('Negative', 0, complete_cases_receptors$HER2.Final.Status)
    
    receptor_statusNum = as.matrix(as.data.frame(lapply(subset(complete_cases_receptors[13:15]), as.numeric)))
    receptor_status$sum <- rowSums(receptor_statusNum)
    
    # adding subtype classification:
    
    ### Luminal     ###   ER+PR+HER+ or ER+PR+HER- , so 3 or 2
    ### HER2 over   ###   ER-PR-HER+  , so 1
    ### Triple Neg  ###   ER-PR-HER-  , so 0
    
    receptSubtype <- vector(mode="character", length=0)
    for (i in 1:length(receptor_status$sum)){
      if (receptor_status$sum[i] > 1) {
          receptSubtype[i]<-"Luminal"
      } else if (receptor_status$sum[i] == 1) {
          receptSubtype[i]<-"HER2"
      } else {
          receptSubtype[i]<-"TNR"
      }
    }
    receptor_status$receptSubtype <- receptSubtype
    # add receptSubtype to the main sample matrix and remove sepate receptors
    
    complete_cases_receptors<-subset(complete_cases_receptors[, c(1:12)])
    complete_cases_receptors$receptSubtype <- receptSubtype
    
    return (complete_cases_receptors)
}

  samples.matrixRecep<-addBRCAReceptorStatus(samplesMatrix)
  # if used addBRCAReceptorStatus DO THIS:
  newdataSE<- dataSE[,c(samples.matrixRecep$barcode)]
  dim(newdataSE)
  
  
addPAM50removePairedCancer<-function(samples.matrix){
  
  # get information on subtype/pation details
  dataSubt <- TCGAquery_subtype(tumor = "BRCA") 
  
  #add extra patient information, but some have NAs; if FALSE- no NAs but lose ~200 patients #~735
  subdiagnosis <- merge(samples.matrix, dataSubt, by="patient", all.x=TRUE) 
  subdiagnosis <- subdiagnosis[order(subdiagnosis$myorder), ]
  
  # get barcodes of samples with no subtype label (Luminal etc.)
  no_diag <- subdiagnosis[is.na(subdiagnosis$PAM50.mRNA), ] 
  no_diag <- no_diag$barcode
  
  # get barcodes of patients with paired samples
  paired <- samples.matrix[ duplicated(samples.matrix$participant, fromLast=FALSE) | duplicated(samples.matrix$participant, fromLast=TRUE),] #get patients with 2 entries
  paired <- as.character(paired$participant)
  paired <- samples.matrix[samples.matrix$participant %in% paired & samples.matrix$condition == "cancer", ]  #get cancer samples from paired data
  paired <- as.character(paired$barcode)
  
  # concatnate barcodes to remove 
  
  #join NA and paired cancer samples
  # <- unique(sort(c(as.character(no_diag), as.character(paired)))) 
  # OTHER OPTION dont remove paired
  to_remove <- unique(sort(as.character(no_diag)))
  print(length(to_remove))
                        
  # remove samples with no subdtype label and paired tumor sample from dataframe, subtype vector and ID-information.
  samples.matrix <- samples.matrix[!samples.matrix$barcode %in% to_remove, ]    
  subdiagnosis <- subdiagnosis[!subdiagnosis$barcode %in% to_remove, ]
  
  # make normal samples from pairs "normal" while retaining info on subtype.
  subdiagnosis$true_status <- ifelse(subdiagnosis$condition == "normal", "normal", as.character(subdiagnosis$PAM50.mRNA))
  
  return(list(samples.matrix=subdiagnosis, samplesToremove=to_remove))
}  

  PAMNPout<-addPAM50removePairedCancer(samplesMatrix)
  samples.matrixPAM50_NP<-PAMNPout$samples.matrix #515/1192
  sampleToremoveSE <- PAMNPout$samplesToremove
  
  # if used addPAM50removePairedCancer DO THIS:
  newdataSE<- dataSE[,!colnames(newdataSE) %in% sampleToremoveSE]
  dim(newdataSE) #7449x515 or 579 if include paired

  
  ######## part for seeing the difference between paired/nonpaired samples ######
  dataWpaired<-newdataSE
  dataWOpaired<-newdataSE
  dataW<-colnames(dataWpaired)
  dataWO<-colnames(dataWOpaired)
  
  r=0
  pairedS <- vector(mode="character", length=0)
  for (i in 1:length(dataW)){
    if (dataW[i] %in% dataWO){
      #print ("yes")
    } else{
      r=r+1
      print (dataW[i])
      pairedS[r]<-dataW[i]
    }}
  length(pairedS)
  cancersSE<-dataSE[, pairedS]
  cancerSEMat<-samples.matrixPAM50_NP[samples.matrixPAM50_NP$barcode %in% pairedS,]
  
  
  
  
  colnames(samples.matrixPAM50_NP)


############################################################
 

  
  # -----------------------------------------------------------------------------------------------------------------------------
  # DIFFERENTIAL EXPRESSION ANALYSIS - EDGER
  # -----------------------------------------------------------------------------------------------------------------------------
  
  
  # edgeR object
  y <- DGEList(counts=newdataSE)
  
  
  # model matrix info
  y$samples$batch <- as.factor(as.integer(as.factor(as.character(samples.matrixPAM50_NP$plate))))
  y$samples$condition <- as.factor(samples.matrixPAM50_NP$condition)
  y$samples$true_status <- as.factor(samples.matrixPAM50_NP$true_status)
  y$samples$participant <- as.factor(as.integer(as.factor(samples.matrixPAM50_NP$participant))) 
  
  
  # correct for library size
  y <- calcNormFactors(y)
  plotMDS(y)
  
  # relevel data
  y$samples$true_status = relevel(y$samples$true_status, ref="normal")
  
  # design matrix - IMPORTANT, model on subtype and batch
  design.mat <- model.matrix(~true_status+batch, data=y$samples)
  
  # estimate dispersion
  y <- estimateDisp(y,design.mat)
  
  
  # fit overdispersed poisson model
  my.fit <- glmFit(y, design.mat)
  
  # Performing likelihood ratio tests.
  N_BL_coef <- glmLRT(my.fit, coef = 2)
  N_H_coef <- glmLRT(my.fit, coef = 3)
  N_A_coef <- glmLRT(my.fit, coef = 4)
  N_B_coef <- glmLRT(my.fit, coef = 5)
  N_NL_coef <- glmLRT(my.fit, coef = 6)
  BL_H_coef <- glmLRT(my.fit, contrast = c(0,-1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0))
  BL_A_coef <- glmLRT(my.fit, contrast = c(0,-1,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0))
  BL_B_coef <- glmLRT(my.fit, contrast = c(0,-1,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0))
  BL_NL_coef <- glmLRT(my.fit, contrast = c(0,-1,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0))
  H_A_coef <- glmLRT(my.fit, contrast = c(0,0,-1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0))
  H_B_coef <- glmLRT(my.fit, contrast = c(0,0,-1,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0))
  H_NL_coef <- glmLRT(my.fit, contrast = c(0,0,-1,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0))
  A_B_coef <- glmLRT(my.fit, contrast = c(0,0,0,-1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0))
  A_NL_coef <- glmLRT(my.fit, contrast = c(0,0,0,-1,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0))
  B_NL_coef <- glmLRT(my.fit, contrast = c(0,0,0,0,-1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0))
  
  # -----------------------------------------------------------------------------------------------------------------------------
  DE_edgeR <- function(my.lrt, my.data, coLFC, coFDR) {
    my.tags <- topTags(my.lrt, n=nrow(my.data$counts))
    my.tags <- my.tags$table
    
    index.up <- which(my.tags$logFC >= coLFC & my.tags$FDR < coFDR)
    index.down <- which(my.tags$logFC <= -coLFC & my.tags$FDR < coFDR)
    direction <- c()
    direction[index.up] <- "up"
    direction[index.down] <- "down"
    direction[!(1:nrow(my.tags) %in% union(index.up,index.down))] <- "no DE"
    my.tags <- cbind(my.tags,direction)
    

    return(my.tags)
  }
  
  # Filter for significance - set log fold change and fdr cutoff
  N_BL_E <- DE_edgeR(N_BL_coef, y, 5, 0.001)
  N_H_E <- DE_edgeR(N_H_coef, y, 5, 0.001)
  N_A_E <- DE_edgeR(N_A_coef, y, 5, 0.001)
  N_B_E <- DE_edgeR(N_B_coef, y, 5, 0.001)
  N_NL_E <- DE_edgeR(N_NL_coef, y, 2, 0.001)
  BL_H_E <- DE_edgeR(BL_H_coef, y, 5, 0.001)
  BL_A_E <- DE_edgeR(BL_A_coef, y, 5, 0.001)
  BL_B_E <- DE_edgeR(BL_B_coef, y, 5, 0.001)
  BL_NL_E <- DE_edgeR(BL_NL_coef, y, 5, 0.001)
  H_A_E <- DE_edgeR(H_A_coef, y, 1, 0.05)
  H_B_E <- DE_edgeR(H_B_coef, y, 1, 0.05)
  H_NL_E <- DE_edgeR(H_NL_coef, y, 5, 0.001)
  A_B_E <- DE_edgeR(A_B_coef, y, 1, 0.05)
  A_NL_E <- DE_edgeR(A_NL_coef, y, 5, 0.001)
  B_NL_E <- DE_edgeR(B_NL_coef, y, 5, 0.001)
  
  
  
  # -----------------------------------------------------------------------------------------------------------------------------
  # VISUALIZATION OF EDGER RESULTS
  # -----------------------------------------------------------------------------------------------------------------------------
  
  
  
  # Get all genes #####THIS HAS TO BE FIXED!
  all_E <- unique(sort(c(as.character(N_NL_E[[1]]$up), as.character(N_NL_E[[2]]$down),
                         as.character(N_BL_E[[1]]$up), as.character(N_BL_E[[2]]$down),
                         as.character(N_A_E[[1]]$up), as.character(N_A_E[[2]]$down),
                         as.character(N_B_E[[1]]$up), as.character(N_B_E[[2]]$down),
                         as.character(N_H_E[[1]]$up), as.character(N_H_E[[2]]$down),
                         as.character(NL_BL_E[[1]]$up), as.character(NL_BL_E[[2]]$down),
                         as.character(NL_A_E[[1]]$up), as.character(NL_A_E[[2]]$down),
                         as.character(NL_B_E[[1]]$up), as.character(NL_B_E[[2]]$down),
                         as.character(NL_H_E[[1]]$up), as.character(NL_H_E[[2]]$down),
                         as.character(BL_A_E[[1]]$up), as.character(BL_A_E[[2]]$down), 
                         as.character(BL_B_E[[1]]$up), as.character(BL_B_E[[2]]$down), 
                         as.character(BL_H_E[[1]]$up), as.character(BL_H_E[[2]]$down),
                         as.character(A_B_E[[1]]$up), as.character(A_B_E[[2]]$down),
                         as.character(A_H_E[[1]]$up), as.character(A_H_E[[2]]$down),
                         as.character(B_H_E[[1]]$up), as.character(B_H_E[[2]]$down))))
  
  
  # Get expression values
  all_E <- samples.matrixPAM50_NP[rownames(samples.matrixPAM50_NP) %in% all_E, ]
  
  write.table(all_E, paste0("BRCA_DE_EdgeR.txt"), sp = "\t", quote = FALSE)
  
  # Scale expression values for plotting
  all_log_E <- log2(all_E+1)
  all_scaled_E <- scale(all_log_E, center = TRUE, scale = FALSE)
  
  
  # Heatmap
  my.col <- get_colors(subdiagnosis$true_status)
  heatmap.plus(as.matrix(all_scaled_E), Rowv=NULL, 
               hclustfun=function(d) hclust(d, method="ward.D2"), trace="none", labRow="",
               labCol=subdiagnosis$true_status, ColSideColors=my.col, margins = c(14,8), cexCol=0.4)
  
  
  
  

###################### EDA for T/N  #########################
dge <- DGEList(newdataSE)
dge <- calcNormFactors(object=dge, method="TMM")
EM <- cpm(x=dge, log=TRUE)

# Perfrom PCA
pca <- prcomp(x=t(EM), scale=TRUE, center=TRUE)

# Inspect components
summary(pca)
#type                               
qplot(data=as.data.frame(pca$x), x=PC1, y=PC4, geom=c("point"), color=samplesMatrix$condition)

#true diesease status
qplot(data=as.data.frame(pca$x), x=PC1, y=PC2, geom=c("point"), color=samples.matrixPAM50_NP$true_status)

ggplot(data=as.data.frame(pca$x),aes(x=PC1,y=PC2,col=samples.matrixPAM50_NP$true_status))+ #, shape =samples.matrixPAM50_NP$tumourType
  geom_point(size=2.5,alpha=0.9)+ #Size and alpha just for fun
  scale_colour_brewer(palette = "Set2")+ #your colors here
  theme_classic()+
  theme(legend.position="bottom")
#### exploring PCs

suppressMessages(library(magrittr))
suppressMessages(library(tidyr))
suppressMessages(library(dplyr))

# T/N
d <- data.frame(samplesMatrix, pca$x)
p <- d %>%  tbl_df %>%  gather(key="PC", value="score", contains("PC")) %>% filter(PC %in% paste0("PC", 1:9))

ggplot(p, aes(x=condition, y=score, fill=condition)) + 
  geom_boxplot(outlier.colour=NaN) + 
  #geom_jitter(alpha=0.5) + 
  facet_wrap(~PC) + 
  scale_fill_brewer(palette="Set2")

#PAM50
d <- data.frame(samples.matrixPAM50_NP, pca$x)
p <- d %>%  tbl_df %>%  gather(key="PC", value="score", contains("PC")) %>% filter(PC %in% paste0("PC", 1:9))

ggplot(p, aes(x=true_status, y=score, fill=true_status)) + 
  geom_boxplot(outlier.colour=NaN) + 
  #geom_jitter(alpha=0.5) + 
  facet_wrap(~PC) + 
  scale_fill_brewer(palette="Set2")


##### heatmaps

# Normalize without log
dgeH <- DGEList(newdataSE)
dgeH <- calcNormFactors(object=dgeH, method="TMM")
EM_nolog <- cpm(x=dgeH, log=FALSE)

# Calculate distance with manhattan
dm <- dist(t(scale(EM_nolog)), method="manhattan")

# Calculate eucliedan distances
distmat <- dist(t(EM))
cormat <- cor(EM)

corrplot::corrplot(corr=cormat, method="ellipse", diag=FALSE, order="AOE")

pheatmap::pheatmap(mat = as.matrix(dm), color = brewer.pal(name = "RdPu", n = 9),
                   clustering_distance_rows = dm, clustering_distance_cols = dm,
                   #annotation_col=design,
                   cluster_cols = TRUE, cluster_rows = FALSE, 
                   show_rownames = FALSE,show_colnames = FALSE)



#now plot using correlation instead of distance, and then annotation by both
#condition and PC-scores of the first 3 PCs. 
#Here we could also have highlighted the result of a clustering analysis.

#Let's first quickly rescale correlations to distances:
cor_as_dist <- as.dist((1-cormat) / 2) # rescale correlations to 01 range distances 


###### classes setup #####

#get only interesting cols
names(samples.matrixPAM50_NP)
design<-subset(samples.matrixPAM50_NP[,c('tumourType','true_status')])#, 'tss'
row.names(design)<-samples.matrixPAM50_NP$barcode

#recorder annotations
design$true_status = factor(design$true_status, levels = c("Basal-like", "HER2-enriched", "Luminal A", "Luminal B", "Normal-like", "normal"))
design$tumourType = factor(design$tumourType, levels = c("type1", "type2"))
true_statusCol= c( "turquoise4", "darkorange2", "slateblue2", "palevioletred2", "gold2", "olivedrab3")
tumourTypCol = c( "orange2","royalblue4")
names(true_statusCol)<-levels(design$true_status)
names(tumourTypCol)<-levels(design$tumourType)


annColour <-list(
  true_status=true_statusCol,
  tumourType=tumourTypCol
)
######


# heatmap of the correlation matrix 
pheatmap::pheatmap(mat = t(cormat), color = brewer.pal(name = "YlGnBu", n = 9), 
                   clustering_distance_rows = cor_as_dist, clustering_distance_cols = cor_as_dist, 
                   annotation_col=design,
                   annotation_colors = annColour,
                   cluster_cols = TRUE, cluster_rows = FALSE,
                   show_rownames = FALSE, show_colnames = FALSE)


pheatmap::pheatmap(mat=as.matrix(distmat), color=brewer.pal(name="YlOrRd", n=9), 
                   clustering_distance_rows=distmat, clustering_distance_cols=distmat,
                   annotation_row=design,
                   cluster_cols = FALSE, cluster_rows = TRUE, 
                   show_rownames = FALSE,show_colnames = FALSE)


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

################## EDA for subtypes ###################

dge <- DGEList(newdataSE)
dge <- calcNormFactors(object=dge, method="TMM")
EM <- cpm(x=dge, log=TRUE)


# Perfrom PCA
pca <- prcomp(x=t(EM), scale=TRUE, center=TRUE)

# Inspect components
summary(pca)




fortyColours<-c("#3579b4","#c8c049","#8996da","#ee5345","#c84297","#43ea3b","#50a376","#281340","#6e5c41","#94f5d2","#fd0d32","#f19832","#b1f555","#d727b1","#f27456","#4bfe9f","#61789b","#2896be","#db1453","#c7a233","#d9a5c8","#1e785f","#3183e5","#82117f","#e5cbb0","#2dc194","#8f2ccf","#4e8fec","#e7ad8a","#234220","#4cee30","#d7b51c","#c96629","#472134","#36d1c8","#9f6f63","#ac8d3c","#a63dbd","#1db9d9","#10c399")
thirteenColours<-c("#3579b4","#c8c049","#8996da","#ee5345","#c84297","#43ea3b","#50a376","#281340","#6e5c41","#fd0d32","#f19832","#94f5d2","#b1f555")
#TSS
qplot(data=as.data.frame(pca$x), x=PC1, y=PC2, geom=c("point"), color=out$tss)

qplot(data=as.data.frame(pca$x), x=PC1, y=PC2, geom=c("text"), color=out$tss, label = out$tss) +scale_color_manual(values=fortyColours)

ggplot(data=as.data.frame(pca$x),aes(x=PC1,y=PC2,col=out$tss))+
  geom_point(size=2,alpha=0.5)+ #Size and alpha just for fun
  scale_color_manual(values = fortyColours)+ #your colors here
  theme_classic()

#TSS after reduction 
qplot(data=as.data.frame(pca$x), x=PC1, y=PC3, geom=c("point"), color=samples.matrixTSS$tss)
ggplot(data=as.data.frame(pca$x),aes(x=PC1,y=PC5,col=samples.matrixTSS$tss))+
  geom_point(size=2,alpha=0.9)+ #Size and alpha just for fun
  #scale_color_manual(values = thirteenColours)+ #your colors here
  theme_classic()

#type                               #or 5
qplot(data=as.data.frame(pca$x), x=PC1, y=PC3, geom=c("point"), color=out$tumourType)

#portion (nothin)
qplot(data=as.data.frame(pca$x), x=PC3, y=PC2, geom=c("point"), color=out$portion)

# receptor
ggplot(data=as.data.frame(pca$x),aes(x=PC1,y=PC3,col=samples.matrixRecep$receptSubtype))+
  geom_point(size=2)+ #Size and alpha just for fun
  scale_fill_brewer(palette="Set3")+ #your colors here
  theme(legend.position="bottom")



ggplot(data=as.data.frame(pca$x),aes(x=PC3,y=PC7,col=out$portion))+
  geom_point(size=2,alpha=0.7)+ #Size and alpha just for fun
  scale_color_manual(values = thirteenColours)+ #your colors here
  theme_classic()

ggplot(data=as.data.frame(pca$x),aes(x=PC1,y=PC3,col=out$plate))+
  geom_point(size=2,alpha=0.7)+ #Size and alpha just for fun
  scale_color_manual(values = fortyColours)+ #your colors here
  theme_classic()



#### exploring PCs

suppressMessages(library(magrittr))
suppressMessages(library(tidyr))
suppressMessages(library(dplyr))

d <- data.frame(samples.matrixTSS, pca$x)
p <- d %>%  tbl_df %>%  gather(key="PC", value="score", contains("PC")) %>% filter(PC %in% paste0("PC", 1:9))

ggplot(p, aes(x=tumourType, y=score, fill=tumourType)) + 
  geom_boxplot(outlier.colour=NaN) + 
  #geom_jitter(alpha=0.5) + 
  facet_wrap(~PC) + 
  scale_fill_brewer(palette="Set2")

ggplot(p, aes(x=portion, y=score, fill=portion)) + 
  geom_boxplot(outlier.colour=NaN) + 
  #geom_jitter(alpha=0.05) +   #uncomment this to see most dense group
  facet_wrap(~PC) + 
  scale_color_manual(values = thirteenColours)
  #scale_fill_brewer(palette="Set3")

ggplot(p, aes(x=tss, y=score, fill=tss)) + 
  geom_boxplot(outlier.colour=NaN) + 
  #geom_jitter(alpha=0.05) +   #uncomment this to see most dense group
  facet_wrap(~PC) 

ggplot(p, aes(x=sample, y=score, fill=sample)) + 
  geom_boxplot(outlier.colour=NaN) + 
  #geom_jitter(alpha=0.5) + 
  facet_wrap(~PC) + 
  scale_fill_brewer(palette="Set3")

ggplot(p, aes(x=plate, y=score, fill=plate)) + 
  geom_boxplot(outlier.colour=NaN) + 
  #geom_jitter(alpha=0.5) + 
  facet_wrap(~PC) + 
  scale_fill_brewer(palette="Set3")

ggplot(p, aes(x=receptSubtype, y=score, fill=receptSubtype)) + 
  geom_boxplot(outlier.colour=NaN) + 
  #geom_jitter(alpha=0.5) + 
  facet_wrap(~PC) + 
  scale_fill_brewer(palette="Set3")

ggplot(p, aes(x=tss, y=score, fill=tss)) + 
  geom_boxplot(outlier.colour=NaN) + 
  #geom_jitter(alpha=0.5) + 
  facet_wrap(~PC) 
