source("TCGAbiolinks_functions.R")

my.dir <- "."
setwd(my.dir)

# -----------------------------------------------------------------------------------------------------------------------------
# TCGABiolinks - PREPARE DATA
# -----------------------------------------------------------------------------------------------------------------------------

my.pattern <- "BRCA"

cancer <- my.pattern
PlatformCancer <- "IlluminaHiSeq_RNASeqV2"
dataType <- "rsem.genes.results"
pathCancer <- paste0("data",my.pattern,"/")

datQuery <- TCGAquery(tumor = cancer, platform = PlatformCancer, level = "3")  
lsSample <- TCGAquery_samplesfilter(query = datQuery)

# get subtype information
dataSubt <- TCGAquery_subtype(tumor = cancer)

# Which samples are Primary Solid Tumor
dataSmTP <- TCGAquery_SampleTypes(barcode = lsSample$IlluminaHiSeq_RNASeqV2, typesample = "TP")

# Which samples are Solid Tissue Normal
dataSmTN <- TCGAquery_SampleTypes(barcode = lsSample$IlluminaHiSeq_RNASeqV2, typesample ="NT")

# get clinical data
dataClin <- TCGAquery_clinic(tumor = cancer, clinical_data_type = "clinical_patient") 

# download samples
TCGAdownload(data = datQuery, path = pathCancer, type = dataType, samples =c(dataSmTP,dataSmTN))

# prepare dataframe
dataAssy <- TCGAprepare(query = datQuery, dir = pathCancer, type = dataType, save = TRUE, summarizedExperiment = TRUE, samples = c(dataSmTP,dataSmTN), filename = paste0(cancer,"_",PlatformCancer,".rda")) 


#load("BRCA_IlluminaHiSeq_RNASeqV2_all.rda")

# build dataframe
dataPrep <- TCGAanalyze_Preprocessing(object = dataAssy, cor.cut = 0.6)  




# -----------------------------------------------------------------------------------------------------------------------------
# NORMALIZE DATA - CG-content and gene-length.
# -----------------------------------------------------------------------------------------------------------------------------


# normalize for CG content
dataNorm <- TCGAanalyze_Normalization(tabDF = dataPrep, geneInfo = geneInfo, method = "gcContent")

# normalize for gene length
dataNorm <- TCGAanalyze_Normalization(tabDF = dataNorm, geneInfo = geneInfo, method = "geneLength")

# filter out low counts
dataFilt <- TCGAanalyze_Filtering(tabDF = dataNorm, method = "quantile", qnt.cut =  0.25)  




# -----------------------------------------------------------------------------------------------------------------------------
# OBTAIN INFO ON BATCH AND IDs (AND SUBTYPES)
# -----------------------------------------------------------------------------------------------------------------------------

# get IDs, plate information etc.
my_IDs <- get_IDs(dataPrep)

# get information on subtype
subdiagnosis <- merge(my_IDs, dataSubt, by="patient", all.x=TRUE)
subdiagnosis <- subdiagnosis[order(subdiagnosis$myorder), ]


# -----------------------------------------------------------------------------------------------------------------------------
# FILTERING - NO SUBTYPE AND PAIRED CANCER SAMPLE
# -----------------------------------------------------------------------------------------------------------------------------

# get barcodes of samples with no subtype label (Luminal etc.)
no_diag <- subdiagnosis[is.na(subdiagnosis$PAM50.mRNA), ] 
no_diag <- no_diag$barcode

# get barcodes of patients with paired samples
paired <- my_IDs[ duplicated(my_IDs$participant, fromLast=FALSE) | duplicated(my_IDs$participant, fromLast=TRUE),] #get patients with 2 entries
paired <- as.character(paired$participant)
paired <- my_IDs[my_IDs$participant %in% paired & my_IDs$condition == "cancer", ]  #get cancer samples from paired data
paired <- as.character(paired$barcode)

# concatnate barcodes to remove 
remove <- unique(sort(c(as.character(no_diag), as.character(paired)))) #join NA and paired cancer samples

# remove samples with no subdtype label and paired tumor sample from dataframe, subtype vector and ID-information.
my_IDs <- my_IDs[!my_IDs$barcode %in% remove, ]    
subdiagnosis <- subdiagnosis[!subdiagnosis$barcode %in% remove, ]
dataFilt <- dataFilt[,!colnames(dataFilt) %in% remove]

#521/1205 participants



# -----------------------------------------------------------------------------------------------------------------------------
# VISUALIZATION - MDS-PLOTS
# -----------------------------------------------------------------------------------------------------------------------------

# MDS plots to get comparison
#myMDSplot(dataPrep, my_IDs$condition, my_IDs$condition)
#myMDSplot(dataNorm, my_IDs$condition, my_IDs$condition)
myMDSplot(dataFilt, my_IDs$condition, my_IDs$condition)

# make normal samples from pairs "normal" while retaining info on subtype.
subdiagnosis$true_status <- ifelse(subdiagnosis$condition == "normal", "normal", as.character(subdiagnosis$PAM50.mRNA))

# MDS plot
myMDSplot(dataFilt, subdiagnosis$true_status, subdiagnosis$condition)





# -----------------------------------------------------------------------------------------------------------------------------
# DIFFERENTIAL EXPRESSION ANALYSIS - EDGER
# -----------------------------------------------------------------------------------------------------------------------------


# edgeR object
y <- DGEList(counts=dataFilt)


# model matrix info
y$samples$batch <- as.factor(as.integer(as.factor(as.character(subdiagnosis$plate))))
y$samples$condition <- as.factor(subdiagnosis$condition)
y$samples$true_status <- as.factor(subdiagnosis$true_status)
y$samples$participant <- as.factor(as.integer(as.factor(subdiagnosis$participant))) 


# correct for library size
y <- calcNormFactors(y)

# relevel data
y$samples$true_status = relevel(y$samples$true_status, ref="normal")

# design matrix - IMPORTANT, model on subtype and batch
design.mat <- model.matrix(~true_status+batch, data=y$samples)

# estimate dispersion
y <- estimateDisp(y,design.mat)


# fit overdispersed poisson model
my.fit <- glmFit(y, design.mat)

# -----------------------------------------------------------------------------------------------------------------------------

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



# Get all genes
all_L <- unique(sort(c(as.character(N_NL_E[[1]]$up), as.character(N_NL_E[[2]]$down), as.character(N_BL_E[[1]]$up), as.character(N_BL_E[[2]]$down), as.character(N_A_E[[1]]$up), as.character(N_A_E[[2]]$down), as.character(N_B_E[[1]]$up), as.character(N_B_E[[2]]$down), as.character(N_H_E[[1]]$up), as.character(N_H_E[[2]]$down), as.character(NL_BL_E[[1]]$up), as.character(NL_BL_E[[2]]$down), as.character(NL_A_E[[1]]$up), as.character(NL_A_E[[2]]$down), as.character(NL_B_E[[1]]$up), as.character(NL_B_E[[2]]$down), as.character(NL_H_E[[1]]$up), as.character(NL_H_E[[2]]$down), as.character(BL_A_E[[1]]$up), as.character(BL_A_E[[2]]$down), as.character(BL_B_E[[1]]$up), as.character(BL_B_E[[2]]$down), as.character(BL_H_E[[1]]$up), as.character(BL_H_E[[2]]$down), as.character(A_B_E[[1]]$up), as.character(A_B_E[[2]]$down), as.character(A_H_E[[1]]$up), as.character(A_H_E[[2]]$down), as.character(B_H_E[[1]]$up), as.character(B_H_E[[2]]$down))))


# Get expression values
all_E <- dataFilt[rownames(dataFilt) %in% all_E, ]

write.table(all_E, paste0(my.dir, my.pattern, "_DE_EdgeR.txt"), sp = "\t", quote = FALSE)

# Scale expression values for plotting
all_log_E <- log2(all_E+1)
all_scaled_E <- scale(all_log_E, center = TRUE, scale = FALSE)


# Heatmap
my.col <- get_colors(subdiagnosis$true_status)
heatmap.plus(as.matrix(all_scaled_E), Rowv=NULL, hclustfun=function(d) hclust(d, method="ward.D2"), trace="none", labRow="", labCol=subdiagnosis$true_status, ColSideColors=my.col, margins = c(14,8), cexCol=0.4)





# -----------------------------------------------------------------------------------------------------------------------------
# DIFFERENTIAL EXPRESSION ANALYSIS - LIMMA
# -----------------------------------------------------------------------------------------------------------------------------


# Design matrix with true_status for batch correction
mod_design <- model.matrix(~as.factor(subdiagnosis$true_status))

# Voom correction
v <- voom(dataFilt, mod_design, plot=TRUE)

# Batch correction
batch_corr <- ComBat(v$E, as.factor(as.integer(as.factor(as.character(subdiagnosis$plate)))), mod_design, par.prior=TRUE,prior.plots=FALSE)

# Visualization with MDS-plot
myMDSplot(dataFilt, subdiagnosis$true_status, subdiagnosis$condition)
myMDSplot(batch_corr, subdiagnosis$true_status, subdiagnosis$condition)

# -----------------------------------------------------------------------------------------------------------------------------

# Design matrix with true_status for differential expression analysis
true_status_design <- model.matrix(~0+as.factor(subdiagnosis$true_status))
colnames(true_status_design) <- levels(as.factor(subdiagnosis$true_status))
colnames(true_status_design) <- c("Basal", "HER2", "LumA", "LumB", "normalLike", "normal")


# -----------------------------------------------------------------------------------------------------------------------------

# Making group contrasts 
N_BL_cotr <- makeContrasts("Basal-normal", levels= true_status_design)
N_NL_cotr <- makeContrasts("normalLike-Basal", levels= true_status_design)
N_A_cotr <- makeContrasts("LumA-normal", levels= true_status_design)
N_B_cotr <- makeContrasts("LumB-normal", levels= true_status_design)
N_H_cotr <- makeContrasts("HER2-normal", levels= true_status_design)
NL_BL_cotr <- makeContrasts("Basal-normalLike", levels= true_status_design)
NL_A_cotr <- makeContrasts("LumA-normalLike", levels= true_status_design)
NL_B_cotr <- makeContrasts("LumB-normalLike", levels= true_status_design)
NL_H_cotr <- makeContrasts("HER2-normalLike", levels= true_status_design)
BL_A_cotr <- makeContrasts("LumA-Basal", levels= true_status_design)
BL_B_cotr <- makeContrasts("LumB-Basal", levels= true_status_design)
BL_H_cotr <- makeContrasts("HER2-Basal", levels= true_status_design)
A_B_cotr <- makeContrasts("LumB-LumA", levels= true_status_design)
A_H_cotr <- makeContrasts("HER2-LumA", levels= true_status_design)
B_H_cotr <- makeContrasts("HER2-LumB", levels= true_status_design)

# -----------------------------------------------------------------------------------------------------------------------------

# Filter for significance - set log fold change and fdr cutoff
N_BL_L <- DE_limma(N_BL_cotr, batch_corr, true_status_design, 5, 0.001)
N_NL_L <- DE_limma(N_NL_cotr, batch_corr, true_status_design, 2, 0.001)
N_A_L <- DE_limma(N_A_cotr, batch_corr, true_status_design, 5, 0.001)
N_B_L <- DE_limma(N_B_cotr, batch_corr, true_status_design, 5, 0.001)
N_H_L <- DE_limma(N_H_cotr, batch_corr, true_status_design, 5, 0.001)
NL_BL_L <- DE_limma(NL_BL_cotr, batch_corr, true_status_design, 5, 0.001)
NL_A_L <- DE_limma(NL_A_cotr, batch_corr, true_status_design, 5, 0.001)
NL_B_L <- DE_limma(NL_B_cotr, batch_corr, true_status_design, 5, 0.001)
NL_H_L <- DE_limma(NL_H_cotr, batch_corr, true_status_design, 5, 0.001)
BL_A_L <- DE_limma(BL_A_cotr, batch_corr, true_status_design, 5, 0.001)
BL_B_L <- DE_limma(BL_B_cotr, batch_corr, true_status_design, 5, 0.001)
BL_H_L <- DE_limma(BL_H_cotr, batch_corr, true_status_design, 5, 0.001)
A_B_L <- DE_limma(A_B_cotr, batch_corr, true_status_design, 1, 0.05)
A_H_L <- DE_limma(A_H_cotr, batch_corr, true_status_design, 1, 0.05)
B_H_L <- DE_limma(B_H_cotr, batch_corr, true_status_design, 1, 0.05)



# -----------------------------------------------------------------------------------------------------------------------------
# VISUALIZATION OF LIMMA RESULTS
# -----------------------------------------------------------------------------------------------------------------------------

# Get all genes
all_L <- unique(sort(c(as.character(N_NL_L[[1]]$up), as.character(N_NL_L[[2]]$down), as.character(N_BL_L[[1]]$up), as.character(N_BL_L[[2]]$down), as.character(N_A_L[[1]]$up), as.character(N_A_L[[2]]$down), as.character(N_B_L[[1]]$up), as.character(N_B_L[[2]]$down), as.character(N_H_L[[1]]$up), as.character(N_H_L[[2]]$down), as.character(NL_BL_L[[1]]$up), as.character(NL_BL_L[[2]]$down), as.character(NL_A_L[[1]]$up), as.character(NL_A_L[[2]]$down), as.character(NL_B_L[[1]]$up), as.character(NL_B_L[[2]]$down), as.character(NL_H_L[[1]]$up), as.character(NL_H_L[[2]]$down), as.character(BL_A_L[[1]]$up), as.character(BL_A_L[[2]]$down), as.character(BL_B_L[[1]]$up), as.character(BL_B_L[[2]]$down), as.character(BL_H_L[[1]]$up), as.character(BL_H_L[[2]]$down), as.character(A_B_L[[1]]$up), as.character(A_B_L[[2]]$down), as.character(A_H_L[[1]]$up), as.character(A_H_L[[2]]$down), as.character(B_H_L[[1]]$up), as.character(B_H_L[[2]]$down))))

# Get expression values
all_L <- batch_corr[rownames(batch_corr) %in% all_L, ]

write.table(all_L, paste0(my.dir, my.pattern, "_DE_Limma.txt"), sp = "\t", quote = FALSE)

# Scale expression values for plotting
all_log_L <- log2(all_L+1)
all_scaled_L <- scale(all_L, center = TRUE, scale = FALSE)


# Heatmap
my.col <- get_colors(subdiagnosis$true_status)
heatmap.plus(as.matrix(all_scaled_L), Rowv=NULL, hclustfun=function(d) hclust(d, method="ward.D2"), trace="none", labRow="", labCol=subdiagnosis$true_status, ColSideColors=my.col, margins = c(14,8), cexCol=0.4)





# -----------------------------------------------------------------------------------------------------------------------------
# COMPARISON EDGER AND LIMMA
# -----------------------------------------------------------------------------------------------------------------------------
png(paste0(my.dir, my.pattern, "_EdgeR_Limma.png"), height= 400, width=600)

venn <- venn.diagram(list(A=rownames(all_E), B=rownames(all_L)), category.names = c("EdgeR", "Limma"), filename=NULL, lwd = 0.7, fill=rainbow(2), sub.cex = 2, cat.cex= 2, cex=1.5)
grid.draw(venn)

dev.off()

# -----------------------------------------------------------------------------------------------------------------------------
# OVERLAP EDGER AND LIMMA
# -----------------------------------------------------------------------------------------------------------------------------

overlap <- intersect(rownames(all_E), rownames(all_L))
overlap <- dataFilt[rownames(dataFilt) %in% overlap, ]

write.table(overlap, paste0(my.dir, my.pattern, "_overlap_EdgeR_Limma.txt"), sp = "\t", quote = FALSE)
