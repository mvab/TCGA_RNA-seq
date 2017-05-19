source("/data/user/marta/pipeline/DE/limma/TCGAbiolinks_functions.R")
library(TCGAbiolinks)
library(SummarizedExperiment)
library(limma)


#---------------------------------------------------------------------

#  LUAD - paired samples

#-------------------------------------------------------------------


SE_LUAD<-get(load("/data/user/marta/pipeline/DE/limma/LUAD/LUAD_Illumina_HiSeq_paired.rda"))
dataframe_LUAD<-get(load("/data/user/marta/pipeline/DE/limma/LUAD/LUAD_Illumina_HiSeq_PreprocessedData_paired.rda"))

length(which(colData(SE_LUAD)$shortLetterCode =="TP"))
length(which(colData(SE_LUAD)$shortLetterCode =="NT"))


my_IDs <- get_IDs(dataframe_LUAD)
#check the tumor samples number
length(which(my_IDs$condition=="cancer"))
#check the normal samples number
length(which(my_IDs$condition=="normal"))

#-------------------------------------------------------------------------

# DIFFERENTIAL EXPRESSION ANALYSIS- LIMMA

#-------------------------------------------------------------------------

dataframe_LUAD <- mean.duplicated.tumor(dataframe_LUAD,my_IDs)
my_IDs <- get_IDs(dataframe_LUAD)

#check the tumor and normal samples number after
#removing tumor technical replicates
length(which(my_IDs$condition=="cancer"))
length(which(my_IDs$condition=="normal"))

#Design matrix

condition <- as.factor(my_IDs$condition)
patientID <- as.factor(my_IDs$participant)

design.matrix <- model.matrix(~0+condition+patientID)
colnames(design.matrix)[c(1,2)] <- c("cancer","normal")


# Making group contrasts
N_C_cotr <- makeContrasts("cancer-normal", levels= design.matrix)


# Filter for significance - set log fold change and fdr cutoff
N_C_L <- DE_limma(N_C_cotr, dataframe_LUAD, design.matrix, 1, 0.01)

# differentially expressed genes file
write.csv(N_C_L, "limma_paired_LUAD.csv", quote = FALSE)

#number of up-regulated genes cancer vs normal
length(which(N_C_L$direction == "up"))

#number of down-regulated genes cancer vs normal
length(which(N_C_L$direction== "down"))

# total mRNA number
nrow(dataframe_LUAD)

# up and down regulated genes
up <- data.frame(rownames(N_C_L[N_C_L$direction == "up", ]))
down <- data.frame(rownames(N_C_L[N_C_L$direction == "down", ]))

write.table(up, "up_limma_paired_LUAD.txt", sep = "\t",col.names = FALSE,row.names = FALSE, quote = FALSE)
write.table(down, "down_limma_paired_LUAD.txt", sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)


#---------------------------------------------------------------

#                         VOLCANO PLOT

#---------------------------------------------------------------

dataDEG <- read.csv("limma_paired_LUAD.csv", sep=",", quote = "\n")

TCGAVisualize_volcano(dataDEG$logFC,dataDEG$adj.P.Val,
                      filename = "volcanoplot_limma_paired_LUAD.png",
                      x.cut = 5,y.cut = 10^-7,
                      names = rownames(dataDEG),
                      color = c("black","red","darkgreen"),
                      names.size = 2,
                      xlab = " Gene expression fold change (Log2)",
                      legend = "State",
                      title = "Volcano plot",width = 6,height = 4)

#-------------------------------------------------------------------------

#DIFFERENTIAL EXPRESSION ANALYSIS - edgeR

#--------------------------------------------------------------------------


SE1_LUAD <- get(load("/data/user/marta/pipeline/DE/limma/LUAD/LUAD_Illumina_HiSeq_paired_edger.rda"))
dataframe1_LUAD <- get(load("/data/user/marta/pipeline/DE/limma/LUAD/LUAD_Illumina_HiSeq_PreprocessedData_paired_edger.rda"))


my_IDs <- get_IDs(dataframe1_LUAD)
dataframe1_LUAD <- mean.duplicated.tumor(dataframe1_LUAD,my_IDs)
#update my_IDs dataframe
my_IDs <- get_IDs(dataframe1_LUAD)
#check the number of normal and tumor samples after tumor replicates removing
length(which(my_IDs$condition=="cancer"))
length(which(my_IDs$condition=="normal"))

# edgeR object
y <- DGEList(counts=dataframe1_LUAD)

# model matrix info
#y$samples$batch <-as.factor(my_IDs$plate)
y$samples$condition <- as.factor(my_IDs$condition)
y$samples$patientID <- as.factor(my_IDs$participant)

y <- calcNormFactors(y)

# relevel data
y$samples$condition = relevel(y$samples$condition, ref="normal")

#0)
design.mat <- model.matrix(~condition+patientID, data=y$samples)

#1)remove dependent column
#design.mat <- model.matrix(~condition+patientID+batch, data=y$samples)
#design.mat <- design.mat[,!(colnames(design.mat)=="batch1949")]
#length(which(colnames(design.mat)== "batch1949"))


# estimate dispersion
y <- estimateDisp(y,design.mat)

# fit overdispersed poisson model
my.fit <- glmFit(y, design.mat)

# Performing likelihood ratio test
N_C_coef <- glmLRT(my.fit, coef = 2)

# Filter for significance - set log fold change and fdr cutoff
N_C_E <- DE_edgeR(N_C_coef, y, 1, 0.01)

# differentially expressed genes file
output.dir <- "/data/user/marta/pipeline/DE/edgeR/LUAD/"
write.csv(N_C_E, paste0(output.dir,"edgeR_paired_LUAD.csv"), quote = FALSE)

#number of up-regulated genes cancer vs normal
length(which(N_C_E$direction == "up"))

#number of down-regulated genes cancer vs normal
length(which(N_C_E$direction== "down"))

# total mRNA number
nrow(dataframe1_LUAD)

# up and down regulated genes
up <- data.frame(rownames(N_C_E[N_C_E$direction == "up", ]))
down <- data.frame(rownames(N_C_E[N_C_E$direction == "down", ]))

write.table(up, paste0(output.dir,"up_edgeR_paired_LUAD.txt"), sep = "\t",col.names = FALSE,row.names = FALSE, quote = FALSE)
write.table(down,paste0(output.dir,"down_edgeR_paired_LUAD.txt"), sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)


