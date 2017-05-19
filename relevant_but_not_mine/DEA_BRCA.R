#-------------------------------------------------------------------------

#DIFFERENTIAL EXPRESSION ANALYSIS - edgeR

#--------------------------------------------------------------------------
setwd("~/Bioinformatics MSc UCPH/0_MasterThesis/TCGAbiolinks/CBL_scripts")
library(edgeR)
library(plyr)
library(ggplot2)


#femaly only T/N data upload
dataSE<-get(load("brcaExp_PreprocessedData_wo_batch_updatedSE_allFemale.rda"))
samplesMatrix<-get(load("brcaExp_PreprocessedData_wo_batch_sampleMatrix_allFemale.rda"))


#function to remove unpaired samples
removeUnpaired <-  function (samples.matrix){

    # get barcodes of patients with paired samples
    paired <- samples.matrix[ duplicated(samples.matrix$participant, fromLast=FALSE) | duplicated(samples.matrix$participant, fromLast=TRUE),] #get patients with 2 entries
    pairedSorted <- paired[order(paired$participant), ]
    pairedParticipants <- as.character(paired$barcode)
    
    return (pairedParticipants)
}

pairedSamples<-removeUnpaired(samplesMatrix)
dataSEpaired <- dataSE[,pairedSamples]
dim(dataSEpaired)

samplesMatrixPaired <- samplesMatrix[samplesMatrix$barcode %in% pairedSamples, ]
#remove unsed levels (otherwise design matrix gets confused)
samplesMatrixPaired$barcode<-factor(samplesMatrixPaired$barcode)
samplesMatrixPaired$patient<-factor(samplesMatrixPaired$patient)
samplesMatrixPaired$participant<-factor(samplesMatrixPaired$participant)
samplesMatrixPaired$tss<-factor(samplesMatrixPaired$tss)


#check the number of normal and tumor samples TODO - has to be an error catcher
length(which(samplesMatrixPaired$condition=="normal"))
length(which(samplesMatrixPaired$condition=="cancer"))


# edgeR object
y <- DGEList(counts=dataSE) #by condition
w <- DGEList(counts=dataSEpaired) # condition + participant =paired


###########################                                         TODO: test the effect of this!
# Trim : Here we only keep genes with at least 2 counts in at least 4 samples:
above_one.Y <- rowSums(dataSE > 1)
trimmed_em.Y <- subset(dataSE, above_one.Y > 600) #half
above_one.W <- rowSums(dataSEpaired > 1)
trimmed_em.W <- subset(dataSEpaired, above_one.W > 111) #half

y_filt <- DGEList(counts=trimmed_em.Y) #by condition
w_filt <- DGEList(counts=trimmed_em.W) # condition + participant =paired
#############################


# model matrix info
y$samples$batch <-as.factor(samplesMatrix$tss)
y$samples$condition <- as.factor(samplesMatrix$condition)
y$samples$patientID <- as.factor(samplesMatrix$participant)


w$samples$batch <-as.factor(samplesMatrixPaired$tss)
w$samples$condition <- as.factor(samplesMatrixPaired$condition)
w$samples$patientID <- as.factor(samplesMatrixPaired$participant)

y <- calcNormFactors(y)
w <- calcNormFactors(w)

#plotMDS(y) #(similar to pca)
#plotMDS(w)



# relevel data
y$samples$condition = relevel(y$samples$condition, ref="normal")
w$samples$condition = relevel(w$samples$condition, ref="normal")

#design matrix for conditions only
#design.mat.Y <- model.matrix(~condition, data=y$samples)
design.mat.Y <- model.matrix(~condition+batch, data=y$samples)

#design matrix for paired samples
#design.mat.W <- model.matrix(~condition+patientID, data=w$samples)
design.mat.W <- model.matrix(~condition+batch, data=w$samples)


# estimate dispersion
y <- estimateDisp(y,design.mat.Y)
w <- estimateDisp(w,design.mat.W)


plotBCV(y)
plotBCV(w)

# fit overdispersed poisson model
my.fit.Y <- glmFit(y, design.mat.Y)
my.fit.W <- glmFit(w, design.mat.W) 



# Performing likelihood ratio test
N_C_coef.Y <- glmLRT(my.fit.Y, coef = 2) #T/N only 

N_C_coef.W <- glmLRT(my.fit.W, coef = 2) # T/N + patient


# Filter for significance - set log fold change and fdr cutoff
filterSignificance <- function (lrtOut, setFDR=0.01, setlogFC=1){

    toptags<-topTags(lrtOut,n=Inf)
    tags<-toptags$table
    
    index.up <- which(tags$logFC >= setlogFC & tags$FDR < setFDR) #get upregulated genes index
    index.down <- which(tags$logFC <= -setlogFC & tags$FDR < setFDR) #get downregulated index

    print ( paste0( "Total number of genes:", dim(tags)[1]))
    print ( paste0( "Upregulated genes:", length(index.up)))
    print ( paste0( "Downregulated genes:", length(index.down)))
    
    # via index create aa vector up/down/noDE and add column to table
    direction <- c()
    direction[index.up] <- "up"
    direction[index.down] <- "down"
    direction[!(1:nrow(tags) %in% union(index.up,index.down))] <- "no DE"
    tags <- cbind(tags,direction)
    
    return(list(annotatedDEA = tags, toptags = toptags))
}

N_C_coef.Y_signif<-filterSignificance(N_C_coef.Y,setFDR=0.01, setlogFC=1)
Y_DEA <- N_C_coef.Y_signif$annotatedDEA
Y_topTags<-N_C_coef.Y_signif$toptags

N_C_coef.W_signif<-filterSignificance(N_C_coef.W,setFDR=0.01, setlogFC=1)
W_DEA <- N_C_coef.W_signif$annotatedDEA
W_topTags<-N_C_coef.W_signif$toptags


# Y : T/N  only 

### volcano plot

ggplot(Y_topTags$table, aes(x=logFC, y=-log10(PValue), col=p.adjust(PValue, method="BH")<0.01)) +
  geom_point(alpha=0.25) + 
  geom_vline(aes(xintercept=2), col="blue") + 
  geom_vline(aes(xintercept=-2), col="blue") + 
  scale_color_manual(values=c("black", "red")) + 
  theme(legend.position="none")

## pvalue distibution

qplot(data=Y_topTags$table, x=PValue, geom="histogram",
      color=I("black"), fill=I("hotpink"), binwidth=0.05)


# W : T/N paired

### volcano plot

ggplot(W_topTags$table, aes(x=logFC, y=-log10(PValue), col=p.adjust(PValue, method="BH")<0.01)) +
  geom_point(alpha=0.25) + 
  geom_vline(aes(xintercept=2), col="blue") + 
  geom_vline(aes(xintercept=-2), col="blue") + 
  scale_color_manual(values=c("black", "red")) + 
  theme(legend.position="none")

## pvalue distibution

qplot(data=W-topTags$table, x=PValue, geom="histogram",
      color=I("black"), fill=I("hotpink"), binwidth=0.05)






###### OUTPUT 



# differentially expressed genes file
output.dir <- "~/Bioinformatics MSc UCPH/0_MasterThesis/TCGAbiolinks/CBL_scripts/output_files/"
write.csv(Y_DEA, paste0(output.dir,"edgeR_conditionOnly_BRCA.csv"), quote = FALSE)
write.csv(W_DEA, paste0(output.dir,"edgeR_paired_BRCA.csv"), quote = FALSE)

# up and down regulated genes
up.Y <- data.frame(rownames(Y_DEA[Y_DEA$direction == "up", ]))
down.Y <- data.frame(rownames(Y_DEA[Y_DEA$direction == "down", ]))
up.W <- data.frame(rownames(W_DEA[W_DEA$direction == "up", ]))
down.W <- data.frame(rownames(W_DEA[W_DEA$direction == "down", ]))

write.table(up.Y, paste0(output.dir,"up_edgeR_contionOnly_BRCA.txt"), sep = "\t",col.names = FALSE,row.names = FALSE, quote = FALSE)
write.table(down.Y,paste0(output.dir,"down_edgeR_contionOnly_BRCA.txt"), sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)

write.table(up.W, paste0(output.dir,"up_edgeR_paired_BRCA.txt"), sep = "\t",col.names = FALSE,row.names = FALSE, quote = FALSE)
write.table(down.W,paste0(output.dir,"down_edgeR_paired_BRCA.txt"), sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)


###
upBoth<-union(up.Y, up.W)
length(upBoth)

# venn?
# unioin stuff? 

m
