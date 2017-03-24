############################################################



# -----------------------------------------------------------------------------------------------------------------------------
# DIFFERENTIAL EXPRESSION ANALYSIS - EDGER
# -----------------------------------------------------------------------------------------------------------------------------


# edgeR object
y <- DGEList(counts=newdataSE)


# model matrix info
y$samples$batch <- as.factor(as.integer(as.factor(as.character(samples.matrixPAM50_NP$plate))))
y$samples$condition <- as.factor(samples.matrixPAM50_NP$condition)
y$samples$PAM50 <- as.factor(samples.matrixPAM50_NP$PAM50)
y$samples$participant <- as.factor(as.integer(as.factor(samples.matrixPAM50_NP$participant))) 


# correct for library size
y <- calcNormFactors(y)
plotMDS(y)

# relevel data
y$samples$PAM50 = relevel(y$samples$PAM50, ref="normal")

# design matrix - IMPORTANT, model on subtype and batch
design.mat <- model.matrix(~PAM50+batch, data=y$samples)

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
my.col <- get_colors(subdiagnosis$PAM50)
heatmap.plus(as.matrix(all_scaled_E), Rowv=NULL, 
             hclustfun=function(d) hclust(d, method="ward.D2"), trace="none", labRow="",
             labCol=subdiagnosis$PAM50, ColSideColors=my.col, margins = c(14,8), cexCol=0.4)


