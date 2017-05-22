#rm(list=ls(all=TRUE))
setwd("~/Bioinformatics MSc UCPH/0_MasterThesis/TCGAbiolinks/CBL_scripts/data/DEA/EnrichmentAnalysis/")

class<- "PAM50"
class<- "stages"
class<- "morph"
class<-"Group1"
class<-"Group2"


DEgenes<- get(load(paste0(class, "_DE_genes_numbers3.rda")))
autoALL_DEgenes<- get(load(paste0(class, "_DE_AUTO_genes_numbers3.rda")))
autoCORE_DEgenes<- get(load(paste0(class, "_DE_AUTOCORE_genes_numbers3.rda")))
autoTF_DEgenes<- get(load(paste0(class, "_DE_AUTOTF_genes_numbers3.rda")))

N <- 15804 # The total number of genes :17372 ##### THIS WILL CHANGE

# add columns for enrichment, p.val, p.val.adj in autophagy df
extra_colnames<-c("up_oddsratio", "down_oddsratio", "both_oddsratio",
                    "up_pval","down_pval","both_pval",
                    "up_pval.adj","down_pval.adj","both_pval.adj")


fisher_test<-function(DEgenes,DEgenes_group,N,K){
  #going to iterate every contrast, up, down, both 
  for (i in 1:dim(DEgenes)[1]){
    for (j in 1:3){
      n<-DEgenes[i,j]   # The list of interesting genes - DEGs in this comparison
      k<-DEgenes_group[i,j]  #  The number of  autophagy genes in in DEGs of comparison
      
      m <- matrix(c( N - K - n + k,  K - k,
                     n - k,     k ),
                  nrow=2, ncol=2, 
                  dimnames = list(c("not Auto", "Auto"),
                                  c("not DE", "DE")))
      
      print(rownames(DEgenes)[i])
      print (m)
      cat("\n")
      
      #Fisher's exact test:
      res <- fisher.test(x=m, alternative="two.sided")
      DEgenes_group[i,j+3]<-res$estimate
      DEgenes_group[i,j+6]<-res$p.value
    }
  }
   #Correct for multiple testing
    DEgenes_group$up_pval.adj <- p.adjust(DEgenes_group$up_pval, method="BH")
    DEgenes_group$down_pval.adj <- p.adjust(DEgenes_group$down_pval, method="BH")
    DEgenes_group$both_pval.adj <- p.adjust(DEgenes_group$both_pval, method="BH")
    
    return(DEgenes_group)
}


#ALL
autoALL_DEgenes<-cbind(autoALL_DEgenes,matrix(data=NA, nrow=nrow(autoALL_DEgenes), ncol=9))
names(autoALL_DEgenes)[4:12]<-extra_colnames
K <- 1090 # The number of gene belonging to a gene famility---- AUTOPHAGY ALL
autoALL_DEgenes<-fisher_test(DEgenes,autoALL_DEgenes,N,K)

#CORE
autoCORE_DEgenes<-cbind(autoCORE_DEgenes,matrix(data=NA, nrow=nrow(autoCORE_DEgenes), ncol=9))
names(autoCORE_DEgenes)[4:12]<-extra_colnames
K <- 155 # The number of gene belonging to a gene famility---- AUTOPHAGY CORE
autoCORE_DEgenes<-fisher_test(DEgenes,autoCORE_DEgenes,N,K)


#TF
autoTF_DEgenes<-cbind(autoTF_DEgenes,matrix(data=NA, nrow=nrow(autoTF_DEgenes), ncol=9))
names(autoTF_DEgenes)[4:12]<-extra_colnames
K <- 97 # The number of gene belonging to a gene famility---- AUTOPHAGY TF
autoTF_DEgenes<-fisher_test(DEgenes,autoTF_DEgenes,N,K)


# the resuts df
autoALL_DEgenes
autoCORE_DEgenes
autoTF_DEgenes

# saving data for easy Viewing 

save(autoALL_DEgenes, file=paste0("EA_", class, "_All_autophagy.rda"))
save(autoCORE_DEgenes, file=paste0("EA_", class, "_Core_autophagy.rda"))
save(autoTF_DEgenes, file=paste0("EA_", class, "_TF_autophagy.rda"))

all<-get(load(paste0("EA_", class, "_All_autophagy.rda")))
core<-get(load(paste0("EA_", class, "_Core_autophagy.rda")))
tf<-get(load(paste0("EA_", class, "_TF_autophagy.rda")))

View(all)
View(core)
View(tf)


stop()

# plots investigation

#specify data!
this_data<-autoALL_DEgenes
this_data<-autoCORE_DEgenes
this_data<-autoTF_DEgenes

par(mfrow=c(1,3))

volcanoUp<-function(this_data){

    #### UP ####
    with(this_data, plot(log2(up_oddsratio), -log10(up_pval), pch=20, cex =2,
                         main="Enrichment vs Significance \n volcano plot for upregulated genes" ))#, xlim=c(-3,4), ylim=c(0,5)) )
    
    # Add colored points: red if padj<0.05, orange of log2FC>1, green if both)
    with(subset(this_data, up_pval.adj<0.05 ), points(log2(up_oddsratio), -log10(up_pval), pch=20, cex =2, col="blue"))
    
    with(subset(this_data, abs(up_oddsratio)>=1), points(log2(up_oddsratio), -log10(up_pval), pch=20, cex =2, col="orange"))
    with(subset(this_data, up_pval.adj<0.05 & abs(up_oddsratio)>=1), points(log2(up_oddsratio), -log10(up_pval), pch=20, cex =2, col="red"))
    
    abline(h = -log10(0.02), col = "green3", lty = 2) # adj Pvalue
    abline(v = 0, col = "blue", lty = 2) #enriched is >0
    mtext("adj pval\n = 0.05", side = 2, at = -log10(0.01), cex = 0.6, line = 0.5, las = 1)
    mtext(c(paste("enriched and significant")), side = 3, at = c(2), cex = 0.6, line = 0.2)
    
    # Label points with the textxy function from the calibrate plot
    #library(calibrate)
    this_data$names<-rownames(this_data)
    with(subset(this_data, up_pval.adj<0.05 & abs(up_oddsratio)>=1), textxy(log2(up_oddsratio),-log10(up_pval), labs=names, cex=.5))
}    
volcanoDown<-function(this_data){
      
       ##### DOWN ####
      # Make a basic volcano plot
      with(this_data, plot(log2(down_oddsratio), -log10(down_pval), pch=20, cex =2,
                           main="Enrichment vs Significance \n volcano plot for downregulated genes" ))#, xlim=c(-3,4), ylim=c(0,5)) )
      
      # Add colored points: red if padj<0.05, orange of log2FC>1, green if both)
      with(subset(this_data, down_pval.adj<0.05 ), points(log2(down_oddsratio), -log10(down_pval), pch=20, cex =2, col="blue"))
      
      with(subset(this_data, abs(down_oddsratio)>=1), points(log2(down_oddsratio), -log10(down_pval), pch=20, cex =2, col="orange"))
      with(subset(this_data, down_pval.adj<0.05 & abs(down_oddsratio)>=1), points(log2(down_oddsratio), -log10(down_pval), pch=20, cex =2, col="red"))
      
      abline(h = -log10(0.02), col = "green3", lty = 2) # adj Pvalue
      abline(v = 0, col = "blue", lty = 2) #enriched is >0
      mtext("adj pval\n = 0.05", side = 2, at = -log10(0.01), cex = 0.6, line = 0.5, las = 1)
      mtext(c(paste("enriched and significant")), side = 3, at = c(2), cex = 0.6, line = 0.2)
      
      # Label points with the textxy function from the calibrate plot
      #library(calibrate)
      this_data$names<-rownames(this_data)
      with(subset(this_data, down_pval.adj<0.05 & abs(down_oddsratio)>=1), textxy(log2(down_oddsratio)-1,-log10(down_pval), labs=names, cex=.5))
}
volcanoBoth<-function(this_data){

    #### BOTH ####
    with(this_data, plot(log2(both_oddsratio), -log10(both_pval), pch=20, cex =2,
                         main="Enrichment vs Significance \n volcano plot for bothregulated genes" ))#, xlim=c(-3,4), ylim=c(0,5)) )
    
    # Add colored points: red if padj<0.05, orange of log2FC>1, green if both)
    with(subset(this_data, both_pval.adj<0.05 ), points(log2(both_oddsratio), -log10(both_pval), pch=20, cex =2, col="blue"))
    
    with(subset(this_data, abs(both_oddsratio)>=1), points(log2(both_oddsratio), -log10(both_pval), pch=20, cex =2, col="orange"))
    with(subset(this_data, both_pval.adj<0.05 & abs(both_oddsratio)>=1), points(log2(both_oddsratio), -log10(both_pval), pch=20, cex =2, col="red"))
    
    abline(h = -log10(0.02), col = "green3", lty = 2) # adj Pvalue
    abline(v = 0, col = "blue", lty = 2) #enriched is >0
    mtext("adj pval\n = 0.05", side = 2, at = -log10(0.01), cex = 0.6, line = 0.5, las = 1)
    mtext(c(paste("enriched and significant")), side = 3, at = c(1), cex = 0.6, line = 0.2)
    
    # Label points with the textxy function from the calibrate plot
    #library(calibrate)
    this_data$names<-rownames(this_data)
    with(subset(this_data, both_pval.adj<0.05 & abs(both_oddsratio)>=1), textxy(log2(both_oddsratio)-1,-log10(both_pval), labs=names, cex=.5))
    
}

volcanoUp(this_data)
volcanoDown(this_data)
volcanoBoth(this_data)

