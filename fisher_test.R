#rm(list=ls(all=TRUE))
setwd("~/Bioinformatics MSc UCPH/0_MasterThesis/TCGAbiolinks/CBL_scripts/data/DEA/testing")

DEgenes<- get(load("stages_DE_genes_numbers3.rda"))
autoALL_DEgenes<- get(load("stages_DE_AUTO_genes_numbers3.rda"))
autoCORE_DEgenes<- get(load("stages_DE_AUTOCORE_genes_numbers3.rda"))
autoTF_DEgenes<- get(load("stages_DE_AUTOTF_genes_numbers3.rda"))

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

stop()
View(autoALL_DEgenes)
View(autoCORE_DEgenes)
View(autoTF_DEgenes)
