#rm(list=ls(all=TRUE))
setwd("~/Bioinformatics MSc UCPH/0_MasterThesis/TCGAbiolinks/CBL_scripts/TCGA-BRCA_analysis_methods/data/EnrichmentAnalysisData/")


# 1)   select the model that you tested in DEA, and have the results for (run one!)

class<- "PAM50"
class<- "stages"
class<- "morph"
class<-"Group1" #PAM50 and stage
class<-"Group2" #morphology and stage
class<-"Group3" #PAM50 and morphology 

# load that model (run all!)
DEgenes<- get(load(paste0(class, "_DE_genes_numbers.rda")))
autoALL_DEgenes<- get(load(paste0(class, "_DE_AUTO_genes_numbers.rda")))
autoCORE_DEgenes<- get(load(paste0(class, "_DE_AUTOCORE_genes_numbers.rda")))
autoTF_DEgenes<- get(load(paste0(class, "_DE_AUTOTF_genes_numbers.rda")))


#### 2)   SETTING UP FOR FISHER TEST (enrichment)



# defining the main function for fisher test
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



## 3)  creating addition columns in the input files to hold fisher test results, running fisher test


# add columns for enrichment (oddsration), p.val, p.val.adj in autophagy df (fisher_test function calculates that for every contast)
extra_colnames<-c("up_oddsratio", "down_oddsratio", "both_oddsratio",
                  "up_pval",      "down_pval",      "both_pval",
                  "up_pval.adj",  "down_pval.adj",  "both_pval.adj")

N <- 15784 # The total number of genes (this number comes from cpm filtering in DEA >2>=19)

# running fisher test on all genes sets (all, core, TFs only)

#ALL
autoALL_DEgenes<-cbind(autoALL_DEgenes,matrix(data=NA, nrow=nrow(autoALL_DEgenes), ncol=9))
names(autoALL_DEgenes)[4:12]<-extra_colnames
K <- 1090 # The number of gene belonging to a gene famility---- AUTOPHAGY ALL (in those 15784 genes)
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

# the resuts df:
#autoALL_DEgenes
#autoCORE_DEgenes
#autoTF_DEgenes



#  4)    saving data for easy Viewing 

save(autoALL_DEgenes, file=paste0("EA_", class, "_All_autophagy.rda"))
save(autoCORE_DEgenes, file=paste0("EA_", class, "_Core_autophagy.rda"))
save(autoTF_DEgenes, file=paste0("EA_", class, "_TF_autophagy.rda"))

all<-get(load(paste0("EA_", class, "_All_autophagy.rda")))
core<-get(load(paste0("EA_", class, "_Core_autophagy.rda")))
tf<-get(load(paste0("EA_", class, "_TF_autophagy.rda")))

View(all)
View(core)
View(tf)


# 5) finding contrasts in which gene sets are significantly enriched (oddsratio >1 , p.adj/p.val < 0.05)

# define function hat will do it for the results df from fisher test
getEnrichedAndSignificant<- function(my_table, adj=TRUE){
  
    if (adj==TRUE){
  #by p.adj
  table_up <- my_table[,c(1,4,10)]
  table_down <- my_table[,c(2,5,11)]
  table_both <- my_table[,c(3,6,12)]
  out_up<- table_up[table_up$up_oddsratio >=1 & table_up$up_pval.adj <=0.05,]
  out_down<- table_down[table_down$down_oddsratio >=1 & table_down$down_pval.adj <=0.05,]
  out_both<- table_both[table_both$both_oddsratio >=1 & table_both$both_pval.adj <=0.05,]
  
  } else{  #adj ==FALSE
    
  #by just pval
  table_up <- my_table[,c(1,4,7)]
  table_down <- my_table[,c(2,5,8)]
  table_both <- my_table[,c(3,6,9)]
  out_up<- table_up[table_up$up_oddsratio >=1 & table_up$up_pval <=0.05,]
  out_down<- table_down[table_down$down_oddsratio >=1 & table_down$down_pval <=0.05,]
  out_both<- table_both[table_both$both_oddsratio >=1 & table_both$both_pval <=0.05,]  
  }
  
  
  return(list(up = out_up, down = out_down, both = out_both))
  
}


# get contasts signisicantly enriched at pvalue adjusted for multiple testing

#with adj p.val
print(getEnrichedAndSignificant(all))
print(getEnrichedAndSignificant(core))
print(getEnrichedAndSignificant(tf))

# <0 rows> means no contasts are significantly enriched


# get contasts signisicantly enriched at pvalue NOT adjusted for multiple testing

#with regulat p.val
print(getEnrichedAndSignificant(all, adj =FALSE))
print(getEnrichedAndSignificant(core,adj =FALSE))
print(getEnrichedAndSignificant(tf,adj =FALSE))


