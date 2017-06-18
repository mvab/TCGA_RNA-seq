setwd("~/Bioinformatics MSc UCPH/0_MasterThesis/TCGAbiolinks/CBL_scripts/data/mfuzz_for_EA/")
source("../../functions.R") 


# putting mfuzz cluster results into format for Enrichment analysis

#loading the data, e.g.
data<-get(load("all_samples_6clusters_genesforDEA_cluster_data.rda"))
data<-get(load("all_samples_genesforDEA_cluster_data.rda"))
data<-get(load("stages_genesforDEA_cluster_data.rda"))
data<-get(load("group1_genesforDEA_cluster_data.rda"))
data<-get(load("lumB_genesforDEA_cluster_data.rda"))



#explore dataset (run this loop)
gene_sum=0
for ( i in names(data)){
  print (paste( i, ":" , length(data[[i]]) ))
  gene_sum=gene_sum+length(data[[i]])
}
gene_sum
##

#create a df to hold DE genes counts 
genesDEdf<-data.frame(matrix( NA,  nrow = length(names(data)), ncol = 1))
rownames(genesDEdf)<-names(data)
colnames(genesDEdf)<-c("genes")
head(genesDEdf)


#create a df to hold DE AUTOPHAGY genes counts 
genesDE_AUTOdf<-data.frame(matrix( NA,  nrow = length(names(data)), ncol = 1))
rownames(genesDE_AUTOdf)<-names(data)
colnames(genesDE_AUTOdf)<-c("genes")
head(genesDE_AUTOdf)


#create a df to hold DE AUTOPHAGY CORE genes counts 
genesDE_AUTOCOREdf<-data.frame(matrix( NA,  nrow = length(names(data)), ncol = 1))
rownames(genesDE_AUTOCOREdf)<-names(data)
colnames(genesDE_AUTOCOREdf)<-c("genes")

head(genesDE_AUTOCOREdf)

#create a df to hold DE AUTOPHAGY CORE genes counts 
genesDE_AUTOTFdf<-data.frame(matrix( NA,  nrow = length(names(data)), ncol = 1))
rownames(genesDE_AUTOTFdf)<-names(data)
colnames(genesDE_AUTOTFdf)<-c("genes")
head(genesDE_AUTOTFdf)


all_genes <- c() # for reference, count how many genes were DE in all comparisons
autophagy_genes <- c()
autophagyCORE_genes <- c()
autophagyTF_genes <- c()


setwd("~/Bioinformatics MSc UCPH/0_MasterThesis/TCGAbiolinks/CBL_scripts/data")

# fill in df with numbers
for (i in rownames(genesDEdf)) {
  print(i)

  # get numbers for all DE
  genesDEdf[i,"genes"]<- length(data[[i]])
    #get numbers for autophagy DE only 
  genesDE_AUTOdf[i,"genes"]<- length(sharedWithAuto(data[[i]]))
  #get numbers for autophagy CORE DE only 
  genesDE_AUTOCOREdf[i,"genes"]<- length(sharedWithAutoCORE(data[[i]]))
  #get numbers for autophagy TF DE only 
  genesDE_AUTOTFdf[i,"genes"]<- length(sharedWithAutoTF(data[[i]]))

  
  all_genes <- unique(sort(append(all_genes, data[[i]])))
  autophagy_genes <- unique(sort(append(autophagy_genes, sharedWithAuto(data[[i]]))))
  autophagyCORE_genes <- unique(sort(append(autophagyCORE_genes, sharedWithAutoCORE(data[[i]]))))
  autophagyTF_genes <- unique(sort(append(autophagyTF_genes, sharedWithAutoTF(data[[i]]))))
  
}

length(all_genes)
length(autophagy_genes)
length(autophagyCORE_genes)
length(autophagyTF_genes)



setwd("~/Bioinformatics MSc UCPH/0_MasterThesis/TCGAbiolinks/CBL_scripts/data/mfuzz_for_EA/")
class<-"all_samples_6clusters_genesforDEA"
class<-"stages_genesforDEA" 
class<-"all_samples_genesforDEA" 
class<-"group1_genesforDEA" 
class<-"lumB_genesforDEA" 

#all_samples_genesforDEA
save(genesDEdf, file=paste0(class, "_DE_genes_numbers3.rda"))
save(genesDE_AUTOdf, file=paste0(class, "_DE_AUTO_genes_numbers3.rda"))
save(genesDE_AUTOCOREdf, file=paste0(class, "_DE_AUTOCORE_genes_numbers3.rda"))
save(genesDE_AUTOTFdf, file=paste0(class, "_DE_AUTOTF_genes_numbers3.rda"))



#### fisher test 
DEgenes<- get(load(paste0(class, "_DE_genes_numbers3.rda")))
autoALL_DEgenes<- get(load(paste0(class, "_DE_AUTO_genes_numbers3.rda")))
autoCORE_DEgenes<- get(load(paste0(class, "_DE_AUTOCORE_genes_numbers3.rda")))
autoTF_DEgenes<- get(load(paste0(class, "_DE_AUTOTF_genes_numbers3.rda")))

N <- 13656 #10035#12684# 13656 # The total number of genes :17372 ##### THIS WILL CHANGE

# add columns for enrichment, p.val, p.val.adj in autophagy df
extra_colnames<-c("oddsratio", 
                  "pval",
                  "pval.adj")


fisher_test<-function(DEgenes,DEgenes_group,N,K){
  #going to iterate every contrast, up, down, both 
  for (i in rownames(DEgenes)){
      n<-DEgenes[i,]   # The list of interesting genes - total no of genes in a cluster
      k<-DEgenes_group[i,"genes"]  #  The number of  autophagy genes in a cluster
      
      m <- matrix(c( N - K - n + k,  K - k,
                     n - k,     k ),
                  nrow=2, ncol=2, 
                  dimnames = list(c("not Auto", "Auto"),
                                  c("not in cluster", "in this cluster")))
      
      print(i)
      print (m)
      cat("\n")
      
      #Fisher's exact test:
      res <- fisher.test(x=m, alternative="two.sided")
      DEgenes_group[i,"oddsratio"]<-res$estimate
      DEgenes_group[i,"pval"]<-res$p.value
    
  }
  #Correct for multiple testing
  DEgenes_group$pval.adj <- p.adjust(DEgenes_group$pval, method="BH")

  return(DEgenes_group)
}


#ALL
autoALL_DEgenes<-cbind(autoALL_DEgenes,matrix(data=NA, nrow=nrow(autoALL_DEgenes), ncol=3))
names(autoALL_DEgenes)[2:4]<-extra_colnames
K <- 1028#595#980#595#224#1090 #1035#198#535 #1001 # The number of gene belonging to a gene famility---- AUTOPHAGY ALL
autoALL_DEgenes<-fisher_test(DEgenes,autoALL_DEgenes,N,K);autoALL_DEgenes

#CORE
autoCORE_DEgenes<-cbind(autoCORE_DEgenes,matrix(data=NA, nrow=nrow(autoCORE_DEgenes), ncol=3))
names(autoCORE_DEgenes)[2:4]<-extra_colnames
K <-149 #68#142#39#68#20#155 #150#16#63# 145 # The number of gene belonging to a gene famility---- AUTOPHAGY CORE
autoCORE_DEgenes<-fisher_test(DEgenes,autoCORE_DEgenes,N,K);autoCORE_DEgenes


#TF
autoTF_DEgenes<-cbind(autoTF_DEgenes,matrix(data=NA, nrow=nrow(autoTF_DEgenes), ncol=3))
names(autoTF_DEgenes)[2:4]<-extra_colnames
K <- 909 #57#82##33#97 #89#27#47#85 # The number of gene belonging to a gene famility---- AUTOPHAGY TF
autoTF_DEgenes<-fisher_test(DEgenes,autoTF_DEgenes,N,K);autoTF_DEgenes

autoALL_DEgenes
autoCORE_DEgenes
autoTF_DEgenes

