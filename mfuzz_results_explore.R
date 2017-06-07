

setwd("~/mfuzz_with_genedata/all_samples/")


#loading the data, e.g.
data<-get(load("autophagy_genes_all_samples_cluster_data.rda"))


#explore dataset (run this loop)
gene_sum=0
for ( i in names(data)){
  print (paste( i, ":" , length(data[[i]]) ))
  gene_sum=gene_sum+length(data[[i]])
}
#NB genes in each clister are non-overlapping

# total number of genes after filtering and clustering in this data group (std=0.3, membership  >0.6)
gene_sum

#get genes, e.g.
cluster1<-data[[1]]