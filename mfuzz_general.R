library(Mfuzz)
library(edgeR)


setwd("")


dataSE<- # has to be a count matrix N x M
samples.matrix<-  # has to be sample annotation matix M x k , where k is columns with annotations, one being time points
dim(samples.matrix)
dim(dataSE)


## get gene lengths for count data

# function 
getGeneLenght<-function(dataSE){
  
  
  #### some ad hoc  steps to remove genes that are mislabelled or not found in either dataSE or gene lengths df
  
  #removing rows where genes are not found
  remove<-c("100133144", "155060","HLA-DRB1","WBP1")
  dataSE<-dataSE[!rownames(dataSE) %in% remove, ]
  dim(dataSE)
  
  #rename rows to their other names
  rownames(dataSE)[match('8225', rownames(dataSE))]<-'GTPBP6'
  rownames(dataSE)[match('DHRS4L2', rownames(dataSE))]<-'DHRS4'
  rownames(dataSE)[match('90288', rownames(dataSE))]<-'C3orf25'
  rownames(dataSE)[match('SLC35E2.728661', rownames(dataSE))]<-'SLC35E2B'
  rownames(dataSE)[match('SLC35E2.9906', rownames(dataSE))]<-'SLC35E2'
  
  my_genes<-rownames(dataSE)
  #length(my_genes)
  
  # gene lengths df
  load("TCGAgeneExons.Rdata") #in /data
  gene_exons<-as.data.frame(geneExons)
  dim(gene_exons)
  colnames(gene_exons)
  
  ##  need to split gene SLC35E2 into 2 parts by location
  gene_exons[gene_exons$gene_name=='SLC35E2',]
  #SLC35E2B from  1590990  to 1624243
  #SLC35E2  from 1656279 to 1677431
  #so rowindex 198235:198244 has to be SLC35E2B in gene_name
  gene_exons[198235:198244,]$gene_name <- as.character('SLC35E2B')
  
  
  genes<-unique(gene_exons$gene_name)
  #length(genes) #25508
  
  #compare!
  shared<-intersect(genes,my_genes)
  #length(shared) #17368
  
  
  #get exon data only for my genes 
  my_gene_exons<-gene_exons[gene_exons$gene_name %in% shared, ]
  #length(unique(my_gene_exons$gene_name))
  
  # column 'width' is gene length
  # add gene length fro all matching genes
  gene_lengths<-aggregate(. ~ gene_name, data=my_gene_exons[,c("gene_name", "width")], FUN=sum)
  #dim(gene_lengths)
  colnames(gene_lengths)[2]<-c("gene_length")
  #head(gene_lengths)
  
  
  #reorder to match order in dataSE
  my_order<-cbind(rownames(dataSE), matrix(0,ncol=1,nrow=length(rownames(dataSE))))
  #head(my_order)
  
  for (i in 1:nrow(my_order)){
    g<-as.character(my_order[i,1])
    my_order[i,2]<-as.numeric(gene_lengths[gene_lengths$gene_name==g,]$gene_length)
  }
  gene_lengths_ordered<-as.data.frame(my_order)
  colnames(gene_lengths_ordered)<- c("gene_name", "gene_length")
  #head(gene_lengths_ordered)
  
  return(list(dataSE=dataSE, gene_lengths=gene_lengths_ordered))
}

#get gene lenghts for genes in the analysis 
getGeneLenght_out<-getGeneLenght(dataSE)
gene_lengths<-getGeneLenght_out$gene_lengths
dataSE<-getGeneLenght_out$dataSE #updated (some genes get removed)
dim(dataSE)
dim(gene_lengths)


# edgeR object
dge <- DGEList(counts=dataSE, genes = data.frame(Length= as.numeric(gene_lengths$gene_length)) ) 
#head(dge$genes)

### this step is need if you have a predifined list of genes you want to cluster
## load gene list that you want to be in the analysis
#genelist<-get(load("~/genes_for_DEA.rda"))
#length(genelist)
#dge<- dge[rownames(dge) %in% genelist, ]
#dim(dge) 


#### RPKM
yr <- calcNormFactors(dge) #default is TMM
# rpkm will use the normalized effective library sizes to compute rpkm instead of the raw library sizes. 

RPKM <- rpkm(yr, log=TRUE, prior.count = 1)
dim(RPKM) 




# remove samples whose stage is unknowm
known_stage<-samples.matrix[samples.matrix$tumourStages!='unknown',]$barcode    #### modify here to suit your case
RPKM<- RPKM[,colnames(RPKM) %in% known_stage]
dim(RPKM)


# create a rpkm matrix that has a column per 'time point', i.e. stage

#function
getAverageExpSampleSet<- function(dataSE, set){
  
  #get list of samples for this set
  set_samples<-samples.matrix[samples.matrix$tumourStages==set,]$barcode
  
  #get their exp vals
  dataSE<- dataSE[,colnames(dataSE) %in% set_samples]
  print(paste0(set, " has ", dim(dataSE)[2], " samples."))
  
  #caclulate average expression 
  rpkm_set_mean <- as.matrix(rowMeans(dataSE))
  colnames(rpkm_set_mean) <- set
  
  return(rpkm_set_mean)
}

#calculate mean exp for each stage
stage_averages<- matrix(nrow=dim(RPKM)[1], ncol=5)                          ##### modify here
colnames(stage_averages)<-c('Normal', 'stage1','stage2','stage3','stage4')
rownames(stage_averages)<- rownames(RPKM)

for (i in colnames(stage_averages)){
  this_average<-getAverageExpSampleSet(RPKM, i)
  stage_averages[,i]<- this_average}
head(stage_averages)






####### MFUZZ #######





#we create the Mfuzz object (ExpressionSet)
exprSet=ExpressionSet(assayData=stage_averages)
#As a first step,  we exclude genes with more than 25% of the measurements missing
#-> genes with 0 RPKM in 25% of the conditions
exprSet.r <- filter.NA(exprSet, thres=0.25) 

#Fuzzy  c-means  like  many  other  cluster  algorithms,  does  not  allow  for  missing  values.
#Thus,  we  timelace  remaining  missing  values  by  the  average  values  expression  value
#of  the corresponding gene.

#replace missing values with means
exprSet.f <- fill.NA(exprSet.r,mode="mean") 

#As soft clustering is noise robust, pre-filtering can usually be avoided.
#However, if the number of genes with small expression changes is large, such pre-filtering #may be necessary to reduce noise.
#This function can be used to exclude genes with low standard deviation.
#min.std : threshold for minimum standard deviation.
#If the standard deviation of a gene's expression is smaller than min.std the corresponding #gene will be excluded.

#Calculation the standard deviation shows, however,
#that the transition between low and high values for variation in 
#gene expression is smooth and no particular cut-off point is indicated 
tmp <- filter.std(exprSet.f,min.std=0.0, visu=F) #genes with lower std will be excluded
dim(tmp)
#setdiff(rownames(exprSet.f), rownames(tmp))


#Since the clustering is performed in Euclidian space, the expression values of 
#genes were standardised to have a mean value of zero and a standard deviation of one. 
#Importantly, Mfuzz assumes that the given expression data are fully preprocessed 
#including  any  data  normalisation.
#The  function standardise does  not  replace  the  normalisation step (eg RPKN #normalization).
#This step ensures that vectors of genes with similar changes in expression are close
#in Euclidean space:

exprSet.s <- standardise(tmp)  #NB this step makes expression between genes COMPARABLE
dim(exprSet.s)
# c is the number of clusters 
# m is the fuzziness parameter ->

# we would like to choose a value which prevents clustering of random data. 
#Note, that fuzzy clustering can be tuned in such manner, that random data is not 
#clustered. This is a clear advantage to hard clustering (such as k-means)


m1=mestimate(exprSet.s) 
m1

 # clusrer!
cl <- mfuzz(exprSet.s,c=6,m=m1) 
mfuzz.plot2(exprSet.s,cl=cl,mfrow=c(2,3),min.mem=0.6, cex.axis=0.9, time.labels=colnames(stage_averages), new.window = T)



#To extract list of genes belonging to the cluster cores, the acore function can be used.
cores<-acore(exprSet.s,cl=cl,min.acore=0.6)

#make a list object with gene names for each cluster
cluster_data<-list()
# get total number of genes after filetering : std=0.3, membership  >0.6
gene_sum=0

for ( i in 1:dim(summary(cores))[1]){
  print (paste( "Cluster", i, ":" , dim(cores[[i]])[1] ))
  this_cluster<-as.character(paste0("cluster",i))
  cluster_data[[this_cluster]]<-rownames(cores[[i]])
  gene_sum=gene_sum+dim(cores[[i]])[1]
}

#cluster_data
gene_sum
