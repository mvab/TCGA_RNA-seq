library('Mfuzz')
library('GenomicFeatures')
library('DESeq')
library('edgeR')

#Nomalization of RNA-seq data:
  
#create the object with all count files, count_files_folder is the path to the folder #containing all the count tables
files=list.files(count_files_folder)
raw=readDGE(files, path=count_files_folder, group=c(1:length(files)), columns=c(1,2))

#create transcripts database from your gtf or gff.

#annotation is the absolute path to your annotation file

#change gtf to gff if you have a gff file.

txdb=makeTxDbFromGFF(annotation,format='gtf')
#then collect the exons per gene id

#if in your annotation file, the gene attribute is not 'gene', change it here
exons.list.per.gene=exonsBy(txdb,by='gene')
#then for each gene, reduce all the exons to a set of non overlapping exons, calculate their #lengths (widths) and sum then
exonic.gene.sizes=as.data.frame(sum(width(reduce(exons.list.per.gene))))
colnames(exonic.gene.sizes)='gene_length_bp'
#our raw data table containing all the samples
datafile=raw$counts

#remove all the genes with a 0 count in all samples
data = datafile[apply(datafile,1,sum)!=0,]
#determine number of studied samples
nblib= dim(data)[2]

#create a factice vector for DESeq normalization
conds = factor(1:nblib)
#normalize raw read counts by library size with DESeq method
cds = newCountDataSet(data, conds)
cds = estimateSizeFactors(cds)
datanorm = t(t(data)/sizeFactors(cds))
colnames(datanorm)=paste(colnames(data),'normalized_by_DESeq', sep='_')
#merge the raw read counts and the read counts normalized by DESeq
alldata = merge(data, datanorm, by='row.names', all=T)
alldata_tmp=alldata[,-1]

rownames(alldata_tmp)=alldata[,1]
alldata=alldata_tmp
#merge the raw and normalized read counts with the genes length (in bp)
alldata_tmp = merge(alldata, exonic.gene.sizes, by='row.names', all.x=T)
alldata=alldata_tmp[,-1]
rownames(alldata)=alldata_tmp[,1]
#recover of the start of normalized read count columns and the end
start=length(files)+1
end=length(files)*2
#calcul of the RPKM
data_norm=merge(datanorm,exonic.gene.sizes, by='row.names', all.x=T)
data_norm_tmp=data_norm[,-1]
rownames(data_norm_tmp)=data_norm[,1]
data_norm=data_norm_tmp
rpkm=rpkm(data_norm[,1:dim(data_norm)[2]-1],gene.length=data_norm$'gene_length_bp', normalized.lib.sizes=FALSE, log=FALSE)
colnames(rpkm)=paste(colnames(data),'normalized_by_DESeq_and_divided_by_gene_length', sep='_')

#merge of rpkm with the table containing raw and normalized read counts and gene #length
alldata_tmp = merge(alldata, rpkm, by='row.names', all.x=T)
alldata=alldata_tmp[,-1]
rownames(alldata)=alldata_tmp[,1]
#write of this table containing the raw read counts and the different normalizations

output=getwd()
write.table(alldata, paste(output,'normalized_counts_all_genes.txt',sep='/'), sep='\t', quote=F, row.names=T, dec='.')



#MFuzz analysis :
  
#determine the first RPKM column in the alldata object ->
  
#RPKM are used by Mfuzz for genes clustering
first_rpkm_column=dim(alldata)[2]-length(files)+1

#here we create a matrix containing all the RPKM columns, used by Mfuzz for the #clustering
exprs=as.matrix(alldata[,first_rpkm_column:dim(alldata)[2]])
#and for each time value containing replicates, we calculate the RPKM count means
#if there are no replicates, we keep the initial RPKM

#here you have to create a 'time' vector, like that :

#time='time1,time1,time1,time2,time2,time2,time3???

#-> give the time value of each file by respecting the same order in the vector than the files #in the folder.
count=1
for ( i in unique(time) ){
  
  if ( dim(as.data.frame(exprs[,which(time==i)]))[2] == 1 ){
    mean_rpkm=data.frame(exprs[,which(time==i)])
  } else {
    mean_rpkm=data.frame(rowMeans(exprs[,which(time==i)]))
  }
  colnames(mean_rpkm)=i
    
  if (count == 1){
    mean_rpkm_ok=mean_rpkm
  } else {
    mean_rpkm_ok=merge(mean_rpkm_ok,mean_rpkm,by='row.names')
    rownames(mean_rpkm_ok)=mean_rpkm_ok[,1]
    mean_rpkm_ok=mean_rpkm_ok[,-1]
  }
  count=count+1
}
#here we have a RPKM matrix containing one column per time value (and not one column #per sample)
exprs_with_time=as.matrix(mean_rpkm_ok, header=TRUE, sep='\t',row.names=1,as.is=TRUE)
#we create the Mfuzz object (ExpressionSet)
exprSet=ExpressionSet(assayData=exprs_with_time)
#----------------------------------------
#As a first step,  we exclude genes with more than 25% of the measurements missing
#-> genes with 0 RPKM in 25% of the conditions
exprSet.r=filter.NA(exprSet, thres=0.25)
#----------------------------------------
#----------------------------------------
#Fuzzy  c-means  like  many  other  cluster  algorithms,  does  not  allow  for  missing  values.
#Thus,  we  timelace  remaining  missing  values  by  the  average  values  expression  value
#of  the corresponding gene.
#Methods for replacement of missing values. Missing values should be indicated by NA in #the expression matrix.
#Mode method for replacement of missing values:
#mean- missing values will be replaced by the mean expression value of the gene,
#median- missing values will be replaced by the median expression value of the gene,
#knn- missing values will be replaced by the averging over the corresponding expression    #values of the k-nearest neighbours,
#knnw-same replacement method as knn, but the expression values averaged are #weighted by the distance to the corresponding neighbour
exprSet.f=fill.NA(exprSet.r,mode='mean')
#----------------------------------------
#----------------------------------------
#As soft clustering is noise robust, pre-filtering can usually be avoided.
#However, if the number of genes with small expression changes is large, such pre-filtering #may be necessary to reduce noise.
#This function can be used to exclude genes with low standard deviation.
#min.std : threshold for minimum standard deviation.
#If the standard deviation of a gene's expression is smaller than min.std the corresponding #gene will be excluded.
tmp=filter.std(exprSet.f,min.std=0, visu=FALSE)
#----------------------------------------
#----------------------------------------
#Since  the  clustering  is  performed  in  Euclidian  space,  the  expression  values  of  genes  #were standardised to have a mean value of zero and a standard deviation of one.  This #step ensures
#that vectors of genes with similar changes in expression are close in Euclidean space
#Importantly, Mfuzz assumes that the given expression data are fully preprocessed #including  any  data  normalisation.
#The  function standardise does  not  replace  the  normalisation step (eg RPKN #normalization).
exprSet.s=standardise(tmp)
#----------------------------------------
#----------------------------------------
#clustering

#here you can choose the number of clusters, as example here, 4 clusters

m1=mestimate(exprSet.s)
cl=mfuzz(exprSet.s,c=4,m=m1)

#----------------------------------------

#here you have to choose the membership cut-off to use with Mfuzz -> see the Mfuzz #paper : http://www.bioinformation.net/002/000200022007.pdf
#you can give one value, eg '0.7' or several values separated by ',', eg '0.5,0.7'

for (membership in c(0.5,0.7)) {
  membership=as.numeric(membership)
  #create one output folder per membership
  
  dir=paste(output,paste('cluster_with_membership',membership, sep='_'),sep='/')
  dir.create(dir, showWarnings = FALSE)
  #----------------------------------------
  #membership cut-off part and plot clusters
  pdf(paste(dir,paste(paste('clusters_Mfuzz_membership_equals_',membership,sep=''),'.pdf',sep=''), sep='/'))
  mfuzz.plot2(exprSet.s,cl=cl,time.labels=unique(time),min.mem=membership, colo='fancy', x11=FALSE)
  
  dev.off()
  #generates one genes list per cluster
  acore.list=acore(exprSet.s,cl=cl,min.acore=membership)
  print(paste('Membership',membership,sep=' : '))
  for (cluster in 1:nb_clusters){
    print(paste(paste('Number of genes in cluster', cluster, sep=' '),dim(acore.list[[cluster]])[1], sep=' : '))
    cluster_table=merge(alldata,acore.list[[cluster]][2], by='row.names', all.y=TRUE)
    write.table(cluster_table,paste(dir,paste(paste('list_of_genes_in_cluster',cluster,sep='_'),'.txt'),sep='/'), sep='\t',row.names=F, dec='.')
  }