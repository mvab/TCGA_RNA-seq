https://www.biostars.org/p/62583/
  https://www.biostars.org/p/112735/
  https://support.bioconductor.org/p/78966/
  https://support.bioconductor.org/p/64585/
  https://www.biostars.org/p/61523/



library(TxDb.Hsapiens.UCSC.hg19.knownGene)

hg19GeneLengths <- function(symbols)
{
  #require(TxDb.Hsapiens.UCSC.hg19.knownGene) 
  #require(org.Hs.eg.db)
  exons.db = exonsBy(TxDb.Hsapiens.UCSC.hg19.knownGene, by='gene')    
  egs    = unlist(  mget(symbols[ symbols %in% keys(org.Hs.egSYMBOL2EG) ],org.Hs.egSYMBOL2EG) )
  sapply(egs,function(eg)
  {
    exons = exons.db[[eg]]
    exons = reduce(exons)
    sum( width(exons) )
  })
}
# 1:3781, 3783:4310, 4312:6698, 6700:8278, 8280:9212, 9214:14294,14296:14946, 14948: 17372
geneLenghts<-hg19GeneLengths( rownames(dataSE)[c(1:3781, 3783:4310,  4312:6698, 6700:8278,  8280:9212, 9214:14294,14296:14946, 14948: 17372) ] )

#not able to find genes in the hg index:
#CSAG2 (ind 3782) 
#DHX40P1 (4311)
#HBD 6699
#LGTN 8279
#MELK 9212
# SPDYE8P 14295
# TEC 14947

gLen_df=data.frame(Length=geneLenghts)

genes_to_remove=c("CSAG2","DHX40P1","HBD", "LGTN", "MELK","SPDYE8P","TEC")
dataSE2<- dataSE[rownames(dataSE) %in% rownames(gLen_df),]

dim(dataSE2)

dim(gLen_df)

y <- DGEList(counts=dataSE2,genes=gLen_df)
y <- calcNormFactors(y)
RPKM <- rpkm(y)










#### adding genes: SKIP THIS ######

dge$genes <- data.frame(SYMBOL=rownames(dge))

#library(mygene)
#q <-queryMany(dge$genes$SYMBOL, scopes="symbol", fields=c("entrezgene"), species="human")

#che k if got all the genes
#'%nin%' <- Negate('%in%')
#dge$genes$SYMBOL[dge$genes$SYMBOL%nin% q$query]

#qtest<-queryMany(c("CSDAP1", "MGC2752","NCRNA00171"), scopes="symbol", fields=c("entrezgene","ensembl.gene"), species="human")
# 3 genes taht dont have entrez
#1 ENSG00000261614  117.6303      CSDAP1 ENSG00000261614
#2 ENSG00000268784  117.6677     MGC2752 ENSG00000268784
#3 ENSG00000229653  117.6849  NCRNA00171 ENSG00000229653
#q[q$query=='CSDAP1',]$`_id`<-NA
#q[q$query=='MGC2752',]$`_id`<-NA
#q[q$query=='NCRNA00171',]$`_id`<-NA

#remove duplicated with ENSembl id
#q<-q[!(startsWith(q$`_id`, "ENSG0000") & !is.na(q$`_id`)),]
#remove just duplicates
#q<- subset(q, !duplicated(query))
#update
#dge$genes$ENTREZID <- q$entrezgene





#add genes annotation ----- FIGURE THIS OUT LATER
#library("BSgenome.Hsapiens.UCSC.hg19")
#geneid <- rownames(y) 
#genes <- select(BSgenome.Hsapiens.UCSC.hg19, keys=geneid, columns=c("SYMBOL", "TXCHROM"), keytype="ENTREZID")
#dim(genes)
#genes <- genes[!duplicated(genes$ENTREZID),]


#editing genes names in dataSE: these 3 genes have entrzid instaead of geneame
#2    155060 LOC155060
#3      8225    GTPBP6
#4     90288   EFCAB12