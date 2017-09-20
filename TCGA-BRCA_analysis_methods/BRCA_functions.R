getAutophagyGenes <- function(dataSE){
  
  autophagy_genes<- as.vector(read.table("autopahagy_genes.txt", as.is = T, header = FALSE))$V1
  all_genes<-rownames(dataSE)
  shared <- intersect(autophagy_genes,all_genes)
  newdataSE<- dataSE[c(shared),]
  
  print(paste0("Total number of genes: ", length(all_genes)))
  print(paste0("Autophagy genes: ", length(shared)))
  
  return(newdataSE)
} 

## getting autophagy genes DE in a supplies gene list
sharedWithAuto <-function (gene_list){
  autophagy_genes<- as.vector(read.table("~/Bioinformatics MSc UCPH/0_MasterThesis/TCGAbiolinks/CBL_scripts/data/autopahagy_genes.txt", as.is = T, header = FALSE))$V1
  shared <- intersect(autophagy_genes,gene_list)
  print(paste0("Total number of genes in the list : ", length(gene_list)))
  print(paste0("Autophagy genes among these: ", length(shared)))
  return(as.vector(unlist(shared)))
  
}
sharedWithAutoCORE <-function (gene_list){
  autophagy_genes<- as.vector(read.table("~/Bioinformatics MSc UCPH/0_MasterThesis/TCGAbiolinks/CBL_scripts/data/autophagic_core.txt", as.is = T, header = FALSE))$V1
  shared <- intersect(autophagy_genes,gene_list)
  print(paste0("Total number of genes in the list : ", length(gene_list)))
  print(paste0("Autophagy CORE genes among these: ", length(shared)))
  return(shared)
  
}
sharedWithAutoTF <-function (gene_list){
  autophagy_genes<- as.vector(read.table("~/Bioinformatics MSc UCPH/0_MasterThesis/TCGAbiolinks/CBL_scripts/data/transcription_factors.txt", as.is = T, header = FALSE))$V1
  shared <- intersect(autophagy_genes,gene_list)
  print(paste0("Total number of genes in the list : ", length(gene_list)))
  print(paste0("Autophagy TF genes among these: ", length(shared)))
  return(shared)
  
}

renameMorph<- function(samples.matrix){
  # this functions renames mophologilical codes to actual names
  
  morph<-as.vector(samples.matrix$tumourTypes)
  morph_renm <- vector(mode="character", length=0)
  
  for (i in 1:length(morph)){
    
    if (is.na(morph[i])){
      morph_renm[i]<-  "Normal" 
    } else if (morph[i] == "unknown") {
      morph_renm[i]<- "Other"  
    } else if (morph[i] == "female_84803") {
      morph_renm[i]<- "Mucinous carcinoma"
    } else if (morph[i] =="female_85003" ){
      morph_renm[i]<- "Ductal carcinoma"
    } else if  (morph[i] =="female_85203"){
      morph_renm[i]<- "Lobular carcinoma"
    } else if (morph[i] =="female_85233"){
      morph_renm[i] <- "Ductual mixed with others"
    } else if (morph[i] =="female_85223" ){
      morph_renm[i]<- "Ductal/lobular mixed"  
    } else if (morph[i] =="female_85753") {
      morph_renm[i]<- "Metaplastic carcinoma"    
    }  
  }

  #length(morph)
  #length(morph_renm)
  samples.matrix$tumourTypes<-morph_renm
  return(samples.matrix)
}



addClinData<-function(sample.matrix){
  # this functiona add metadata to all samples 
  
  clindata <- get(load("cleaned_metadata.rda")) # new from 2 sources provided by Kristoffer
  
  samples.matrix_clinincal <- merge(samples.matrix, clindata, by="patient", all.x=TRUE) 
  samples.matrix_clinincal <- samples.matrix_clinincal[order(samples.matrix_clinincal$myorder), ]
  
  # replacing NAs with Normal in relevant categories
  samples.matrix_clinincal$substage <- ifelse(samples.matrix_clinincal$condition == "normal", "Normal", as.character(samples.matrix_clinincal$substage))
  samples.matrix_clinincal$histology <- ifelse(samples.matrix_clinincal$condition == "normal", "Normal", as.character(samples.matrix_clinincal$histology))
  samples.matrix_clinincal$morphology <- ifelse(samples.matrix_clinincal$condition == "normal", "Normal", as.character(samples.matrix_clinincal$morphology))
  samples.matrix_clinincal$tumour_size <- ifelse(samples.matrix_clinincal$condition == "normal", "Normal", as.character(samples.matrix_clinincal$tumour_size))
  samples.matrix_clinincal$nodes <- ifelse(samples.matrix_clinincal$condition == "normal", "Normal", as.character(samples.matrix_clinincal$nodes))
  samples.matrix_clinincal$metastasis <- ifelse(samples.matrix_clinincal$condition == "normal", "Normal", as.character(samples.matrix_clinincal$metastasis))
  
  samples.matrix_clinincal$tumourStages <- ifelse(samples.matrix_clinincal$condition == "normal", "Normal", as.character(samples.matrix_clinincal$tumourStages))
  
  return(samples.matrix_clinincal)
}


addPAM50annotation<-function(samples.matrix){
  #add PAM50 annotation for all samples used in Marina's project  

  #this file contains cutated PAM50 annotation from 2 studies: TCGA Network 2012 and Ciriello 2015. Only samples in  agreement are included
  extradata<-get(load("PAM50_annotation.rda")) #872 samples

  #add extra patient information 
  subdiagnosis <- merge(samples.matrix, extradata, by="patient", all.x=T) 
  subdiagnosis <- subdiagnosis[order(subdiagnosis$myorder), ]
  
  # get barcodes of samples with PAM50 recorded
  diag <- subdiagnosis[!is.na(subdiagnosis$PAM50new), ] 
  diag <- diag$barcode 
  
  # get barcodes of patients with paired samples

  #get patients with 2 entries == 224
  paired <- samples.matrix[ duplicated(samples.matrix$participant, fromLast=FALSE) | duplicated(samples.matrix$participant, fromLast=TRUE),] 
  paired <- as.character(paired$participant)
  #get normal samples from paired data == 112
  paired <- samples.matrix[samples.matrix$participant %in% paired & samples.matrix$condition == "normal", ]  
  paired <- as.character(paired$barcode)
  
  # concatnate barcodes to keep
  #join NA and paired cancer samples
  to_keep <- unique(sort(c(as.character(diag), as.character(paired)))) 
  
  # remove samples with no subdype label and paired tumor sample from dataframe, subtype vector and ID-information.
  samples.matrix <- samples.matrix[samples.matrix$barcode %in% to_keep, ]    
  subdiagnosis <- subdiagnosis[subdiagnosis$barcode %in% to_keep, ]
  
  # make normal samples from pairs "normal" while retaining info on subtype.
  subdiagnosis$PAM50 <- ifelse(subdiagnosis$condition == "normal", "Normal", as.character(subdiagnosis$PAM50new))
  
  print ("Keeping PAM50 samples and all normals:")
  print( table(subdiagnosis$PAM50))
  subdiagnosis$PAM50new<-NULL
  
  return(list(samples.matrix=subdiagnosis, samplesToKeep=to_keep, normal_barcodes = paired))
}  

addPAM50<-function(samples.matrix){
  # add PAM50 annotation that is available from TCGAbiolinks only (2012 TCGA study)
  
  
  # get information on subtype/pation details
  #library(TCGAbiolinks)
  #dataSubt <- TCGAquery_subtype(tumor = "BRCA") 
    ####OR##### (the same, but you don't need to call TCGAbilonks)
  dataSubt<-get(load("dataSubt.rda"))
  
  #add extra patient information #578
  subdiagnosis <- merge(samples.matrix, dataSubt, by="patient", all.x=T) 
  subdiagnosis <- subdiagnosis[order(subdiagnosis$myorder), ]
  
  # get barcodes of samples with PAM50 recorded
  diag <- subdiagnosis[!is.na(subdiagnosis$PAM50.mRNA), ] 
  diag <- diag$barcode #578
  
  # get barcodes of patients with paired samples
  
  #get patients with 2 entries == 224
  paired <- samples.matrix[ duplicated(samples.matrix$participant, fromLast=FALSE) | duplicated(samples.matrix$participant, fromLast=TRUE),] 
  paired <- as.character(paired$participant)
  #get normal samples from paired data == 112
  paired <- samples.matrix[samples.matrix$participant %in% paired & samples.matrix$condition == "normal", ]  
  paired <- as.character(paired$barcode)
  
  # concatnate barcodes to keep
  
  #join NA and paired cancer samples
  to_keep <- unique(sort(c(as.character(diag), as.character(paired)))) #578+112=626 unique (some normals have subtype label)
  
  #only keep diag
  #to_keep <- (c(as.character(diag)))
  
  # remove samples with no subdtype label and paired tumor sample from dataframe, subtype vector and ID-information.
  samples.matrix <- samples.matrix[samples.matrix$barcode %in% to_keep, ]    
  subdiagnosis <- subdiagnosis[subdiagnosis$barcode %in% to_keep, ]
  
  # make normal samples from pairs "normal" while retaining info on subtype.
  subdiagnosis$PAM50 <- ifelse(subdiagnosis$condition == "normal", "normal", as.character(subdiagnosis$PAM50.mRNA))
  
  print ("Keeping PAM50 samples and all normals:")
  print( table(subdiagnosis$PAM50))
  
  return(list(samples.matrix=subdiagnosis, samplesToKeep=to_keep, normal_barcodes = paired))
}  

removeDuctalOther <- function(dataSE, samples.matrix){
  samples.matrix[samples.matrix$tumourTypes=="Ductual mixed with others",]$barcode -> to_remove
  to_remove<- as.character(to_remove)
  
  samples.matrix<-samples.matrix[!samples.matrix$barcode %in% to_remove,]
  dim(samples.matrix)
  dataSE<-dataSE[, !colnames(dataSE) %in% to_remove ]
  dim(dataSE)
  
  return(list(dataSE=dataSE, samples.matrix=samples.matrix))
}

addSampleData<-function(y, samples.matrix) {
  
  # adding samples information
  #y$samples$condition <- as.factor(samples.matrix$condition)
  y$samples$PAM50 <- as.factor(samples.matrix$PAM50)
  y$samples$morphology <- as.factor(samples.matrix$tumourTypes) # from sample lists!
  y$samples$stages <- as.factor(samples.matrix$tumourStages) # from sample lists!
  y$samples$year <- as.factor(samples.matrix$year_diagnosed)
  y$samples$tss <- factor(samples.matrix$tss)
  y$samples$age <- as.factor(samples.matrix$ageGroups)
  
  
  #fixing NAs
  y$samples$age <- as.character(y$samples$age)
  y$samples[is.na(y$samples$age),]$age<-"UnknownAge"
  y$samples[y$samples$age=="70+",]$age<-"70andUp"
  y$samples[y$samples$age=="< 40",]$age<-"40andDown"
  y$samples$age <- as.factor(y$samples$age)
  
  y$samples$year <- as.character(y$samples$year)
  y$samples[is.na(y$samples$year),]$year<-"UnknownYear"
  y$samples$year <- as.factor(y$samples$year)
  
  
  #stage + PAM50
  y$samples$Group1 <- factor(paste( gsub(" ", "", samples.matrix$PAM50) , samples.matrix$tumourStages,sep="."))
  #morphology +stage
  y$samples$Group2 <- factor(paste( gsub(" ", "", samples.matrix$tumourTypes) , samples.matrix$tumourStages,sep="."))
  #pam50+morphology
  y$samples$Group3 <- factor(paste( gsub(" ", "", samples.matrix$PAM50) , samples.matrix$tumourTypes,sep="."))
  #pam50+age
  y$samples$Group4 <- factor(paste( gsub(" ", "", samples.matrix$PAM50) , y$samples$age,sep="."))
  
  
  
  ## making normal the baselayer
  #y$samples$condition = relevel(y$samples$condition, ref="normal")
  y$samples$PAM50 = relevel(y$samples$PAM50, ref="Normal")
  y$samples$morphology = relevel(y$samples$morphology, ref="Normal")
  y$samples$stages = relevel(y$samples$stages, ref="Normal")
  y$samples$age = relevel(y$samples$age, ref="40andDown") ###########--> should i make one for normal??
  y$samples$Group1 = relevel(y$samples$Group1, ref = "Normal.Normal")
  y$samples$Group2 = relevel(y$samples$Group2, ref = "Normal.Normal")
  y$samples$Group3 = relevel(y$samples$Group3, ref = "Normal.Normal")
  
  return (y)
}





getGeneLenghtOLD<-function(dataSE){
  
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
  #length(my_genes)# 17369
  
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
  
  # col width already has each exon length
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

getGeneLenghtNEW<-function(dataSE){
  
  my_genes<-rownames(dataSE)
  #length(my_genes)# 17366
  # my_genes[duplicated(my_genes) ]
  # [1] "SLC35E2"  is duplicated
  
  # so remove duplicate
  my_genes<-my_genes[!duplicated(my_genes) ]
  dataSE<-dataSE[!duplicated(rownames(dataSE)), ]
  
  
  load("TCGAgeneExons.Rdata") #in /data
  gene_exons<-as.data.frame(geneExons)
  dim(gene_exons)
  colnames(gene_exons)
  
  genes<-unique(gene_exons$gene_name)
  #length(genes) #25508
  
  #compare!
  shared<-intersect(genes,my_genes)
  #length(shared) #17365
  
  
  #get exon data only for my genes 
  my_gene_exons<-gene_exons[gene_exons$gene_name %in% shared, ]
  #length(unique(my_gene_exons$gene_name))
  
  # col width already has each exon length
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



## funtion to just remove unknown stage and other morphologysamples 
removeUnknownOther <- function(dataSE, samples.matrix){
  samples.matrix[samples.matrix$tumourTypes=="OtherMorph",]$barcode -> other
  samples.matrix[samples.matrix$tumourStages=="unknown",]$barcode -> unknown
  
  print(paste( length(other), "other morphology samples to remove"))
  print(paste( length(unknown), "unknown stage samples to remove"))
  to_remove<- (c(as.character(other), as.character(unknown)))
  
  samples.matrix<-samples.matrix[!samples.matrix$barcode %in% to_remove,]
  dim(samples.matrix)
  dataSE<-dataSE[, !colnames(dataSE) %in% to_remove ]
  dim(dataSE)
  
  return(list(dataSE=dataSE, samples.matrix=samples.matrix))
}


