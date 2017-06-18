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
  return(shared)
  
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


addClinData<-function(sample.matrix){
  
  #clindata<-get(load("clinData_prepared.rda")) #old from elena' data
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



renameMorph<- function(samples.matrix){
  
  morph<-as.vector(samples.matrix$tumourTypes)
  morph_renm <- vector(mode="character", length=0)
  
  for (i in 1:length(morph)){
    
    if (is.na(morph[i])){
      morph_renm[i]<-  "Normal" 
    } else if (morph[i] == "unknown") {
      morph_renm[i]<- "Other"  
    } else if (morph[i] == "female_84803") {
      morph_renm[i]<- "Mucinous adenocarcinoma"
    } else if (morph[i] =="female_85003" ){
      morph_renm[i]<- "Ductal carcinoma"
    } else if  (morph[i] =="female_85203"){
      morph_renm[i]<- "Lobular carcinoma"
    } else if (morph[i] =="female_85233"){
      morph_renm[i] <- "Ductual mixed with others"
    } else if (morph[i] =="female_85223" ){
      morph_renm[i]<- "Ductal and lobular mixed"  
    } else if (morph[i] =="female_85753") {
      morph_renm[i]<- "Metaplastic carcinoma"    
    }  
  }
  #length(morph)
  #length(morph_renm)
  samples.matrix$tumourTypes<-morph_renm
  return(samples.matrix)
}

addXtraPAMandNormal<-function(samples.matrix){
  
  #extradata<-get(load("extradata_var_w_old.rda")) #991
  extradata<-get(load("newPAM50list.rda")) #872

  #add extra patient information #578
  subdiagnosis <- merge(samples.matrix, extradata, by="patient", all.x=T) 
  subdiagnosis <- subdiagnosis[order(subdiagnosis$myorder), ]
  
  # get barcodes of samples with PAM50 recorded
  diag <- subdiagnosis[!is.na(subdiagnosis$PAM50new), ] 
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
  
  # remove samples with no subdtype label and paired tumor sample from dataframe, subtype vector and ID-information.
  samples.matrix <- samples.matrix[samples.matrix$barcode %in% to_keep, ]    
  subdiagnosis <- subdiagnosis[subdiagnosis$barcode %in% to_keep, ]
  
  # make normal samples from pairs "normal" while retaining info on subtype.
  subdiagnosis$PAM50 <- ifelse(subdiagnosis$condition == "normal", "Normal", as.character(subdiagnosis$PAM50new))
  
  print ("Keeping PAM50 samples and all normals:")
  print( table(subdiagnosis$PAM50))
  subdiagnosis$PAM50new<-NULL
  
  return(list(samples.matrix=subdiagnosis, samplesToKeep=to_keep, normal_barcodes = paired))
}  


addPAM50andNormal<-function(samples.matrix){
  
  # get information on subtype/pation details
  #dataSubt <- TCGAquery_subtype(tumor = "BRCA") 
  
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

addOnlyTopTSS <- function(samples.matrix){
  
  #filtering out patient with neither Positive nor Negative
  #tss_above50samples <-c("BH","A2","E2","A8","D8","E9", "AR", "B6","AC" ) #660 patients
  tss <-c("BH") 
  
  #samples.matrix <- samples.matrix[samples.matrix$tss %in% tss_above50samples,] 
  samples.matrix <- samples.matrix[samples.matrix$tss %in% tss,] 
  
  return (samples.matrix)
}  


getGeneLenght<-function(dataSE){
  
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



## funtion to just remove unknoen anf other samples 
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





addBRCAReceptorStatus <- function(samples.matrix){
  
  # get information on subtype/pation details
  dataSubt <- TCGAquery_subtype(tumor = "BRCA") 
  
  #add extra patient information, but some have NAs; if FALSE- no NAs but lose ~200 patients #~735
  subdiagnosis <- merge(samples.matrix, dataSubt, by="patient", all.x=TRUE) 
  subdiagnosis <- subdiagnosis[order(subdiagnosis$myorder), ]
  
  # we only have female samples, so NA in gender can be replaced with FEMALE
  subdiagnosis[c("Gender")][is.na(subdiagnosis[c("Gender")])] <- "FEMALE"
  
  # note:  the are many thing not relevant - remove them and then do complete cases
  
  # testing how many pations how no NAs for ER PR HER cols
  complete_cases_receptors<- subdiagnosis[complete.cases(subdiagnosis[, c(1:12, 15:17)]),]# 704 patients 
  
  #filtering out patient with neither Positive nor Negative
  other_recp_status <-c("Indeterminate", "Not Performed", "Performed but Not Available", "Equivocal", "Not Available")
  complete_cases_receptors$ER.Status[complete_cases_receptors$ER.Status %in% other_recp_status] <- NA
  complete_cases_receptors$PR.Status[complete_cases_receptors$PR.Status %in% other_recp_status] <- NA
  complete_cases_receptors$HER2.Final.Status[complete_cases_receptors$HER2.Final.Status %in% other_recp_status] <- NA
  
  complete_cases_receptors<- complete_cases_receptors[complete.cases(complete_cases_receptors[, c(1:12, 15:17)]),]
  complete_cases_receptors<- subset(complete_cases_receptors[, c(1:12, 15:17)]) # 646 patients
  dim(complete_cases_receptors)
  
  
  # replacing Positive and Negative in receptor column with 0 and 1 
  
  complete_cases_receptors$ER.Status <- gsub('Positive', 1, complete_cases_receptors$ER.Status)
  complete_cases_receptors$ER.Status <- gsub('Negative', 0, complete_cases_receptors$ER.Status)
  complete_cases_receptors$PR.Status <- gsub('Positive', 1, complete_cases_receptors$PR.Status)
  complete_cases_receptors$PR.Status <- gsub('Negative', 0, complete_cases_receptors$PR.Status)
  complete_cases_receptors$HER2.Final.Status <- gsub('Positive', 1, complete_cases_receptors$HER2.Final.Status)
  complete_cases_receptors$HER2.Final.Status <- gsub('Negative', 0, complete_cases_receptors$HER2.Final.Status)
  
  receptor_statusNum = as.matrix(as.data.frame(lapply(subset(complete_cases_receptors[13:15]), as.numeric)))
  receptor_status$sum <- rowSums(receptor_statusNum)
  
  # adding subtype classification:
  
  ### Luminal     ###   ER+PR+HER+ or ER+PR+HER- , so 3 or 2
  ### HER2 over   ###   ER-PR-HER+  , so 1
  ### Triple Neg  ###   ER-PR-HER-  , so 0
  
  receptSubtype <- vector(mode="character", length=0)
  for (i in 1:length(receptor_status$sum)){
    if (receptor_status$sum[i] > 1) {
      receptSubtype[i]<-"Luminal"
    } else if (receptor_status$sum[i] == 1) {
      receptSubtype[i]<-"HER2"
    } else {
      receptSubtype[i]<-"TNR"
    }
  }
  receptor_status$receptSubtype <- receptSubtype
  # add receptSubtype to the main sample matrix and remove sepate receptors
  
  complete_cases_receptors<-subset(complete_cases_receptors[, c(1:12)])
  complete_cases_receptors$receptSubtype <- receptSubtype
  
  return (complete_cases_receptors)
}

samples.matrixRecep<-addBRCAReceptorStatus(samples.matrix)
# if used addBRCAReceptorStatus DO THIS:
dataSE<- dataSE[,c(samples.matrixRecep$barcode)]
dim(dataSE)


keepOnly<-function(dataSE, samples.matrix){
  
  #morphology
  samples.matrix[samples.matrix$tumourTypes=="Lobular carcinoma" | samples.matrix$tumourTypes== "Ductal carcinoma"| samples.matrix$tumourTypes== "Normal",]$barcode -> morph
  length(morph)
  samples.matrix<-samples.matrix[samples.matrix$barcode %in% morph,]
  dim(samples.matrix)
  
  #stages
  samples.matrix[samples.matrix$tumourStages!="unknown" & samples.matrix$tumourStages!="stage4",]$barcode -> stage
  length(stage)
  samples.matrix<-samples.matrix[samples.matrix$barcode %in% stage,]
  dim(samples.matrix)
  
  #pam50
  samples.matrix[samples.matrix$PAM50!="Normal-like",]$barcode -> pam50
  length(pam50)
  samples.matrix<-samples.matrix[samples.matrix$barcode %in% pam50,]
  dim(samples.matrix)
  
  to_keep <- samples.matrix$barcode
  length(to_keep)
  
  dataSE<-dataSE[, colnames(dataSE) %in% to_keep ]
  dim(dataSE)
  
  return(list(dataSE=dataSE, samples.matrix=samples.matrix))
  
}

onlysamples<-keepOnly(dataSE, samples.matrix)
dataSE<-onlysamples$dataSE
samples.matrix<-onlysamples$samples.matrix



## renaming differing PAM labels (514, 39 renamed)
addAltPAM50<-function(samples.matrix){
  
  PAM50_upd<-get(load("potentially_wrong_PAM.rda"))
  names(PAM50_upd)[2]<-"PAM50upd"
  
  samples.matrix <- merge(samples.matrix, PAM50_upd, by="patient", all.x=TRUE) 
  samples.matrix <- samples.matrix[order(samples.matrix$myorder), ]
  samples.matrix$PAM50upd <- ifelse(samples.matrix$condition == "normal", "normal", as.character(samples.matrix$PAM50upd))
  
  for (i in 1:length(samples.matrix$patient)){
    if (is.na(samples.matrix$PAM50upd[i])){
      samples.matrix$PAM50upd[i]<-samples.matrix$PAM50[i]
    }
  }
  return(samples.matrix)
}
samples.matrix<-addAltPAM50(samples.matrix)
names(samples.matrix)
