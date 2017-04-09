getAutophagyGenes <- function(dataSE){
  
  autophagy_genes<- as.vector(read.table("autopahagy_genes.txt", as.is = T, header = FALSE))$V1
  all_genes<-rownames(dataSE)
  shared <- intersect(autophagy_genes,all_genes)
  newdataSE<- dataSE[c(shared),]
  
  print(paste0("Total number of genes: ", length(all_genes)))
  print(paste0("Autophagy genes: ", length(shared)))
  
  return(newdataSE)
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