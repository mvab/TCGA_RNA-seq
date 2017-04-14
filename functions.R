getAutophagyGenes <- function(dataSE){
  
  autophagy_genes<- as.vector(read.table("autopahagy_genes.txt", as.is = T, header = FALSE))$V1
  all_genes<-rownames(dataSE)
  shared <- intersect(autophagy_genes,all_genes)
  newdataSE<- dataSE[c(shared),]
  
  print(paste0("Total number of genes: ", length(all_genes)))
  print(paste0("Autophagy genes: ", length(shared)))
  
  return(newdataSE)
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