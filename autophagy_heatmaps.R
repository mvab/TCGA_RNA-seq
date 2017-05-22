setwd("~/Bioinformatics MSc UCPH/0_MasterThesis/TCGAbiolinks/CBL_scripts/data/DEA/")

samples.matrix<-get(load("samplesMatrix_full.rda"))
v<- get(load("v_expression_for_hm.rda"))
autophagy_genesDE<- get(load("autophagy_genesDE.rda"))
dim(v)


## getting autophagy lists
setwd("~/Bioinformatics MSc UCPH/0_MasterThesis/TCGAbiolinks/CBL_scripts/autophagy_genes_lists")

#### autophagy genes cathegories #####
lipid<- as.vector(read.table("lipid.txt", as.is = T, header = FALSE))$V1
phosphatidyl<- as.vector(read.table("phosphatidylinositol.txt", as.is = T, header = FALSE))$V1
endo_exosomes<- as.vector(read.table("endo_exosomes.txt", as.is = T, header = FALSE))$V1
transport<- as.vector(read.table("transport.txt", as.is = T, header = FALSE))$V1
rabs<- as.vector(read.table("rabs_and_effectors.txt", as.is = T, header = FALSE))$V1
docking_and_fusion<- as.vector(read.table("docking_and_fusion.txt", as.is = T, header = FALSE))$V1
mito<- as.vector(read.table("mito.txt", as.is = T, header = FALSE))$V1
autophagic_core<- as.vector(read.table("autophagic_core.txt", as.is = T, header = FALSE))$V1
transcription_factors<- as.vector(read.table("transcription_factors.txt", as.is = T, header = FALSE))$V1
receptors_and_ligands<- as.vector(read.table("receptors_and_ligands.txt", as.is = T, header = FALSE))$V1
mTOR_induction<- as.vector(read.table("induction_mTOR.txt", as.is = T, header = FALSE))$V1
lysosome<- as.vector(read.table("lysosome.txt", as.is = T, header = FALSE))$V1

autophagy_genes <- list (
  Lipid=lipid,
  Phosphatidyl = phosphatidyl,
  Endo_exosomes = endo_exosomes,
  Transport = transport,
  RABs = rabs,
  Docking_and_fusion = docking_and_fusion,
  Mitophagy = mito,
  Autophagic_core = autophagic_core,
  Transcription_factors = transcription_factors,
  Receptors_and_ligands = receptors_and_ligands,
  mTOR_induction = mTOR_induction,
  Lysosome = lysosome
)

######

setwd("~/Bioinformatics MSc UCPH/0_MasterThesis/TCGAbiolinks/CBL_scripts/data/DEA/")








#### genes classes setup #####
###### patient classes setup #####

#get only interesting cols
names(samples.matrix)
hm.design<-subset(samples.matrix[,c('tumourTypes','tumourStages','PAM50')])#, 'tss'
row.names(hm.design)<-samples.matrix$barcode

#recorder annotations
hm.design$PAM50 = factor(hm.design$PAM50, levels = c("Basal-like", "HER2-enriched", "Luminal A", "Luminal B", "Normal-like", "Normal"))
hm.design$tumourStages = factor(hm.design$tumourStages, levels = c("stage1", "stage2", "stage3", "stage4", "unknown")) 

hm.design$tumourTypes = factor(hm.design$tumourTypes,  
                               levels = c("Normal",  "Mucinous adenocarcinoma", 
                                          "Ductal carcinoma","Lobular carcinoma", "Ductual mixed with others ",
                                          "Ductal and lobular mixed", "Metaplastic carcinoma" ,"Other")) #other is a mix of few samples of different ones

# colours
PAM50Col <- c( "#E69F00", "#0072B2", "#F0E442","red3", "#CC79A7","#009E73")##D55E00->altred
tumourTypCol = c("#009E73",  "#0072B2", "#D55E00","#56B4E9","#E69F00","#F0E442","#CC79A7", "white")
tumourStgCol = c( "red","green","blue", "yellow2", "white")

names(PAM50Col)<-levels(hm.design$PAM50)
names(tumourTypCol)<-levels(hm.design$tumourTypes)
names(tumourStgCol)<-levels(hm.design$tumourStages)

annColour <-list(
  PAM50=PAM50Col,
  tumourTypes=tumourTypCol,
  tumourStages=tumourStgCol
  )



######


#extract only genes in auto group of interest
for (i in 1:length(names(autophagy_genes))){
  name<-(names(autophagy_genes)[i])
  group<-autophagy_genes[[name]]
  
  geneExp <- v$E[rownames(v$E) %in% group , ]
  print (paste0( name,"; Genes: " dim(geneExp)[1], " --- making a heatmap.."))
  
  pheatmap::pheatmap(mat = as.matrix(geneExp), color = brewer.pal(name = "YlGnBu", n = 9),
                   clustering_distance_rows = 'manhattan', 
                   clustering_distance_cols = 'manhattan', 
                   #scale="row",
                   annotation_col=hm.design,
                   annotation_colors = annColour,
                   #annotation_row=ok_g_functions,
                   cluster_cols = T, cluster_rows = T, 
                   show_rownames = T,show_colnames = F,
                   fontsize = 5,
                   main=paste0("Autophagy-related genes in group: ", name))
}



