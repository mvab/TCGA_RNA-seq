
setwd("~/Bioinformatics MSc UCPH/0_MasterThesis/TCGAbiolinks/CBL_scripts/autophagy_genes_lists")

#autophagy genes cathegories
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
  lipid=lipid,
  phosphatidyl = phosphatidyl,
  endo_exosomes = endo_exosomes,
  transport = transport,
  rabs = rabs,
  docking_fusion = docking_and_fusion,
  mito = mito,
  autoph_core = autophagic_core,
  transcr_factors = transcription_factors,
  receptors_ligands = receptors_and_ligands,
  mTOR_induction = mTOR_induction,
  lysosome = lysosome
)

all_auto_genes <- vector(mode="character", length=0)
for (i in 1:length(autophagy_genes)){
  print(paste0(names(autophagy_genes[i]), ":       ",length(unique(autophagy_genes[[i]])),"/",length(autophagy_genes[[i]])," unique in this group" ))
  all_auto_genes<-append(all_auto_genes,unique(autophagy_genes[[i]]))
}
length(all_auto_genes) #1392
length(unique(all_auto_genes))#1183
all_unique_auto_genes<-unique(all_auto_genes)

write(unique(all_unique_auto_genes), file = "autopahagy_genes.txt")




gene_functions <- vector(mode="character", length=0)

for (name in names(autophagy_genes)) {
    print(name)
    for (i in all_unique_auto_genes){
      
      if(i %in% autophagy_genes[[name]]){
          if (is.na(gene_functions[i])) {
            gene_functions[i]<- name
          }else{
            gene_functions[i]<- "multifunction"
          }
      }
      
    }
}  

save(gene_functions, file="autophagy_genes_functions.rda")


auto_matrix<-as.matrix(gene_functions)
auto_matrix$genes <-rownames(auto_matrix)
rownames(auto_matrix)<-NULL
dim(auto_matrix)
head(auto_matrix)

save(auto_matrix, file="autophagy_functions.rda")


#####



dataSE<-get(load("~/Bioinformatics MSc UCPH/0_MasterThesis/TCGAbiolinks/CBL_scripts/data/BRCA_Illumina_HiSeqnew_updatedSE_allFemale_PreprocessedData_wo_batch_GC_010_updatedSE_allTypes_allStages_Female.rda"))

# current all genes
all_genes<-rownames(dataSE)
# all possible autophagy genes
length(all_unique_auto_genes)
# autophagy that are not in current genes
removed_autophagy <- all_unique_auto_genes[which(!all_unique_auto_genes %in% all_genes)]
length(removed_autophagy)


## print what genes were removed in each functonal group with filtering

for (name in names(autophagy_genes)) {
  print(name)
  removed <- vector(mode="character", length=0)

  for (i in 1:length(removed_autophagy)){

    if(removed_autophagy[i] %in% c(autophagy_genes[[name]])){
      removed<-append(removed, removed_autophagy[i])
    }
  }
  print(paste0("Removed ",length(removed), " genes in group ", name, ", they are:"))
  print(removed)
  cat("\n")
  removed<-NULL
}  


