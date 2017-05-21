
setwd("~/Bioinformatics MSc UCPH/0_MasterThesis/TCGAbiolinks/CBL_scripts/data/DEA")


LumAvsLumB_up       <- as.character(read.table("LumAvsLumB_up.txt")$V1)
LumAvsBasal_up       <- as.character(read.table("LumAvsBasal_up.txt")$V1)
LumAvsHER2_up       <- as.character(read.table("LumAvsHER2_up.txt")$V1)
LumAvsNormLike_up       <- as.character(read.table("LumAvsNormLike_up.txt")$V1)
LumAvsNormal_up       <- as.character(read.table("LumAvsNormal_up.txt")$V1)
LumBvsBasal_up       <- as.character(read.table("LumBvsBasal_up.txt")$V1)
LumBvsHER2_up       <- as.character(read.table("LumBvsHER2_up.txt")$V1)
LumBvsNormLike_up       <- as.character(read.table("LumBvsNormLike_up.txt")$V1)
LumBvsNormal_up       <- as.character(read.table("LumBvsNormal_up.txt")$V1)
BasalvsLumA_up       <- as.character(read.table("BasalvsLumA_up.txt")$V1)
BasalvsLumB_up       <- as.character(read.table("BasalvsLumB_up.txt")$V1)
BasalvsHER2_up       <- as.character(read.table("BasalvsHER2_up.txt")$V1)
BasalvsNormLike_up      <- as.character(read.table("BasalvsNormLike_up.txt")$V1)
BasalvsNormal_up       <- as.character(read.table("BasalvsNormal_up.txt")$V1)
HER2vsNormLike_up       <- as.character(read.table("HER2vsNormLike_up.txt")$V1)
HER2vsNormal_up       <- as.character(read.table("HER2vsNormal_up.txt")$V1)
NormalLikevsNormal_up       <- as.character(read.table("NormalLikevsNormal_up.txt")$V1)

LumAvsLumB_down       <- as.character(read.table("LumAvsLumB_down.txt")$V1)
LumAvsBasal_down       <- as.character(read.table("LumAvsBasal_down.txt")$V1)
LumAvsHER2_down       <- as.character(read.table("LumAvsHER2_down.txt")$V1)
LumAvsNormLike_down       <- as.character(read.table("LumAvsNormLike_down.txt")$V1)
LumAvsNormal_down       <- as.character(read.table("LumAvsNormal_down.txt")$V1)
LumBvsBasal_down       <- as.character(read.table("LumBvsBasal_down.txt")$V1)
LumBvsHER2_down       <- as.character(read.table("LumBvsHER2_down.txt")$V1)
LumBvsNormLike_down       <- as.character(read.table("LumBvsNormLike_down.txt")$V1)
LumBvsNormal_down       <- as.character(read.table("LumBvsNormal_down.txt")$V1)
BasalvsLumA_down       <- as.character(read.table("BasalvsLumA_down.txt")$V1)
BasalvsLumB_down       <- as.character(read.table("BasalvsLumB_down.txt")$V1)
BasalvsHER2_down       <- as.character(read.table("BasalvsHER2_down.txt")$V1)
BasalvsNormLike_down      <- as.character(read.table("BasalvsNormLike_down.txt")$V1)
BasalvsNormal_down       <- as.character(read.table("BasalvsNormal_down.txt")$V1)
HER2vsNormLike_down       <- as.character(read.table("HER2vsNormLike_down.txt")$V1)
HER2vsNormal_down       <- as.character(read.table("HER2vsNormal_down.txt")$V1)
NormalLikevsNormal_down       <- as.character(read.table("NormalLikevsNormal_down.txt")$V1)


up_vs_normal <- list(LumA = LumAvsNormal_up,
                     LumB = LumBvsNormal_up,
                     Basal = BasalvsNormal_up,
                     HER2 = HER2vsNormal_up,
                     NormLike = NormalLikevsNormal_up)

down_vs_normal <- list(LumA = LumAvsNormal_down,
                     LumB = LumBvsNormal_down,
                     Basal = BasalvsNormal_down,
                     HER2 = HER2vsNormal_down,
                     NormLike = NormalLikevsNormal_down)
library(VennDiagram)

#### PRETTY DIAGRAM ####
venn <- venn.diagram(up_vs_normal[1:4], 
                     category.names = names(up_vs_normal[1:4]),
                     filename = NULL, lwd = 0.5,
                     main = "Venn diagram Upregulated VS normal", main.cex=1.5,
                     col = "transparent", 
                     fill=c("cornflowerblue","green","yellow","darkorchid1"),
                     alpha =0.5,
                     label.col = c("orange", "white", "darkorchid4", "white", "white", 
                                   "white",    "white", "white", "darkblue", "white", "white", "white", "white", 
                                   "darkgreen", "white"), 
                     fontfamily = "serif",
                     fontface = "bold",
                     cat.col = c("darkblue", "darkgreen", "orange", "darkorchid4"), 
                     sub.cex = 1.5,
                     cat.cex=0.8,
                     cat.fontfamily = 'sans',
                     #rotation.degree = 270,
                     #margin = 0.1,
                     cex=1)

venn <- venn.diagram(up_vs_normal[1:4], 
                     category.names = names(up_vs_normal[1:4]),
                     filename = NULL, lwd = 0.5,
                     main = "Venn diagram Upregulated VS normal", main.cex=1.5,
                     col = "transparent", 
                     fill=c("cornflowerblue","green","yellow","darkorchid1"),
                     alpha =0.5,
                     label.col = c("orange", "white", "darkorchid4", "white", "white", 
                                   "white",    "white", "white", "darkblue", "white", "white", "white", "white", 
                                   "darkgreen", "white"), 
                     labels.default = c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15),
                     fontfamily = "serif",
                     fontface = "bold",
                     cat.col = c("darkblue", "darkgreen", "orange", "darkorchid4"), 
                     sub.cex = 1.5,
                     cat.cex=0.8,
                     cat.fontfamily = 'sans',
                     #rotation.degree = 270,
                     #margin = 0.1,
                     cex=1)
grid.draw(venn)
dev.off()
########

all_up <- unique(sort(c(LumAvsLumB_up, LumAvsBasal_up,
                        LumAvsHER2_up, LumAvsNormLike_up ,
                        LumAvsNormal_up, LumBvsBasal_up, 
                        LumBvsHER2_up, LumBvsNormLike_up ,
                        LumBvsNormal_up, BasalvsHER2_up,
                        BasalvsNormLike_up,BasalvsNormal_up,
                        HER2vsNormLike_up ,  HER2vsNormal_up,  NormalLikevsNormal_up )))
all_down <- unique(sort(c(LumAvsLumB_down, LumAvsBasal_down,
                        LumAvsHER2_down, LumAvsNormLike_down ,
                        LumAvsNormal_down, LumBvsBasal_down, 
                        LumBvsHER2_down, LumBvsNormLike_down ,
                        LumBvsNormal_down, BasalvsHER2_down,
                        BasalvsNormLike_down,BasalvsNormal_down,
                        HER2vsNormLike_down ,  HER2vsNormal_down,  NormalLikevsNormal_down )))



# Comparing autophagy to all DE genes
autophagy_genes<- as.vector(read.table("../autopahagy_genes.txt", as.is = T, header = FALSE))$V1
sharedUP <- intersect(autophagy_genes,all_up)
sharedDOWN <- intersect(autophagy_genes,all_down)
print(paste0("Total number of UP genes : ", length(all_up)))
print(paste0("Autophagy genes in UP: ", length(sharedUP)))
print(paste0("Total number of DOWN genes : ", length(all_down)))
print(paste0("Autophagy genes in DOWN: ", length(sharedDOWN)))
autoUPandDOWN <- intersect (sharedUP, sharedDOWN)
print(paste0("Autophagy genes that are both UP and DOWN: ", length(autoUPandDOWN )))

#comparing autphagy to only DE genes vs normal 
up_norm <- as.vector(unlist(up_vs_normal))
down_norm <- as.vector(unlist(down_vs_normal))
sharedUP <- intersect(autophagy_genes,up_norm)
sharedDOWN <- intersect(autophagy_genes,down_norm)
print(paste0("Total number of UP genes : ", length(unique(up_norm))))
print(paste0("Autophagy genes in UP: ", length(sharedUP)))
print(paste0("Total number of DOWN genes : ", length(unique(down_norm))))
print(paste0("Autophagy genes in DOWN: ", length(sharedDOWN)))
autoUPandDOWN <- intersect (sharedUP, sharedDOWN)
print(paste0("Autophagy genes that are both UP and DOWN: ", length(autoUPandDOWN )))



## getting autophagy genes DE in each subtype
sharedWithAuto <-function (gene_list){
  autophagy_genes<- as.vector(read.table("../autopahagy_genes.txt", as.is = T, header = FALSE))$V1
  shared <- intersect(autophagy_genes,gene_list)
  print(paste0("Total number of genes in the list : ", length(gene_list)))
  print(paste0("Autophagy genes among these: ", length(shared)))
  return(shared)
  
}

auto_LumAvsLumB_up       <- sharedWithAuto(LumAvsLumB_up       )
auto_LumAvsBasal_up      <- sharedWithAuto(LumAvsBasal_up      )
auto_LumAvsHER2_up       <- sharedWithAuto(LumAvsHER2_up       )
auto_LumAvsNormLike_up   <- sharedWithAuto(LumAvsNormLike_up   )
auto_LumAvsNormal_up     <- sharedWithAuto(LumAvsNormal_up     )
auto_LumBvsBasal_up      <- sharedWithAuto(LumBvsBasal_up      )
auto_LumBvsHER2_up       <- sharedWithAuto(LumBvsHER2_up       )
auto_LumBvsNormLike_up   <- sharedWithAuto(LumBvsNormLike_up   )
auto_LumBvsNormal_up     <- sharedWithAuto(LumBvsNormal_up     )
auto_BasalvsLumA_up      <- sharedWithAuto(BasalvsLumA_up      )
auto_BasalvsLumB_up      <- sharedWithAuto(BasalvsLumB_up      )
auto_BasalvsHER2_up      <- sharedWithAuto(BasalvsHER2_up      )
auto_BasalvsNormLike_up  <- sharedWithAuto(BasalvsNormLike_up  )
auto_BasalvsNormal_up      <- sharedWithAuto(BasalvsNormal_up      )
auto_HER2vsNormLike_up     <- sharedWithAuto(HER2vsNormLike_up     )
auto_HER2vsNormal_up       <- sharedWithAuto(HER2vsNormal_up       )
auto_NormalLikevsNormal_up <- sharedWithAuto(NormalLikevsNormal_up )

auto_LumAvsLumB_down       <- sharedWithAuto(LumAvsLumB_down       )
auto_LumAvsBasal_down      <- sharedWithAuto(LumAvsBasal_down      )
auto_LumAvsHER2_down       <- sharedWithAuto(LumAvsHER2_down       )
auto_LumAvsNormLike_down   <- sharedWithAuto(LumAvsNormLike_down   )
auto_LumAvsNormal_down     <- sharedWithAuto(LumAvsNormal_down     )
auto_LumBvsBasal_down      <- sharedWithAuto(LumBvsBasal_down      )
auto_LumBvsHER2_down       <- sharedWithAuto(LumBvsHER2_down       )
auto_LumBvsNormLike_down   <- sharedWithAuto(LumBvsNormLike_down   )
auto_LumBvsNormal_down     <- sharedWithAuto(LumBvsNormal_down     )
auto_BasalvsLumA_down      <- sharedWithAuto(BasalvsLumA_down      )
auto_BasalvsLumB_down      <- sharedWithAuto(BasalvsLumB_down      )
auto_BasalvsHER2_down      <- sharedWithAuto(BasalvsHER2_down      )
auto_BasalvsNormLike_down  <- sharedWithAuto(BasalvsNormLike_down  )
auto_BasalvsNormal_down    <- sharedWithAuto(BasalvsNormal_down    )
auto_HER2vsNormLike_down   <- sharedWithAuto(HER2vsNormLike_down   )
auto_HER2vsNormal_down     <- sharedWithAuto(HER2vsNormal_down     )
auto_NormalLikevsNormal_down<- sharedWithAuto(NormalLikevsNormal_down)




auto_up_vs_normal <- list(auto_LumA = auto_LumAvsNormal_up,
                          auto_LumB = auto_LumBvsNormal_up,
                          auto_Basal = auto_BasalvsNormal_up,
                          auto_HER2 = auto_HER2vsNormal_up,
                          auto_NormLike = auto_NormalLikevsNormal_up)

auto_down_vs_normal <- list(auto_LumA = auto_LumAvsNormal_down,
                            auto_LumB = auto_LumBvsNormal_down,
                            auto_Basal = auto_BasalvsNormal_down,
                            auto_HER2 = auto_HER2vsNormal_down,
                            auto_NormLike = auto_NormalLikevsNormal_down)






draw5venn<-function (genes_list, title){
      grid.newpage()
      
      venn <- venn.diagram(genes_list, 
                           category.names = names(genes_list),
                           filename = NULL, lwd = 0.5,
                           main = title, main.cex=1.5,
                           fill=c("turquoise", "salmon", "blue", "yellow", "green"),
                           col = "transparent",
                           label.col =c("grey"),
                           alpha =0.5,
                           label.cex = 2,
                           cat.just=list(c(0.6,1) , c(0,-4) , c(1,0) , c(1,1) , c(0.8,-3)),
                           sub.cex = 5,
                           cat.cex=0.8, #group names
                           cat.fontfamily = 'sans',
                           cex=2)
      
      overlaps<-calculate.overlap(genes_list)
      print (overlaps)
      overlaps <- rev(overlaps)
      
      ovelap_list<-list()
      for (i in 1:length(overlaps)){
          venn[[i+10]]$label <- i #paste(overlaps[[i]], collapse = "\n")
          #print(length(overlaps[[i]]) )
          if (length(overlaps[[i]]) != 0){
            id<-paste0(i)
            ovelap_list[[id]]<-overlaps[[i]]
            print (i)
            print (overlaps[[i]])
            
          }
      }
      
      grid.draw(venn)
      
    
      venn <- venn.diagram(genes_list, 
                         category.names = names(genes_list),
                         filename = NULL, lwd = 0.5,
                         main = title, main.cex=1.5,
                         fill=c("turquoise", "salmon", "blue", "yellow", "green"),
                         col = "transparent",
                         label.col =c("black"),
                         alpha =0.01,
                         label.cex = 2,
                         cat.just=list(c(0.6,1) , c(0,-4) , c(1,0) , c(1,1) , c(0.8,-3)),
                         sub.cex = 5,
                         cat.cex=0.8, #group names
                         cat.fontfamily = 'sans',
                         cex=1)
      grid.draw(venn)
      return(ovelap_list)
}

# down genes vs normal
overlap_list<-draw5venn(auto_down_vs_normal, "Venn diagram Autophagy Down-regulated genes vs Normal")

overlap_list<-draw5venn(auto_up_vs_normal, "Venn diagram Autophagy Up-regulated genes vs Normal")

authopagy_function <- get(load("../autophagy_genes_lists/autophagy_functions.rda"))
head(authopagy_function)



for (name in names(overlap_list)) {
  print(paste0("Overlap ",name, ":"))
  for (gene in overlap_list[[name]]){
    print (authopagy_function[authopagy_function$genes==gene,], row.names = FALSE, col.names= FALSE)
  } 
  cat("\n")
  
}




draw4venn<-function (genes_list, title){
  grid.newpage()
  
  venn <- venn.diagram(genes_list, 
                       category.names = names(genes_list),
                       filename = NULL, lwd = 0.5,
                       main = title, main.cex=1.5,
                       fill=c("turquoise", "salmon", "blue", "purple"),
                       col = "transparent",
                       label.col =c("grey"),
                       alpha =0.5,
                       label.cex = 2,
                       #cat.just=list(c(0.6,1) , c(0,-4) , c(1,0) , c(1,1) , c(0.8,-3)),
                       sub.cex = 5,
                       cat.cex=0.8, #group names
                       cat.fontfamily = 'sans',
                       cex=2)
  
  
  overlaps<-calculate.overlap(genes_list)
  overlaps <- rev(overlaps)
  for (i in 1:length(overlaps)){
    venn[[i+8]]$label <- i #paste(overlaps[[i]], collapse = "\n")
    if (length(overlaps[[i]]) != 0){
      print (i)
      index <-paste0("a", i)
      #print (index)
      print (overlaps[[index]])
    }
  }
  
  grid.draw(venn)
  
  
  venn <- venn.diagram(genes_list, 
                       category.names = names(genes_list),
                       filename = NULL, lwd = 0.5,
                       main = title, main.cex=1.5,
                       fill=c("turquoise", "salmon", "blue", "purple"),
                       col = "transparent",
                       label.col =c("black"),
                       alpha =0.01,
                       label.cex = 2,
                       #cat.just=list(c(0.6,1) , c(0,-4) , c(1,0) , c(1,1) , c(0.8,-3)),
                       sub.cex = 5,
                       cat.cex=0.8, #group names
                       cat.fontfamily = 'sans',
                       cex=1)
  grid.draw(venn)
}
draw4venn(up_vs_normal[1:4], "Venn diagram Upregulated All VS Normal")


draw3venn<-function (genes_list, title){
  grid.newpage()
  
  venn <- venn.diagram(genes_list, 
                       category.names = names(genes_list),
                       filename = NULL, lwd = 0.5,
                       main = title, main.cex=1.5,
                       fill=c("turquoise", "salmon", "blue"),
                       col = "transparent",
                       label.col =c("grey"),
                       alpha =0.5,
                       label.cex = 2,
                       #cat.just=list(c(0.6,1) , c(0,-4) , c(1,0) , c(1,1) , c(0.8,-3)),
                       sub.cex = 5,
                       cat.cex=0.8, #group names
                       cat.col=c("white"),
                       cat.fontfamily = 'sans',
                       cex=2)
  
  overlaps<-calculate.overlap(genes_list)
  overlaps <- rev(overlaps)
  for (i in 1:length(overlaps)){
    venn[[i+6]]$label <- i #paste(overlaps[[i]], collapse = "\n")
   # if (length(overlaps[[i]]) != 0){
      print (i)
      index <-paste0("a", i)
      #print (index)
      print (length(overlaps[[index]]))
   #}
  }
  
  grid.draw(venn)
  
  
  venn <- venn.diagram(genes_list, 
                       category.names = names(genes_list),
                       filename = NULL, lwd = 0.5,
                       main = title, main.cex=1.5,
                       fill=c("turquoise", "salmon", "blue"),
                       col = "transparent",
                       label.col =c("black"),
                       alpha =0.01,
                       label.cex = 2,
                       #cat.just=list(c(0.6,1) , c(0,-4) , c(1,0) , c(1,1) , c(0.8,-3)),
                       sub.cex = 5,
                       cat.cex=0.8, #group names
                       cat.fontfamily = 'sans',
                       cex=1)
  grid.draw(venn)
}
draw3venn(up_vs_normal[c(3,4,5)], "Venn diagram Upregulated All VS Normal") #NB print out in this one does not fully work

draw2venn<-function (genes_list, title){
  grid.newpage()
  
  venn <- venn.diagram(genes_list, 
                       category.names = names(genes_list),
                       filename = NULL, lwd = 0.5,
                       main = title, main.cex=1.5,
                       fill=c( "salmon", "blue"),
                       col = "transparent",
                       label.col =c("grey"),
                       alpha =0.5,
                       label.cex = 2,
                       #cat.just=list(c(0.6,1) , c(0,-4) , c(1,0) , c(1,1) , c(0.8,-3)),
                       sub.cex = 5,
                       cat.cex=0.8, #group names
                       cat.fontfamily = 'sans',
                       cex=2)
  
  overlaps<-calculate.overlap(genes_list)
  overlaps <- rev(overlaps)
  for (i in 1:length(overlaps)){
    venn[[i+4]]$label <- i #paste(overlaps[[i]], collapse = "\n")
    if (length(overlaps[[i]]) != 0){
      print (i)
      #index a1 is the overlap
      if (i == 1 ){
        print (length(overlaps[[2]]) - length(overlaps[[1]]))
      } else if (i == 2 ){
        print (length(overlaps[[3]]) - length(overlaps[[1]]))
      }else {
        print (length(overlaps[[1]]))
      }
    }
  }
  
  grid.draw(venn)
  
  
  venn <- venn.diagram(genes_list, 
                       category.names = names(genes_list),
                       filename = NULL, lwd = 0.5,
                       main = title, main.cex=1.5,
                       fill=c( "salmon", "blue"),
                       col = "transparent",
                       label.col =c("black"),
                       alpha =0.01,
                       label.cex = 2,
                       #cat.just=list(c(0.6,1) , c(0,-4) , c(1,0) , c(1,1) , c(0.8,-3)),
                       sub.cex = 5,
                       cat.cex=0.8, #group names
                       cat.fontfamily = 'sans',
                       cex=1)
  grid.draw(venn)
} #NB print out in this one does not fully work
draw2venn(up_vs_normal[c(1,2)], "Venn diagram Upregulated All VS Normal")


## basal study
 
auto_up_basal <- list(BasalvsLumA = auto_BasalvsLumA_up ,
                      BasalvsLumB = auto_BasalvsLumB_up ,
                      BasalvsHER2 = auto_BasalvsHER2_up  ,  
                      BasalvsNormal = auto_BasalvsNormal_up,
                      BasalvsNormLike = auto_BasalvsNormLike_up  )
  
auto_down_basal <- list(BasalvsLumA = auto_BasalvsLumA_down ,
                      BasalvsLumB = auto_BasalvsLumB_down ,
                      BasalvsHER2 = auto_BasalvsHER2_down  ,  
                      BasalvsNormal = auto_BasalvsNormal_down ,
                      BasalvsNormLike = auto_BasalvsNormLike_down  )

  
draw5venn(auto_up_basal, "Venn diagram Autophagy Up-regulated genes in Basal")
draw5venn(auto_down_basal, "Venn diagram Autophagy Down-regulated genes in Basal")

draw2venn(auto_up_basal, "Venn diagram Autophagy Up-regulated genes in Basal")



x<-calculate.overlap(auto_up_basal)
x
