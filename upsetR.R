library(UpSetR)
library(ggplot2)
library(grid)
library(plyr)
movies <- read.csv(system.file("extdata", "movies.csv", package = "UpSetR"), 
                   header = T, sep = ";")
dim(movies)




#try my
expressionInput <- c(`stage1&LumA` = 96,`stage2&LumA` = 219,`stage3&LumA` =96 ,`stage4&LumA` =4 ,
                     `stage1&LumB` = 22, `stage2&LumB` = 104,`stage3&LumB` = 55 , `stage4&LumB` = 4 ,
                     `stage1&Basal` = 22 ,`stage2&Basal` = 106,`stage3&Basal` = 21, `stage4&Basal` = 2,
                     `stage1&HER2` = 6 ,`stage2&HER2` = 47,`stage3&HER2` = 19, `stage4&HER2` = 2,
                     `stage1&NormLike` = 5 ,`stage2&NormLike` = 13,`stage3&NormLike` = 8,
                     `stage1&Lobular` = 14,`stage2&Lobular` = 78,`stage3&Lobular` =49 ,
                     `stage1&Ductal` = 118, `stage2&Ductal` = 368,`stage3&Ductal` = 129 , `stage4&Ductal` = 11 ,
                     `stage1&LobDuct` = 5 ,`stage2&LobDuct` = 11,`stage3&LobDuct` = 7, 
                     `stage1&DuctOthers` = 3 ,`stage2&DuctOthers` = 8,`stage3&DuctOthers` = 4, 
                     `stage1&Metaplast` = 2 ,`stage2&Metaplast` = 4,`stage3&Metaplast` = 1,
                     `stage1&Mucinous` = 3 ,`stage2&Mucinous` = 5,`stage3&Mucinous` = 4)


#metadata

sets<- c("Basal", "HER2", "LumA", "LumB" ,"NormLike" , "LobDuct",  "Ductal", "DuctOthers", "Lobular", "Metaplast", "Mucinous",
         "stage1"   , "stage2", "stage3", "stage4"  )
groups<-c(rep("PAM50",5), rep("Morph", 6), rep("stage",4))

metadata <- as.data.frame(cbind(sets, groups))
names(metadata) <- c("sets", "groups")
View(metadata)


upset(fromExpression(expressionInput), 
      order.by = "freq", nsets = 38, keep.order = TRUE,
      sets = c("NormLike","HER2", "Basal",  "LumB", "LumA",
              "stage4", "stage3",  "stage2", "stage1",
             "Mucinous", "Metaplast", "DuctOthers", "LobDuct", "Ductal", "Lobular"), 
      mainbar.y.label = "Group intersections sizes", sets.x.label = "Samples per group" ,
      set.metadata = list(data = metadata, 
                         plots = list(list(type = "matrix_rows", column = "groups", #assign = 5,
                         colors = c(PAM50 = "green", Morph = "navy",  stage = "purple"),alpha = 0.4))))
      



stop()


















upset(movies, main.bar.color = "black",
      queries = list(list(query = intersects, params = list("Drama"), active = T)), 
      attribute.plots = list(gridrows = 50))
upset(movies, main.bar.color = "black",
      queries = list(list(query = intersects)), 
      attribute.plots = list(gridrows = 50))

upset(movies, main.bar.color = "black",
  queries = list(list(query = intersects,  params = list("Drama"), color = "red", active = T),
                 list(query = intersects,  params = list("Action", "Drama"), active = T), 
                 list(query = intersects,  params = list("Drama", "Comedy", "Action"), color = "orange", active = T)), 
  attribute.plots = list(gridrows = 20, query.legend = "bottom"))



setwd("~/Bioinformatics MSc UCPH/0_MasterThesis/TCGAbiolinks/CBL_scripts/data")
design<-read.csv("testdesign.csv",header = T, sep = ",")

design2<-read.csv(system.file("testdesign.json", "testdesign.csv", package = "UpSetR"), 
         header = T, sep = ",")

View(design)

upset(design, main.bar.color = "black",
      queries = list(list(query = intersects, params = list("BasalLike"), active = T)), 
      attribute.plots = list(gridrows = 20))
