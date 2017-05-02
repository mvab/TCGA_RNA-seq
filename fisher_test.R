
n_matrix <- matrix(data=0, 6, 6,
                   dimnames = list(c("LumA", "LumB", "Basal", "HER2", "Norm-like", "Normal"),
                                  c("LumA", "LumB", "Basal", "HER2", "Norm-like", "Normal")))
class(n_matrix)<-"numeric"            
k_matrix <- matrix(data=0, 6, 6,
                   dimnames = list(c("LumA", "LumB", "Basal", "HER2", "Norm-like", "Normal"),
                                   c("LumA", "LumB", "Basal", "HER2", "Norm-like", "Normal")))
class(k_matrix)<-"numeric"  

pval_matrix <- matrix(data=0, 6, 6,
                   dimnames = list(c("LumA", "LumB", "Basal", "HER2", "Norm-like", "Normal"),
                                   c("LumA", "LumB", "Basal", "HER2", "Norm-like", "Normal")))
class(pval_matrix)<-"numeric"  

odds_matrix <- matrix(data=0, 6, 6,
                   dimnames = list(c("LumA", "LumB", "Basal", "HER2", "Norm-like", "Normal"),
                                   c("LumA", "LumB", "Basal", "HER2", "Norm-like", "Normal")))
class(odds_matrix)<-"numeric"  





getIndexOfPair<-function(row, col, matrix){
  row_ind<-match(row, rownames(matrix))
  col_ind<-match(col, colnames(matrix))
  return (list(row=row_ind, col=col_ind))
}
index<-getIndexOfPair("LumB", "HER2", n_matrix)

n_matrix["LumA", "LumB"]<-86
n_matrix["LumA", "Basal"]<-1219
n_matrix["LumA", "HER2"]<-448
n_matrix["LumA", "Norm-like"]<-360
n_matrix["LumA", "Normal"]<-1402
n_matrix["LumB", "Basal"]<-1082
n_matrix["LumB", "HER2"]<-372
n_matrix["LumB", "Norm-like"]<-902
n_matrix["LumB", "Normal"]<-2085
n_matrix["Basal", "HER2"]<-720
n_matrix["Basal", "Norm-like"]<-814
n_matrix["Basal", "Normal"]<-2002
n_matrix["HER2", "Norm-like"]<-647
n_matrix["HER2", "Normal"]<-2016
n_matrix["Norm-like", "Normal"]<-486

k_matrix["LumA", "LumB"]<-2
k_matrix["LumA", "Basal"]<-44
k_matrix["LumA", "HER2"]<-16
k_matrix["LumA", "Norm-like"]<-10
k_matrix["LumA", "Normal"]<-45
k_matrix["LumB", "Basal"]<-41
k_matrix["LumB", "HER2"]<-15
k_matrix["LumB", "Norm-like"]<-36
k_matrix["LumB", "Normal"]<-72
k_matrix["Basal", "HER2"]<-25
k_matrix["Basal", "Norm-like"]<-34
k_matrix["Basal", "Normal"]<-73
k_matrix["HER2", "Norm-like"]<-31
k_matrix["HER2", "Normal"]<-79
k_matrix["Norm-like", "Normal"]<-16
  

k_matrix
n_matrix

for (i in 1:6){
  for (j in 1:6){
    x<-n_matrix[i,j]
    y<-k_matrix[i,j]
    
    N <- 15725 # The total number of genes :17372
    K <- 1087 # The number of gene belonging to a gene famility---- AUTOPHAGY
    n <- x # The list of interesting genes - DEGs in this comparison
    k <- y # The number of gene family members in the interesting genes : autophagy in DEGs of comparison
    
    m <- matrix(c(k, K - k, n - k, N - K - n + k),
                2, 2, dimnames = list(c("DE both", "no DE"),
                                      c("Autophagy", "Rest")))
    #Fisher's exact test is included in base R:
    res <- fisher.test(x=m, alternative="two.sided")
    pval_matrix[i,j]<-res$p.value
    odds_matrix[i,j]<-res$estimate
  }
}
options("scipen"=10)
options()$scipen

pval_matrix
odds_matrix #magnitude of difference
