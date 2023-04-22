corr_check <- function(Dataset, threshold){
  matriz_cor <- cor(Dataset)
  matriz_cor
  
  for (i in 1:nrow(matriz_cor)){
    correlations <-  which((abs(matriz_cor[i,i:ncol(matriz_cor)]) > threshold) & (matriz_cor[i,i:ncol(matriz_cor)] != 1))
    
    if(length(correlations)> 0){
      lapply(correlations,FUN =  function(x) (cat(paste(colnames(Dataset)[i], "with",colnames(Dataset)[x]), "\n")))
      
    }
  }
}