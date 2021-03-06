setwd("D:\\UMCG\\Genetica\\Projects\\Covid19\\pgs\\pgs_correlation_output")

library(readr)

## Functions 

corToZ <- function(r, df){
  t = sqrt(df) * abs(r) / sqrt(1 - (r*r))
  p = min(pt(t, df, log.p = T), pt(t, df, lower.tail=FALSE, log.p = T))
  z = qnorm(p, log.p = T)
  if(r > 0){
    z <- z * -1
  }
  return(z)
}


zToCor <- function(z, df){
  t <- qt(pnorm(-abs(z),log.p = T), df, log.p=T)
  r <- t/sqrt(df+t^2)
  if(z > 0){
    r <- r * -1
  }
  return(r)
}


## Main

# ! make sure r files and n files are in the same order.
rFiles <- c("allParti_correctResiduals_cyto/Neuroticism/covid_export_correlations_Neuroticism_correlations_allParti_correctResiduals_cyto_03-03-2021.txt", "allParti_correctResiduals_ugli/Neuroticism/covid_export_correlations_Neuroticism_correlations_allParti_correctResiduals_ugli_03-03-2021.txt")
nFiles <- c("allParti_correctResiduals_cyto/Neuroticism/covid_export_correlations_Neuroticism_n_values_allParti_correctResiduals_cyto_03-03-2021.txt", "allParti_correctResiduals_ugli/Neuroticism/covid_export_correlations_Neuroticism_n_values_allParti_correctResiduals_ugli_03-03-2021.txt")
outputFile <- "test.txt"


rMatrices <- lapply(rFiles, function(f){
  table_tmp <- read_delim(f, delim = "\t", quote = "")
  table <- as.matrix(table_tmp[,-c(1:4)])
  rownames(table) <- table_tmp[,1][[1]]
  rm(table_tmp)
  return(table)
})

nMatrices <- lapply(nFiles, function(f){
  table_tmp <- read_delim(f, delim = "\t", quote = "")
  table <- as.matrix(table_tmp[,-1])
  rownames(table) <- table_tmp[,1][[1]]
  rm(table_tmp)
  return(table)
})

if(length(rFiles != length(nFiles))){
  stop("number r files not equals n files")
}

expectedDim = dim(rMatrices[[1]])
expectedRows = rownames(rMatrices[[1]])
expectedCols = colnames(rMatrices[[1]])


for(i in 1:length(rFiles)){
  if(!all(dim(rMatrices[[i]]) == expectedDim) | !all(rownames(rMatrices[[i]]) == expectedRows) | !all(colnames(rMatrices[[i]]) == expectedCols)){
    stop("Matrix not equal")
  }
}

for(i in 1:length(rFiles)){
  if(!all(dim(nMatrices[[i]]) == expectedDim) | !all(rownames(nMatrices[[i]]) == expectedRows) | !all(colnames(nMatrices[[i]]) == expectedCols)){
    stop("Matrix not equal")
  }
}

zMatrices <- lapply(1:length(rFiles), function(i){
  
  rMatrix <- rMatrices[[i]]
  nMatrix <- nMatrices[[i]]
  zMatrix <- rMatrices[[i]]
  
  #if(!all(dim(rMatrix) == dim(nMatrix)) | !all(rownames(rMatrix) == rownames(nMatrix)) | !all(colnames(rMatrix) == colnames(nMatrix))){
  #  stop("Matrix not equal")
  #}
  
  for(r in 1:nrow(rMatrix)){
    for(c in 1:ncol(rMatrix)){
      if(nMatrix[r,c] <= 2 | is.na(rMatrix[r,c])){
        zMatrix[r,c] <- 0
        nMatrices[[i]][r,c] <<- 0#force to zero so that z-score is not used
      } else {
        zMatrix[r,c] <- corToZ(rMatrix[r,c], nMatrix[r,c]-2)
      }
      
      
    }
  }
  
  return(zMatrix)
  
})


metaRMatrix <- rMatrices[[1]]

for(r in 1:nrow(rMatrices[[1]])){
  for(c in 1:ncol(rMatrices[[1]])){
    
    a=0;#sum(weight_i * Z_i)
    b=0;#sum(weigth_i^2)
    sumN = 0;
    
    for(i in 1:length(rFiles)){
      n = nMatrices[[i]][r,c]
      sumN = sumN + n
      a = a + (n * zMatrices[[i]][r,c])
      b = b + (n*n)
    }
    
    if(sumN <= 2){
      metaRMatrix[r,c] <- NA
    } else {
      z <- a / sqrt(b)
      
      metaRMatrix[r,c] <- zToCor(z, sumN - 2) 
    }
    
  }
}

write.table(metaRMatrix, file = outputFile, sep = "\t", quote = F, col.names = NA)
