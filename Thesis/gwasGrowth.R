setwd("C:\\Users\\patri\\Dropbox\\UMCG\\Thesis\\Chapter 1")

gwas <- read.delim("./gwas_catalog_v1.0.1-associations_e84_r2016-07-17.tsv", stringsAsFactors = F)

str(gwas)

gwas$DATE


library(zoo)
gwas$DATEQ <- as.factor(as.yearqtr(as.Date(gwas$DATE, format = "%Y-%m-%d")))






qCount <- table(factor(gwas$DATEQ[gwas$P.VALUE <= 5e-8 & gwas$DATEQ != "NA QNA" & gwas$DATEQ != "2015 Q4"]))

qCount2 <- qCount

for(i in 2:length(qCount2)){
  qCount2[i] <- qCount2[i] + qCount2[i-1]
}

str(qCount2)

png("./gwasGrowth.png")
barplot(qCount2)
dev.off()