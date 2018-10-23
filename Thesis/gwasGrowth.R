setwd("C:\\Users\\patri\\Dropbox\\UMCG\\Thesis\\")

gwas <- read.delim("./Chapter 1/gwas_catalog_v1.0.2-associations_e93_r2018-08-04.tsv", stringsAsFactors = F)

str(gwas)

gwas$DATE


library(zoo)
gwas$DATEQ <- as.factor(as.yearqtr(as.Date(gwas$DATE, format = "%Y-%m-%d")))





qCount <- table(factor(gwas$DATEQ[gwas$P.VALUE <= 5e-8 & gwas$DATEQ != "NA QNA"]))# & gwas$DATEQ != "2015 Q4"

qCount2 <- qCount

for(i in 2:length(qCount2)){
  qCount2[i] <- qCount2[i] + qCount2[i-1]
}

str(qCount2)

png("./Chapter 1/gwasGrowth.png")
barplot(qCount2)
dev.off()

pdf("./Chapter 1/gwasGrowth.pdf", height = 5, width = 7)
par(xpd = NA, las = 2)
barplot(qCount2)
dev.off()

qCount2
