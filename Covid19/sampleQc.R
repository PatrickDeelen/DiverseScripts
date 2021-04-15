#remoter::server(verbose = T, port = 55556, password = "laberkak", sync = T)

remoter::client("localhost", port = 55556, password = "laberkak")

library(readr)
pheno <- read_delim("/groups/umcg-lifelines/tmp01/projects/ov20_0554/analysis/risky_behaviour/PRS_correlation/questioniare_subset_participants_with_genome_data/questionaire_df_subset_participants_with_genome_data_01-03-2021.txt", delim = "\t", quote = "")
dim(pheno)
pheno2 <- as.data.frame(pheno)
row.names(pheno2) <- pheno2[,1]

#saveRDS(pheno2, "/groups/umcg-lifelines/tmp01/projects/ov20_0554/analysis/risky_behaviour/PRS_correlation/questioniare_subset_participants_with_genome_data/questionaire_df_subset_participants_with_genome_data_01-03-2021.rds")

pheno2 <- readRDS("/groups/umcg-lifelines/tmp01/projects/ov20_0554/analysis/risky_behaviour/PRS_correlation/questioniare_subset_participants_with_genome_data/questionaire_df_subset_participants_with_genome_data_01-03-2021.rds")
dim(x)

c("gender_recent", "age_recent", "age2_recent", "chronic_recent", "household_recent", "have_childs_at_home_recent") %in% colnames(pheno2)


sum(pheno2==8888, na.rm =T)

dim(pheno2)

qOverview <- read.delim("/groups/umcg-lifelines/tmp01/projects/ov20_0554/analysis/pgs_correlations/questionair_time_overview_nl.txt", stringsAsFactors = F, row.names = 1)
vls <- colnames(qOverview)[-c(21,22)]

vl = vls[1]


missingVraagList <- sapply(vls, function(vl){
  
  cat("----", vl, "\n")
  
  vlq <- qOverview[,vl]
  vlq <- vlq[vlq!=""]
  vlq <- vlq[vlq%in%colnames(pheno2)]
  
  cat(length(vlq), "\n")
  
  
  vlDate <- qOverview["responsdatum covid-vragenlijst",vl]
  
  #participants of this vl
  vlp <- !is.na(pheno2[,vlDate])
  
  print(table(vlp))
  
  #calculate percentage only among participants this vl
  vlMissing <- apply(pheno2[vlp,vlq], 2, function(x){sum(is.na(x)*100/length(x))})
  
  print(str(vlMissing))
  
  return(vlMissing)
  
})

missingVraagList[["X14.0"]]

range(sapply(missingVraagList, length))

pdf("/groups/umcg-lifelines/tmp01/projects/ov20_0554/missingPerVl.pdf")
sapply(vls, function(vl){
  hist(missingVraagList[[vl]], main = vl, breaks = 100)
})
dev.off()

qForSampleQc <- do.call("c", sapply(missingVraagList, function(x){names(x)[x<=5]}))

missing <- sapply(vls, function(vl){
  
  vlq <- qOverview[,vl]
  vlq <- vlq[vlq!=""]
  vlq <- vlq[vlq%in%qForSampleQc]
  
  vlMissing <- apply(pheno2[,vlq], 1, function(x){sum(is.na(x)*100/length(x))})
  
  return(vlMissing)
  
})

str(missing)

hist(missing, breaks = 100, ylim = c(0,1000))
dev.off()

ppvl <- apply(missing, 2, function(x){
  sum(x <= 5)
})
barplot(ppvl)
dev.off()

qpp <- apply(missing, 1, function(x){
  sum(x <= 5)
})
table(qpp)



colnames(missing)[1:12]

participatedFirstHalf <- apply(missing[,c(1:12)], 1, function(x){
  any(x <= 5)
})

participatedSecondHalf <- apply(missing[,c(13:19)], 1, function(x){
  any(x <= 5)
})

sum(participatedFirstHalf)
sum(participatedSecondHalf)


participatedBothHalf <- participatedFirstHalf & participatedSecondHalf
sum(participatedBothHalf)



qpp <- apply(missing[participatedBothHalf,], 1, function(x){
  sum(x <= 5)
})
table(qpp)


median(qpp)

inclusionPerVl <- missing[participatedBothHalf,] <= 5
str(inclusionPerVl)
#write.table(inclusionPerVl, file = "/groups/umcg-lifelines/tmp01/projects/ov20_0554/analysis/risky_behaviour/PRS_correlation/inclusionPerVl.txt", quote = F, sep = "\t", col.names = NA)

x <- apply(missing[names(qpp)[qpp==0],], 1, function(x){sum(x < 100)})
table(x)
names(x)[x>10]

p <- names(qpp)[qpp>=1]

missing2 <- missing[p,] 

y <- apply(missing2[,-c(1:4)], 1, function(x){
  sum(x > 5 & x<100)
})

barplot(table(y))
  dev.off()

rpng(width = 1000, height = 1000)
barplot(table(qpp))
dev.off()



validationSet <- apply(missing[,c("X4.0","X9.0", "X14.0", "X17.0")], 1, function(x){
  all(x <= 5)
})
table(validationSet)

write.table(names(validationSet)[validationSet], file = "/groups/umcg-lifelines/tmp01/projects/ov20_0554/analysis/risky_behaviour/PRS_correlation/validationSamples.txt", quote = F, row.names = F, col.names = F);

vl = "X1.0"

pdf("/groups/umcg-lifelines/tmp01/projects/ov20_0554/pcPerWeek.pdf")

sapply(vls, function(vl){


  
vlq <- qOverview[,vl]
vlq <- vlq[vlq!=""]
vlq <- vlq[vlq%in%qForSampleQc]


vlClean <- pheno2[missing[,vl] <= 5, vlq]

types <- sapply(vlClean, class)
vlClean <- as.matrix(vlClean[,types == "numeric" | types == "logical"])


vlCleanColMean <- apply(vlClean, 2, mean, na.rm = T)
for(c in 1:ncol(vlClean)){
  vlClean[is.na(vlClean[,c ]),c] <- vlCleanColMean[c]
}

vlClean <- vlClean[,apply(vlClean, 2, sd)>0]
dim(vlClean)
vlClean <- scale(vlClean)

pcaRes <- prcomp(t(vlClean), center = T, scale. = T)

#plot(pcaRes)
#dev.off()

png(paste0("/groups/umcg-lifelines/tmp01/projects/ov20_0554/pcPerWeek/", vl, ".png"))
plot(pcaRes$rotation[,c(1,2)], bg = adjustcolor("dodgerblue2", alpha.f = 0.1), pch = 21, col=adjustcolor("dodgerblue2", alpha.f = 0.1))
(sd1 <- sd(pcaRes$rotation[,1]))
abline(v=sd1*4, col = "firebrick")
abline(v=-sd1*4, col = "firebrick")
(sd2 <- sd(pcaRes$rotation[,2]))
abline(h=sd2*4, col = "firebrick")
abline(h=-sd2*4, col = "firebrick")
grDevices::dev.off()
})
