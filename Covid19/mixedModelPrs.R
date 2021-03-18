#remoter::server(verbose = T, port = 55556, password = "laberkak", sync = T)

remoter::client("localhost", port = 55556, password = "laberkak")

setwd("/groups/umcg-lifelines/tmp01/projects/ov20_0554/analysis/pgs_correlations/")

library("nlme")
library(heatmap3)
#library("ordinal")
library("GLMMadaptive")
library(readr)


if(FALSE){
  #Run once to conver pheno data to RDS format
  pheno <- read_delim("/groups/umcg-lifelines/tmp01/projects/ov20_0554/analysis/risky_behaviour/PRS_correlation/questtionnaire_data_subset_genome_filter_correct_filled_questionnaires/covid_export_questionnaire_1-17_questionnaire_filter_genome_correct_filled_18-03-2021.txt", delim = "\t", quote = "", guess_max = 100000)
  dim(pheno)
  pheno2 <- as.data.frame(pheno)
  row.names(pheno2) <- pheno2[,1]
  
  saveRDS(pheno2, "/groups/umcg-lifelines/tmp01/projects/ov20_0554/analysis/risky_behaviour/PRS_correlation/questtionnaire_data_subset_genome_filter_correct_filled_questionnaires/covid_export_questionnaire_1-17_questionnaire_filter_genome_correct_filled_18-03-2021.rds")
  
}


#Constants
startdate <- as.Date("30/03/2020","%d/%m/%Y")
confounders <- c("gender_recent", "age_recent", "age2_recent", "chronic_recent", "household_recent", "have_childs_at_home_recent")


##load pheno and prs
pheno2 <- readRDS("/groups/umcg-lifelines/tmp01/projects/ov20_0554/analysis/risky_behaviour/PRS_correlation/questtionnaire_data_subset_genome_filter_correct_filled_questionnaires/covid_export_questionnaire_1-17_questionnaire_filter_genome_correct_filled_18-03-2021.rds")

if(!all(confounders %in% colnames(pheno2))){
  stop("Not all confounders found")
}

prsGsa <- read.delim("/groups/umcg-lifelines/tmp01/projects/ov20_0554/analysis/risky_behaviour/PRS_correlation/input_PGS_data_ugli/PGS_combined_ugli_26-02-2021.txt", stringsAsFactors = F)
prsCyto <- read.delim("/groups/umcg-lifelines/tmp01/projects/ov20_0554/analysis/risky_behaviour/PRS_correlation/input_PGS_data_cyto_v2_duplicate_removed/PGS_combined_cyto_duplicate_from_ugli_removed_18-03-2021.txt", stringsAsFactors = F)

if(!all(colnames(prsGsa) == colnames(prsCyto))){
  stop("Colnames must be equal")
}
if(!all(!row.names(prsGsa$PROJECT_PSEUDO_ID) %in% row.names(prsCyto$PROJECT_PSEUDO_ID))){
  stop("Overlapping samples")
}

prs <- rbind(prsGsa, prsCyto)

colnames(prs)[colnames(prs) == "BMI"] <- "BMI_gwas"

if(!all(pheno2$PROJECT_PSEUDO_ID %in% prs$PROJECT_PSEUDO_ID)){
  stop("Not all pheno have genetics")#easly solved but code makes this assumtion
}

pheno2$array <- factor(as.numeric(pheno2$PROJECT_PSEUDO_ID %in% prsGsa$PROJECT_PSEUDO_ID), levels = 0:1, labels = c("Cyto", "Gsa"))

pheno3 <- merge(pheno2, prs, by = "PROJECT_PSEUDO_ID")


totalPart <- nrow(pheno3)

#na col to use in reshapre for missing questions 
pheno3$naCOl <- NA


##Load and format questions meta data
qOverview <- as.matrix(read.delim("/groups/umcg-lifelines/tmp01/projects/ov20_0554/analysis/pgs_correlations/quest_overview_nl_new_quest17_codes.txt", stringsAsFactors = F, row.names = 1))
vls <- colnames(qOverview)[-c(20,21)]

# mask questions with >75% missing
qOverview[,-c(20,21)] <- apply(qOverview[,-c(20,21)], 1:2, function(x){
  if(x==""){
    return ("")
  } else if(!x %in% colnames(pheno3)) {
    return("")
  } else {
    if(sum(is.na(pheno3[,x])) >= (totalPart/4)){
      return (x)
    } else {
      return ("")
    }
  }
})

qOverview["responsdatum covid-vragenlijst",]
sum(qOverview == "")


# select questions with 4 repeats spread out over first and last half
qPassQc <- sapply(rownames(qOverview), function(q){
  
  qsFirstHalf <- qOverview[q,c(1:12)]
  qsSecondHalf <- qOverview[q,c(13:19)]
  return(sum(qsFirstHalf!="") >= 2 && sum(qsSecondHalf!="") >= 2)
  
})
table(qPassQc)

qOverview2 <- qOverview[qPassQc,-c(20,21)]
dim(qOverview2)


qNameMap <- data.frame(orginal = row.names(qOverview2), new = make.names(row.names(qOverview2), unique = T), stringsAsFactors = F)
row.names(qNameMap) <- row.names(qOverview2)


qList <- lapply(qNameMap[,1], function(q){
  qs <- qOverview2[q,]
  qs[qs==""]="naCOl"
  qs[!qs %in% colnames(pheno3)] <- "naCOl"
  return(qs)
})
names(qList) <- qNameMap[,2]
qList[[qNameMap["responsdatum covid-vragenlijst",2]]]
qList[[qNameMap["als u moet kiezen, denkt u zelf dat u een coronavirus/covid-19 infectie hebt (gehad)?",2]]]


## Reshape to long format and clean some variables

if(any(names(qList) %in% colnames(pheno3))){
  stop("Column name clash after reshape")
}

vragenLong <- reshape(pheno3, direction = "long", idvar = "PROJECT_PSEUDO_ID", varying = qList, v.names = names(qList), times = vls, timevar = "vl")
vragenLong$vl2 <- as.numeric(factor(vragenLong$vl, levels = vls, ordered = T))
vragenLong$vl3 <- factor(vragenLong$vl, levels = vls, ordered = F)
vragenLong$days <- as.numeric(difftime(vragenLong[,qNameMap["responsdatum covid-vragenlijst",2]], startdate ,units="days"))
vragenLong$days2 <- vragenLong$days*vragenLong$days


vragenLong$gender_recent <- factor(vragenLong$gender_recent, levels = 0:1, labels = c("female","male"))
vragenLong$household_recent <- factor(vragenLong$household_recent, levels = 0:1, labels = c("single-person household","multi-person household"))
vragenLong$have_childs_at_home_recent <- factor(vragenLong$have_childs_at_home_recent, levels = 0:1, labels = c("No childeren at home","Childeren at home"))
vragenLong$chronic_recent <- factor(vragenLong$chronic_recent, levels = 0:1, labels = c("Healthy","Chronic disease"))


hist(vragenLong$days, breaks = 330)
dev.off()

str(vragenLong)


## Run models

qLoop <- as.list(qNameMap[,2])
names(qLoop) <- qNameMap[,1]

q=qLoop[[1]]
q<-qNameMap["hoe waardeert u uw kwaliteit van leven over de afgelopen 7 dagen?",2]
zScoreList <- lapply(qLoop, function(q){
  zScores = tryCatch({
    fixedModel <- as.formula(paste(q, "~((gender_recent+age_recent+age2_recent+household_recent+have_childs_at_home_recent+chronic_recent +", paste(colnames(prs)[-1], collapse = " + ") ,")*days*array + days2) "))
    randomModel <- as.formula("~1|PROJECT_PSEUDO_ID")
    res <-  lme(fixed = fixedModel, random=randomModel, data= vragenLong[,c("PROJECT_PSEUDO_ID", q,colnames(prs)[-1],"gender_recent","age_recent","age2_recent","household_recent","have_childs_at_home_recent","chronic_recent", "days", "days2", "array" )],na.action=na.omit)
    tTable <- summary(res)$tTable
    zScores <- qnorm((tTable[,"p-value"]/2)) 
    zScores[is.infinite(zScores)] <- -30
    zScores[tTable[,"Value"] > 0] <- zScores[tTable[,"Value"] > 0] * -1
    return(zScores)
  }, error = function(e){print(q); return(NULL)})
  return(zScores)
})

str(zScoreList)
zScoreList2 <- zScoreList[!sapply(zScoreList, is.null)]

# combine into z-score matrix excluding intercept
zscores <- do.call("rbind", zScoreList2)
zscores <- zscores[,colnames(zscores) != "(Intercept)"]

write.table(zscores, file = "zscoreMatrix.txt", sep = "\t", quote = F, col.names = NA)


#below is testingground

sum(is.na(vragenLong[!is.na(vragenLong[,q]),c(colnames(prs)[-1])]))



cor.test(fitted(res), vragenLong[!is.na(vragenLong[,q]),q])

  write.table(summary(res)$tTable, file = "tmp3.txt", sep = "\t", quote = F, col.names = NA)

dim(zscores)

pdf("mixedModelZscores.pdf", width = 50, height = 50)
heatmap3(zscores, balanceColor = T, scale = "none", margins = c(20,20))
dev.off()

x <- zScoreList[[1]]
x[order(abs(x))]

q <- qLoop[[1]]

zScores[order(abs(zScores))]

###### Below is old testing stuff

q<-qNameMap["hoe waardeert u uw kwaliteit van leven over de afgelopen 7 dagen?",2]

question = qNameMap["op hoeveel momenten van de dag eet u iets?",2]
gwas = "Neuroticism"
fixedModel <- as.formula(paste(question, "~((gender_recent+age_recent+age2_recent+household_recent+have_childs_at_home_recent+chronic_recent +", paste(colnames(prs)[-1], collapse = " + ") ,")*days)"))
randomModel <- as.formula("~1|PROJECT_PSEUDO_ID")

res <-  lme(fixed = fixedModel, random=randomModel, data= vragenLong[,c("PROJECT_PSEUDO_ID", question,colnames(prs)[-1],"gender_recent","age_recent","age2_recent","household_recent","have_childs_at_home_recent","chronic_recent", "days" )],na.action=na.omit)
(x <- summary(res))
x$tTable

str(vragenLong[,question])





fixedModel <- as.formula(paste(q, "~(gender_recent+age_recent+age2_recent+household_recent+have_childs_at_home_recent+chronic_recent)+days+days2 "))
randomModel <- as.formula("~1|PROJECT_PSEUDO_ID")
res <-  lme(fixed = fixedModel, random=randomModel, data= vragenLong[,c("PROJECT_PSEUDO_ID", q,colnames(prs)[-1],"gender_recent","age_recent","age2_recent","household_recent","have_childs_at_home_recent","chronic_recent", "days", "days2" )],na.action=na.omit)
str(res$residuals)

summary(res)

x <- vragenLong[rownames(res$residuals),c("PROJECT_PSEUDO_ID", q,colnames(prs)[-1],"gender_recent","age_recent","age2_recent","household_recent","have_childs_at_home_recent","chronic_recent", "days", "days2" )]
str(x)

x$residuals <- res$residuals


res2 <-  lme(fixed = fixedModel, random=randomModel, data= vragenLong[,c("PROJECT_PSEUDO_ID", q,colnames(prs)[-1],"gender_recent","age_recent","age2_recent","household_recent","have_childs_at_home_recent","chronic_recent", "days", "days2" )],na.action=na.omit)
tTable <- summary(res2)$tTable





fixedModel <- as.formula(paste(question, "~(", gwas ,"*days)"))
randomModel <- as.formula("~1|PROJECT_PSEUDO_ID")
oridnalFit <-  mixed_model(fixed = fixedModel, random=randomModel, data= vragenLong[,c("PROJECT_PSEUDO_ID", question,gwas,"gender_recent","age_recent","age2_recent","household_recent","have_childs_at_home_recent","chronic_recent", "days" )],na.action=na.omit, family = poisson(), iter_EM = 0)
summary(oridnalFit)


#vragenLong[,qNameMap["hoe waardeert u uw kwaliteit van leven over de afgelopen 7 dagen?",2]] <- factor(vragenLong[,qNameMap["hoe waardeert u uw kwaliteit van leven over de afgelopen 7 dagen?",2]], levels = 1:10,  ordered = T)
model <- as.formula(paste(question, "~((gender+age_recent+age2_recent+household_recent+have_childs_at_home_recent+chronic_recent +", gwas ,")*days) + (1|PROJECT_PSEUDO_ID)"))
ordinalFit <- clmm(model, data = vragenLong,na.action=na.omit)

res <- lme(fixed=cijfer~((COVID.19.susceptibility)*days), random=~1|PROJECT_PSEUDO_ID, data= vragenLong,na.action=na.omit)
(x <- summary(res))
x$tTable


res <- lme(fixed=cijfer ~ Neuroticism*days, random=~1|PROJECT_PSEUDO_ID, data= vragenLong[vragenLong$vl>=13,],na.action=na.omit)
(x <- summary(res))

hist(vragenLong$cijfer, breaks = 10)
dev.off()

res <- lme(fixed=cijfer~((gender_recent+age_recent+age2_recent+household_recent+have_childs_at_home_recent + Neuroticism)*days), random=~1|PROJECT_PSEUDO_ID, data= vragenLong,na.action=na.omit)

str(res)

hist(res$residuals)
hist(res$fitted)


plot(ranef(res))
dev.off()

plot(res)
dev.off()

boxplot(pheno3$Neuroticism~pheno3$chronic_recent)
dev.off()
t.test(pheno3$Neuroticism ~pheno3$chronic_recent)

meanCijfer <- sapply(vls, function(vl){
  mean(vragenLong[vragenLong$vl==vl,"cijfer"], na.rm = T)
})

meanCijferHigh <- sapply(vls, function(vl){
  mean(vragenLong[vragenLong$vl==vl & vragenLong$Neuroticism >= 0.2,"cijfer"], na.rm = T)
})

meanCijferLow <- sapply(vls, function(vl){
  mean(vragenLong[vragenLong$vl==vl & vragenLong$Neuroticism < -0.2,"cijfer"], na.rm = T)
})


barplot(rbind(meanCijferHigh, meanCijferLow), beside = T, las = 2)
dev.off()

boxplot(vragenLong$cijfer~vragenLong$vl)
dev.off()

plot(vragenLong$days2 * vragenLong$chronic_recent, vragenLong$Neuroticism, bg = adjustcolor("dodgerblue2", alpha.f = 0.1), pch = 21, col=adjustcolor("dodgerblue2", alpha.f = 0.3))
dev.off()
cor.test(vragenLong$days2 * vragenLong$chronic_recent, vragenLong$Neuroticism)

hist(pheno3$Neuroticism)
dev.off()

boxplot(vragenLong$Neuroticism[!is.na(vragenLong$days)] ~ vragenLong$vl[!is.na(vragenLong$days)])
dev.off()

res2 <- lme(fixed=cijfer~((gender_recent+age_recent+age2_recent+chronic_recent+household_recent+have_childs_at_home_recent + Neuroticism)*days)+ ((gender_recent+age_recent+age2_recent+chronic_recent+household_recent+have_childs_at_home_recent + Neuroticism)*days2), random=~1+Neuroticism|PROJECT_PSEUDO_ID, data= vragenLong,na.action=na.omit)
summary(res2)

anova(res, res2)

res <- lme(fixed=cijfer~ Neuroticism*days, random=~1+gender_recent+age_recent+age2_recent+chronic_recent+household_recent+have_childs_at_home_recent|PROJECT_PSEUDO_ID, data= vragenLong,na.action=na.omit)


















question = qNameMap["hoe waardeert u uw kwaliteit van leven over de afgelopen 7 dagen?",2]
fixedModel <- as.formula(paste(question, "~((gender_recent+age_recent+age2_recent+household_recent+have_childs_at_home_recent+chronic_recent +", paste(colnames(prs)[-1], collapse = " + ") ,"))"))
lmFit <- lm(fixedModel, data = vragenLong[vragenLong$vl == "X2.0",])
summary(lmFit)
