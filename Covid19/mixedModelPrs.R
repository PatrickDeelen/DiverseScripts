#remoter::server(verbose = T, port = 55556, password = "laberkak", sync = T)

remoter::client("localhost", port = 55556, password = "laberkak")

setwd("/groups/umcg-lifelines/tmp01/projects/ov20_0554/analysis/pgs_correlations/")

library("nlme")
library(heatmap3)
library("ordinal")
library("GLMMadaptive")
library(readr)
library(lme4)
library(meta)

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

sampleQc <- read.delim("/groups/umcg-lifelines/tmp01/projects/ov20_0554/analysis/risky_behaviour/PRS_correlation/inclusionPerVl.txt" , stringsAsFactors = F, row.names = 1)

pheno2 <- pheno2[pheno2$PROJECT_PSEUDO_ID %in% row.names(sampleQc),]
dim(pheno2)

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

selectedTraits <- read.delim("/groups/umcg-lifelines/tmp01/projects/ov20_0554/analysis/pgs_correlations/selectedTraits.txt", header = F, stringsAsFactors = F)[,1]
if(!all(selectedTraits %in% colnames(prs))){
  stop("Not all traits found")
}
prs <- prs[,c("PROJECT_PSEUDO_ID",selectedTraits)]
colnames(prs)[colnames(prs) == "BMI"] <- "BMI_gwas"

if(!all(pheno2$PROJECT_PSEUDO_ID %in% prs$PROJECT_PSEUDO_ID)){
  stop("Not all pheno have genetics")#easly solved but code makes this assumtion
}

pheno2$array <- factor(as.numeric(pheno2$PROJECT_PSEUDO_ID %in% prsGsa$PROJECT_PSEUDO_ID), levels = 0:1, labels = c("Cyto", "Gsa"))

arrayList <- as.list(levels(pheno2$array))
names(arrayList) <- levels(pheno2$array)



pheno3 <- merge(pheno2, prs, by = "PROJECT_PSEUDO_ID")

# Exclude some black lised samples

blackList <- read.delim("/groups/umcg-lifelines/tmp01/projects/ov20_0554/PROJECT_IDS_EARLIER_FILLED_IN.txt", row.names = NULL, sep = "")
sum(pheno3$PROJECT_PSEUDO_ID %in% blackList[,1])
pheno3 <- pheno3[!pheno3$PROJECT_PSEUDO_ID %in% blackList[,1],]


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
  } else if(grepl("responsedate_adu", x)){
    #always include data question
    return(x)
  } else {
    if(sum(is.na(pheno3[,x])) >= (totalPart/4)){
      return (x)
    } else {
      return ("")
    }
  }
})



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



##Create ever covid pos variables

testedPos <- qList[[qNameMap["hebt u een coronavirus/covid-19 infectie (gehad)?",2]]]
everPosQNames <- paste0("everC19Pos",1:ncol(qOverview2))
qList[["everC19Pos"]] <- everPosQNames
qNameMap <- rbind(qNameMap, c("everC19Pos","everC19Pos"))
rownames(qNameMap)[nrow(qNameMap)] <- "everC19Pos"
everPos <- matrix(data = 0, nrow = nrow(pheno3), ncol = ncol(qOverview2) , dimnames = list(row.names = pheno3$PROJECT_PSEUDO_ID, col.names = everPosQNames) )

everPos[!is.na(pheno3[,testedPos[1]]) & pheno3[,testedPos[1]] == 1, 1] <- 1

for(i in 2:ncol(qOverview2)){
  everPos[(!is.na(pheno3[,testedPos[i]]) & pheno3[,testedPos[i]] == 1) | everPos[,i-1] == 1, i] <- 1
}

if(any(colnames(everPos) %in% colnames(pheno3))){
  stop("Duplicate col names")
}

pheno3 <- merge(pheno3, everPos, by.x = "PROJECT_PSEUDO_ID", by.y = 0)

## Reshape to long format and clean some variables

if(any(names(qList) %in% colnames(pheno3))){
  stop("Column name clash after reshape")
}


vragenLong <- reshape(pheno3, direction = "long", idvar = "PROJECT_PSEUDO_ID", varying = qList, v.names = names(qList), times = vls, timevar = "vl")
vragenLong$vl2 <- as.numeric(factor(vragenLong$vl, levels = vls, ordered = T))
vragenLong$vl3 <- factor(vragenLong$vl, levels = vls, ordered = F)
vragenLong$days <- as.numeric(difftime(vragenLong[,qNameMap["responsdatum covid-vragenlijst",2]], startdate ,units="days"))
vragenLong$days2 <- vragenLong$days*vragenLong$days
#vragenLong$days3 <- vragenLong$days*vragenLong$days*vragenLong$days
#vragenLong$days4 <- vragenLong$days*vragenLong$days*vragenLong$days*vragenLong$days
#vragenLong$days5 <- vragenLong$days*vragenLong$days*vragenLong$days*vragenLong$days*vragenLong$days


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


q=qLoop[[2]]
q<-qNameMap["hoe waardeert u uw kwaliteit van leven over de afgelopen 7 dagen?",2]
q<-qNameMap["Positive tested cumsum",2]
q<-qNameMap["hebt u een coronavirus/covid-19 infectie (gehad)?",2]
q<-qNameMap["everC19Pos",2]
q<-qNameMap["BMI",2]
zScoreList <- lapply(qLoop, function(q){
  zScores = tryCatch({
    
    resultsPerArray <- lapply(arrayList, function(array){
      fixedModel <- as.formula(paste(q, "~((gender_recent+age_recent+age2_recent+household_recent+have_childs_at_home_recent+chronic_recent +", paste(colnames(prs)[-1], collapse = " + ") ,")*days + days2 ) "))
      randomModel <- as.formula("~1|PROJECT_PSEUDO_ID")
      res <-  lme(fixed = fixedModel, random=randomModel, data= vragenLong[vragenLong$array == array,c("PROJECT_PSEUDO_ID", q,colnames(prs)[-1],"gender_recent","age_recent","age2_recent","household_recent","have_childs_at_home_recent","chronic_recent", "days", "days2", "array" )],na.action=na.omit)#, control = lmeControl(opt = "optim")
      tTable <- summary(res)$tTable
      
      zscore <- apply(tTable, 1, function(x){
        if(x["t-value"] > 0){
          return(-qnorm(pt(x["t-value"], x["DF"], lower.tail=FALSE, log.p = T),log.p = T))
        } else {
          return(qnorm(pt(x["t-value"], x["DF"], lower.tail=TRUE, log.p = T),log.p = T))
        }
      })
      return(cbind(tTable, zscore))
    })
  }, error = function(e){print(q); return(NULL)})
  return(zScores)
})
cor.test(matrix())
str(zScoreList)
zScoreList2 <- zScoreList[!sapply(zScoreList, is.null)]

# combine into z-score matrix excluding intercept
zscores <- do.call("rbind", zScoreList2)
zscores <- zscores[,colnames(zscores) != "(Intercept)"]

write.table(zscores, file = "zscoreMatrix.txt", sep = "\t", quote = F, col.names = NA)


#below is testingground

plot(resultsPerArray[["Cyto"]][,"zscore"], resultsPerArray[["Gsa"]][,"zscore"])
dev.off()


str(anova(res4,res5))


d <- vragenLong[vragenLong$array == array & !is.na(vragenLong[,q]),]

prsTrait = "Neuroticism"

prsRange <- quantile(d[,prsTrait],probs = seq(0,1,0.1))
dayRange <- range(d[,"days"])


str(res3)

str(tTable)

coef <- tTable[,"Value"]

length(coef)


days = seq(dayRange[1], dayRange[2], 1)
days2 = days * days



plot(days, coef["(Intercept)"] + coef["days"] * days + coef["days2"] * days2 + prsRange[10] * days * coef[paste0(prsTrait,":days")] + coef[prsTrait] * prsRange[10], col = "red", type = "l", ylab = q, xlab = "Dagen sinds 30 maart 2020", main = paste0("Red is high ", prsTrait, " PRS\nBlue is low PRS"))
points(days, coef["(Intercept)"] + coef["days"] * days + coef["days2"] * days2 + prsRange[1] * days * coef[paste0(prsTrait,":days")] + coef[prsTrait] * prsRange[1], col = "blue", type = "l")
dev.off()

plot(days, coef["(Intercept)"] + coef["days"] * days + coef["days2"] * days2 + prsRange[10] * days * coef["Neuroticism:days"] + coef["Neuroticism"] * prsRange[5], col = "red", type = "l", ylim = c(2,7), ylab = "Quality of life prediciton by model", xlab = "Dagen sinds 30 maart 2020", main = "Red is high neurotism PRS and blue is low PRS")
points(days, coef["(Intercept)"] + coef["days"] * days + coef["days2"] * days2 + prsRange[1] * days * coef["Neuroticism:days"] + coef["Neuroticism"] * prsRange[5], col = "blue", type = "l")
dev.off()

plot(days, coef["days"] * days + coef["days2"] * days2 + 1 * days * coef["chronic_recentChronic disease:days"] + coef["chronic_recentChronic disease"] * 0, col = "red", type = "l", ylim = c(-4,0.25), ylab = "Quality of life prediciton by model", xlab = "Dagen sinds 30 maart 2020", main = "Red is high neurotism PRS and blue is low PRS")
points(days,  coef["days"] * days + coef["days2"] * days2 + 0 * days * coef["chronic_recentChronic disease:days"] + coef["chronic_recentChronic disease"] * 0, col = "blue", type = "l")
dev.off()


plot(fitted.values(res), d[,q])
dev.off()

plot(d[,"days"], d[,q])
abline(res)
dev.off()

sum(is.na(vragenLong[!is.na(vragenLong[,q]),c(colnames(prs)[-1])]))



cor.test(fitted(res), vragenLong[!is.na(vragenLong[,q]),q])

write.table(summary(res)$tTable, file = "tmp3.txt", sep = "\t", quote = F, col.names = NA)

dim(zscores)

pdf("mixedModelZscores.pdf", width = 50, height = 50)
heatmap3(zscores, balanceColor = T, scale = "none", margins = c(20,20))
dev.off()

x <- zScoreList[[1]]
x[order(abs(x))]

q <- qLoop[[2]]

zScores[order(abs(zScores))]

###### Below is old testing stuff

tTable

lmFit <- glm(fixedModel,data=vragenLong[vragenLong$array == array,c("PROJECT_PSEUDO_ID", q,colnames(prs)[-1],"gender_recent","age_recent","age2_recent","household_recent","have_childs_at_home_recent","chronic_recent", "days", "days2", "array" )])

#### @Robert fit for normal binom below


glmBinomFit <- glm(fixedModel ,family=binomial(link='logit'),data=vragenLong[vragenLong$array == array,c("PROJECT_PSEUDO_ID", q,colnames(prs)[-1],"gender_recent","age_recent","age2_recent","household_recent","have_childs_at_home_recent","chronic_recent", "days", "days2", "array" )])



prsTrait = "COVID.19.susceptibility"
prsRange <- quantile(prs[,prsTrait],probs = seq(0,1,0.1))

dummy <- vragenLong[1:307,c(q,colnames(prs)[-1],"gender_recent","age_recent","age2_recent","household_recent","have_childs_at_home_recent","chronic_recent", "days", "days2")]
dummy$days <- 1:307
dummy$days2 <- dummy$days * dummy$days
for(prsCol in colnames(prs)[-1]){
  dummy[,prsCol] <- mean(prs[,prsCol])
}
dummy[,"age_recent"] <- mean(pheno3$age_recent)
dummy[,"age2_recent"] <- dummy[,"age_recent"] * dummy[,"age_recent"]
dummy[,"household_recent"] <- levels(dummy[,"household_recent"])[1]
dummy[,"have_childs_at_home_recent"] <- levels(dummy[,"have_childs_at_home_recent"])[1]
dummy[,"gender_recent"] <- levels(dummy[,"gender_recent"])[1]
dummy[,"chronic_recent"] <- levels(dummy[,"chronic_recent"])[1]

dummy[,prsTrait] <- prsRange[1]
predictLow <- predict(glmBinomFit, dummy, type = "response")
dummy[,prsTrait] <- mean(prs[,prsTrait])
predictMedium <- predict(glmBinomFit, dummy, type = "response")
dummy[,prsTrait] <- prsRange[10]
predictHigh <- predict(glmBinomFit, dummy, type = "response")


plot.new()
plot.window(xlim = range(dummy$days), ylim = range(predictLow, predictHigh))
axis(side = 1)
axis(side = 2)
title(xlab = "days", ylab = "C19 pos")
points(predictLow, col = "blue", type = "l")
points(predictMedium, col = "green", type = "l")
points(predictHigh, col = "red", type = "l")
dev.off()



######### @Robert mixed model below





model <- as.formula(paste(q, "~((gender_recent+age_recent+age2_recent+household_recent+have_childs_at_home_recent+chronic_recent +", paste(colnames(prs)[-1], collapse = " + ") ,")*days + days2 ) + (1|PROJECT_PSEUDO_ID) "))
d <- vragenLong[vragenLong$array == array & !is.na(vragenLong[,q]) & !is.na(vragenLong[,"days"]),c("PROJECT_PSEUDO_ID", q,colnames(prs),"gender_recent","age_recent","age2_recent","household_recent","have_childs_at_home_recent","chronic_recent", "days", "days2", "array" )]
#d[,q] <- d[,q]-1
d[,q] <- as.factor(d[,q])
table(d[,q] )

apply(d, 2, range)

dMean <- mean(d$days)
dSd <- sd(d$days)

(range(d$days)-dMean)/dSd

d2 <- lapply(d, function(x){
  if(is.numeric(x)){
    return(scale(x))
  } else{
    return(x)
  }
})

d2 <- as.data.frame(d2, stringsAsFactors = F)
str(d2)

glmMerFit <- glmer(model, data = d2, family = binomial, nAGQ=0 )#glmerControl(optimizer = "nloptwrap2") #control = glmerControl(optimizer = "bobyqa")  
summary(glmMerFit)$coefficients



prsTrait = "COVID.19.susceptibility"
prsRange <- quantile(d2[,prsTrait],probs = seq(0,1,0.1))


dummy <- d2[1:307,c("PROJECT_PSEUDO_ID", q,colnames(prs)[-1],"gender_recent","age_recent","age2_recent","household_recent","have_childs_at_home_recent","chronic_recent", "days", "days2")]
dummy$PROJECT_PSEUDO_ID = "A"
dummy$days <- (1:307-dMean)/dSd
dummy$days2 <- dummy$days * dummy$days
dummy[,colnames(prs)] <- 0
dummy[,"age_recent"] <- 0
dummy[,"age2_recent"] <- 0
dummy[,"household_recent"] <- levels(dummy[,"household_recent"])[1]
dummy[,"have_childs_at_home_recent"] <- levels(dummy[,"have_childs_at_home_recent"])[1]
dummy[,"gender_recent"] <- levels(dummy[,"gender_recent"])[1]
dummy[,"chronic_recent"] <- levels(dummy[,"chronic_recent"])[1]

dummy[,prsTrait] <- prsRange[1]
predictLow <- predict(glmMerFit, dummy, type = "response", re.form = NA)
dummy[,prsTrait] <- mean(d2[,prsTrait])
predictMedium <- predict(glmMerFit, dummy, type = "response", re.form = NA)
dummy[,prsTrait] <- prsRange[10]
predictHigh <- predict(glmMerFit, dummy, type = "response", re.form = NA)


plot.new()
plot.window(xlim = c(1,307), ylim = range(predictLow, predictHigh))
axis(side = 1)
axis(side = 2)
title(xlab = "days", ylab = "C19 pos")
points(predictLow, col = "blue", type = "l")
points(predictMedium, col = "green", type = "l")
points(predictHigh, col = "red", type = "l")
dev.off()







#### Other









summary(glmBinomFit)$coefficients  


plot(summary(lmFit)$coefficients[,"t value"], summary(glmBinomFit)$coefficients[,"z value"] )
dev.off()



write.table(cbind(summary(lmFit)$coefficients, summary(glmBinomFit)$coefficients , tTable, summary(glmMerFit)$coefficients), file = "test.txt", sep ="\t",quote = F)





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





fixedModel <- as.formula(paste(q, "~((gender_recent+age_recent+age2_recent+household_recent+have_childs_at_home_recent+chronic_recent +", paste(colnames(prs)[-1], collapse = " + ") ,")*days + days2) "))
randomModel <- as.formula("~1|PROJECT_PSEUDO_ID")
d <- vragenLong[vragenLong$array == array & !is.na(vragenLong[,q]),c("PROJECT_PSEUDO_ID", q,colnames(prs)[-1],"gender_recent","age_recent","age2_recent","household_recent","have_childs_at_home_recent","chronic_recent", "days", "days2", "array" )]
#d[,q] <- d[,q]-1
oridnalFit <-  mixed_model(fixed = fixedModel, random=randomModel, data=d ,na.action=na.omit, family = poisson())
table(d[,q])
summary(oridnalFit)

model <- as.formula(paste(q, "~((gender_recent+age_recent+age2_recent+household_recent+have_childs_at_home_recent+chronic_recent +", paste(colnames(prs)[-1], collapse = " + ") ,")*days + days2 ) + (1|PROJECT_PSEUDO_ID) "))
model <- as.formula(paste(q, "~((gender_recent+age_recent+age2_recent+household_recent+have_childs_at_home_recent+chronic_recent + COVID.19.susceptibility)*days + days2) + (1|PROJECT_PSEUDO_ID) "))
d <- vragenLong[vragenLong$array == array & !is.na(vragenLong[,q]) & !is.na(vragenLong[,"days"]),c("PROJECT_PSEUDO_ID", q,colnames(prs),"gender_recent","age_recent","age2_recent","household_recent","have_childs_at_home_recent","chronic_recent", "days", "days2", "array" )]
#d[,q] <- d[,q]-1
d[,q] <- as.factor(d[,q])
table(d[,q] )

apply(d, 2, range)

str(d)

d2 <- lapply(d, function(x){
  if(is.numeric(x)){
    return(scale(x))
  } else{
    return(x)
  }
})

d2 <- as.data.frame(d2, stringsAsFactors = F)
str(d2)

d$days <- scale(d$days)
d$days2 <- d$days * d$days

d$age_recent <- scale(d$age_recent)
d$age2_recent <- d$age_recent * d$age_recent


library(nloptr)
defaultControl <- list(algorithm="NLOPT_LN_BOBYQA",xtol_rel=1e-6,maxeval=1e5)
nloptwrap2 <- function(fn,par,lower,upper,control=list(),...) {
  for (n in names(defaultControl)) 
    if (is.null(control[[n]])) control[[n]] <- defaultControl[[n]]
    res <- nloptr(x0=par,eval_f=fn,lb=lower,ub=upper,opts=control,...)
    with(res,list(par=solution,
                  fval=objective,
                  feval=iterations,
                  conv=if (status>0) 0 else status,
                  message=message))
}

glmMerFit <- glmer(model, data = d2, family = binomial, nAGQ=0 )#glmerControl(optimizer = "nloptwrap2") #control = glmerControl(optimizer = "bobyqa")  
summary(glmMerFit)$coefficients

lapply(d2, mean, na.rm = T)


#vragenLong[,qNameMap["hoe waardeert u uw kwaliteit van leven over de afgelopen 7 dagen?",2]] <- factor(vragenLong[,qNameMap["hoe waardeert u uw kwaliteit van leven over de afgelopen 7 dagen?",2]], levels = 1:10,  ordered = T)
model <- as.formula(paste(q, "~((gender_recent+age_recent+age2_recent+household_recent+have_childs_at_home_recent+chronic_recent +", paste(colnames(prs)[-1], collapse = " + ") ,")*days + days2) + (1|PROJECT_PSEUDO_ID) "))
d <- vragenLong[vragenLong$array == array & !is.na(vragenLong[,q]),c("PROJECT_PSEUDO_ID", q,colnames(prs)[-1],"gender_recent","age_recent","age2_recent","household_recent","have_childs_at_home_recent","chronic_recent", "days", "days2", "array" )]
d[,q] <- as.factor(d[,q])
table(d[,q] )
ordinalFit <- clmm(model, data = d,na.action=na.omit, link = "logit")

wares <- lme(fixed=cijfer~((COVID.19.susceptibility)*days), random=~1|PROJECT_PSEUDO_ID, data= vragenLong,na.action=na.omit)
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

resultsPerArray[[1]][1,]
resultsPerArray[[2]][1,]

metagen(c(5.864397e+00, 6.849405e+00), c(7.555978e-01, 4.037785e-01))$zval.fixed

x <- (resultsPerArray[[1]][1,1]/(resultsPerArray[[1]][1,2]*resultsPerArray[[1]][1,2])) + (resultsPerArray[[2]][1,1]/(resultsPerArray[[2]][1,2]*resultsPerArray[[2]][1,2]))
y <- (1/(resultsPerArray[[1]][1,2]*resultsPerArray[[1]][1,2])) + (1/(resultsPerArray[[2]][1,2]*resultsPerArray[[2]][1,2]))


a <- x/y


b <- sqrt(1/y)

a/b



table(pheno3$array)

((resultsPerArray[[1]][1,6] * 5287) + (resultsPerArray[[2]][1,6] * 12532)) / sqrt((5287*5287)+(12532*12532))


cat(colnames(prs), sep = "\n")
