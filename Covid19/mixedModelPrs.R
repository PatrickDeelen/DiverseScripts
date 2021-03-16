#remoter::server(verbose = T, port = 55556, password = "laberkak", sync = T)

remoter::client("localhost", port = 55556, password = "laberkak")

library("nlme")
library("ordinal")
library("GLMMadaptive")

qOverview <- as.matrix(read.delim("/groups/umcg-lifelines/tmp01/projects/ov20_0554/analysis/pgs_correlations/questionair_time_overview_nl.txt", stringsAsFactors = F, row.names = 1))
vls <- colnames(qOverview)[-c(20,21)]

pheno2 <- readRDS("/groups/umcg-lifelines/tmp01/projects/ov20_0554/analysis/risky_behaviour/PRS_correlation/questioniare_subset_participants_with_genome_data/questionaire_df_subset_participants_with_genome_data_01-03-2021.rds")

confounders <- c("gender_recent", "age_recent", "age2_recent", "chronic_recent", "household_recent", "have_childs_at_home_recent")

prs <- read.delim("/groups/umcg-lifelines/tmp01/projects/ov20_0554/analysis/risky_behaviour/PRS_correlation/input_PGS_data_ugli/PGS_combined_ugli_26-02-2021.txt", stringsAsFactors = F)

str(prs)

vragen <- qOverview["hoe waardeert u uw kwaliteit van leven over de afgelopen 7 dagen?",vls]
vragen <- vragen[vragen!=""]

pheno3 <- merge(pheno2, prs, by = "PROJECT_PSEUDO_ID")



pheno3$naCOl <- NA

str(as.list(t(qOverview)))


qNameMap <- data.frame(orginal = row.names(qOverview), new = make.names(row.names(qOverview), unique = T), stringsAsFactors = F)
row.names(qNameMap) <- row.names(qOverview)


qList <- lapply(qNameMap[,1], function(q){
  qs <- qOverview[q,-c(20,21)]
  qs[qs==""]="naCOl"
  qs[!qs %in% colnames(pheno3)] <- "naCOl"
  return(qs)
})
names(qList) <- qNameMap[,2]
qList[[qNameMap["responsdatum covid-vragenlijst",2]]]
qList[[qNameMap["als u moet kiezen, denkt u zelf dat u een coronavirus/covid-19 infectie hebt (gehad)?",2]]]

qList[c(qNameMap["als u moet kiezen, denkt u zelf dat u een coronavirus/covid-19 infectie hebt (gehad)?",2],qNameMap["responsdatum covid-vragenlijst",2])]

#vragen <- qOverview["hoe waardeert u uw kwaliteit van leven over de afgelopen 7 dagen?",vls]
#vragen <- qOverview["als u moet kiezen, denkt u zelf dat u een coronavirus/covid-19 infectie hebt (gehad)?",vls]
#date <- qOverview["responsdatum covid-vragenlijst",vls]
date <- date[vragen!=""]
vragen <- vragen[vragen!=""]


#moet ik nu ook centeren?
#in de intereactie neuro tijd nu indicatie dat effect van prs veranderd met de tijd

#tot 10 aug (vl11)
#omcatten naar dagen vanaf #30 maart

startdate <- as.Date("30/03/2020","%d/%m/%Y")

sum(!vragen %in% colnames(pheno3))




vragenLong <- reshape(pheno3[, c("PROJECT_PSEUDO_ID", vragen, date, confounders, "Neuroticism")], direction = "long", idvar = "PROJECT_PSEUDO_ID", varying = list(date, vragen), v.names = c("date", "cijfer"), times = names(vragen), timevar = "vl")


test <- qList[c(qNameMap["als u moet kiezen, denkt u zelf dat u een coronavirus/covid-19 infectie hebt (gehad)?",2],qNameMap["responsdatum covid-vragenlijst",2])]
names(test)
vragenLong <- reshape(pheno3[, c("PROJECT_PSEUDO_ID", qList[[qNameMap["responsdatum covid-vragenlijst",2]]], qList[[qNameMap["als u moet kiezen, denkt u zelf dat u een coronavirus/covid-19 infectie hebt (gehad)?",2]]], confounders, "Neuroticism")], direction = "long", idvar = "PROJECT_PSEUDO_ID", varying = test, v.names = names(test), times = vls, timevar = "vl")



vragenLong <- reshape(pheno3, direction = "long", idvar = "PROJECT_PSEUDO_ID", varying = qList, v.names = names(qList), times = vls, timevar = "vl")
vragenLong$vl2 <- as.numeric(factor(vragenLong$vl, levels = vls, ordered = T))
vragenLong$vl3 <- factor(vragenLong$vl, levels = vls, ordered = F)
vragenLong$days <- as.numeric(difftime(vragenLong[,qNameMap["responsdatum covid-vragenlijst",2]], startdate ,units="days"))
vragenLong$days2 <- vragenLong$days*vragenLong$days


vragenLong$gender_recent <- factor(vragenLong$gender_recent, levels = 0:1, labels = c("female","male"))
vragenLong$household_recent <- factor(vragenLong$household_recent, levels = 0:1, labels = c("single-person household","multi-person household"))
vragenLong$have_childs_at_home_recent <- factor(vragenLong$have_childs_at_home_recent, levels = 0:1, labels = c("No childeren at home","Childeren at home"))
vragenLong$chronic_recent <- factor(vragenLong$chronic_recent, levels = 0:1, labels = c("Healthy","Chronic disease"))


str(vragenLong$vl3)

str(vragenLong)
head(vragenLong, n =1)

vragenLong[vragenLong$PROJECT_PSEUDO_ID=="00067887-50eb-4d61-9f5b-820de4c18c26",]

hist(vragenLong$days, breaks = 330)
dev.off()

str(vragenLong)


modelLm <- lm()

vragenLong[,"gender_recent"]

question = qNameMap["op hoeveel momenten van de dag eet u iets?",2]
gwas = "Neuroticism"
fixedModel <- as.formula(paste(question, "~((gender_recent+age_recent+age2_recent+household_recent+have_childs_at_home_recent+chronic_recent +", paste(colnames(prs)[-1], collapse = " + ") ,")*days)"))
randomModel <- as.formula("~1|PROJECT_PSEUDO_ID")

res <-  lme(fixed = fixedModel, random=randomModel, data= vragenLong[,c("PROJECT_PSEUDO_ID", question,colnames(prs)[-1],"gender_recent","age_recent","age2_recent","household_recent","have_childs_at_home_recent","chronic_recent", "days" )],na.action=na.omit)
(x <- summary(res))
x$tTable

str(vragenLong[,question])


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
