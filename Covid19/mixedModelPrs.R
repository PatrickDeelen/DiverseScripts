#remoter::server(verbose = T, port = 55556, password = "laberkak", sync = T)

remoter::client("localhost", port = 55556, password = "laberkak")

library("nlme")

qOverview <- as.matrix(read.delim("/groups/umcg-lifelines/tmp01/projects/ov20_0554/analysis/pgs_correlations/questionair_time_overview_nl.txt", stringsAsFactors = F, row.names = 1))
vls <- colnames(qOverview)[-c(20,21)]

pheno2 <- readRDS("/groups/umcg-lifelines/tmp01/projects/ov20_0554/analysis/risky_behaviour/PRS_correlation/questioniare_subset_participants_with_genome_data/questionaire_df_subset_participants_with_genome_data_01-03-2021.rds")

confounders <- c("gender_recent", "age_recent", "age2_recent", "chronic_recent", "household_recent", "have_childs_at_home_recent")

prs <- read.delim("/groups/umcg-lifelines/tmp01/projects/ov20_0554/analysis/risky_behaviour/PRS_correlation/input_PGS_data_ugli/PGS_combined_ugli_26-02-2021.txt", stringsAsFactors = F)

str(prs)

vragen <- qOverview["hoe waardeert u uw kwaliteit van leven over de afgelopen 7 dagen?",vls]
vragen <- vragen[vragen!=""]

pheno3 <- merge(pheno2, prs, by = "PROJECT_PSEUDO_ID")

vragenLong <- reshape(pheno3[, c("PROJECT_PSEUDO_ID", vragen, confounders, "Neuroticism")], direction = "long", idvar = "PROJECT_PSEUDO_ID", varying = list(vragen) , v.names = "cijfer", times = names(vragen), timevar = "vl")
vragenLong$vl <- as.numeric(factor(vragenLong$vl, levels = vls, ordered = T))
str(vragenLong)
head(vragenLong, n =1)

vragenLong[vragenLong$PROJECT_PSEUDO_ID=="00067887-50eb-4d61-9f5b-820de4c18c26",]

(res <- lme(fixed=cijfer~gender_recent+age_recent+chronic_recent+household_recent+have_childs_at_home_recent + (Neuroticism*vl), random=~1+vl|PROJECT_PSEUDO_ID, data= vragenLong,na.action=na.omit))


(res <- lme(fixed=cijfer~(Neuroticism*vl), random=~1+vl|PROJECT_PSEUDO_ID, data= vragenLong,na.action=na.omit))


summary(res)


anova(res)




date <- qOverview["responsdatum covid-vragenlijst",vls]

vragen <- qOverview["hoe waardeert u uw kwaliteit van leven over de afgelopen 7 dagen?",vls]
date <- date[vragen!=""]
vragen <- vragen[vragen!=""]



#moet ik nu ook centeren?
#in de intereactie neuro tijd nu indicatie dat effect van prs veranderd met de tijd

#tot 10 aug (vl11)
#omcatten naar dagen vanaf #30 maart

startdate <- as.Date("30/03/2020","%d/%m/%Y")


vragenLong <- reshape(pheno3[, c("PROJECT_PSEUDO_ID", vragen, date, confounders, "Neuroticism")], direction = "long", idvar = "PROJECT_PSEUDO_ID", varying = list(date, vragen), v.names = c("date", "cijfer"), times = names(vragen), timevar = "vl")
vragenLong$vl2 <- as.numeric(factor(vragenLong$vl, levels = vls, ordered = T))
vragenLong$days <- as.numeric(difftime(vragenLong$date, startdate ,units="days"))
vragenLong$days2 <- vragenLong$days*vragenLong$days


str(vragenLong)
head(vragenLong, n =1)

vragenLong[vragenLong$PROJECT_PSEUDO_ID=="00067887-50eb-4d61-9f5b-820de4c18c26",]

hist(vragenLong$days, breaks = 330)
dev.off()

str(vragenLong)


res <- lme(fixed=cijfer~((gender_recent+age_recent+age2_recent+household_recent+have_childs_at_home_recent + Neuroticism)*days), random=~1+vl|PROJECT_PSEUDO_ID, data= vragenLong,na.action=na.omit)
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
