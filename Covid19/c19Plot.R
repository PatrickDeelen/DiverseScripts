workdir <- "/groups/umcg-lifelines/tmp01/projects/ov20_0554/analysis/pgs_correlations/"
phenoPath <- "/groups/umcg-lifelines/tmp01/projects/ov20_0554/analysis/risky_behaviour/PRS_correlation/combined_questionnaires_v22_14-05-2021_genome_fitered_participants_filtered_correct_recoded/questionnaire_subset_participants_filtered_recoded_answers_longitudinal_filtered_15-05-2021"
prsLabelFile <- "prsLables.txt"
preparedDataFile <- "longitudinal.RData"
intermediatesdir <-  "/groups/umcg-lifelines/tmp01/projects/ov20_0554/analysis/pgs_correlations/longiIntermediates/"

setwd(workdir)
load(preparedDataFile)

library(survival)
#library("viridis")
#library(scales)

#show_col(viridis_pal(option= "C")(12))
#dev.off()

qName <- "Positive tested cumsum"
q<-qNameMap[qName,2]




colHigh = "#6300A7"
colMedium = "#D5546E"
colLow = "#FCD225"
colAxis = "grey70"
colMean = "#808080"


startdate <- as.Date("30/03/2020","%d/%m/%Y")
startday <- 0
endday <- startday + 307
axisAt <- c(startday,startday+100,startday+200, endday)


qInfo <- selectedQ[q,]

qLable <- qInfo[,"label_en"]




covid19Events <- lapply(pheno3$PROJECT_PSEUDO_ID, function(id){
  idVragen <- vragenLong[vragenLong[,"PROJECT_PSEUDO_ID"] == id,c(qNameMap["hebt u een coronavirus/covid-19 infectie (gehad)?",2],"vl", qNameMap["responsdatum covid-vragenlijst",2],usedPrs,"gender_recent","age_recent","age2_recent","household_recent","have_childs_at_home_recent","chronic_recent", "array", "days", "days2")]
  
  
  firstVlPos <- which(idVragen[,qNameMap["hebt u een coronavirus/covid-19 infectie (gehad)?",2]] == 1)[1]
  
  
  if(!is.na(firstVlPos)){
    firstPosDate <- idVragen[firstVlPos,qNameMap["responsdatum covid-vragenlijst",2]]
    
    eventInfo <- idVragen[firstVlPos,]
    eventInfo$eventDate <- firstPosDate
    eventInfo$eventDays <- as.numeric(difftime(firstPosDate, startdate ,units="days"))
    
    return(eventInfo)
  } else{
    #now test for censoring 
    
    censorDate <- 0
    eventInfo <- idVragen[max(which(!is.na(idVragen[,"responsdatum.covid.vragenlijst"]))),]
    lastCompletedVl = eventInfo[,"vl"]
    
    if(lastCompletedVl == "X17.0"){
      #completed
      censorDate <- endDate
    } else {
      #censoring
      censorDate <- vlDates[which(rownames(vlDates) == lastCompletedVl)+1,]
    }
    
    eventInfo$eventDate <- censorDate
    eventInfo$eventDays <- as.numeric(difftime(censorDate, startdate ,units="days"))
    eventInfo[,qNameMap["hebt u een coronavirus/covid-19 infectie (gehad)?",2]] <- 0
    
    
    return(eventInfo)
    
    
  }
  
  
  
})


covid19Events2 <- do.call("rbind", covid19Events)

prsRangeKm <- quantile(prs[,effectName],probs = c(0,0.1,0.4,0.6,0.9,1))
covid19Events2$prsQuantile <- cut(covid19Events2[,effectName],breaks = prsRangeKm, include.lowest = T)
length(table(covid19Events2$prsQuantile))


table(covid19Events2[,"eventDays"], useNA="always")
table(covid19Events2[,q], useNA="always")

survModel <- as.formula(paste(" Surv(eventDays, ",q,") ~prsQuantile"))


km_fit <- survfit(survModel, data=covid19Events2)
summary(km_fit, times = c(1,180,max(covid19Events2$eventDays)))          



pdf("c19Figure.pdf", width = 3.5, height = 7, useDingbats = F)
#rpng()
layout(matrix(c(1,2,3,4), nrow=4, byrow = T), heights = c(1,1,1,1))

par(mar = c(3,5,1,1), xpd = NA, las = 1, cex.lab = 0.6, cex.axis = 0.6)


plot.new()
plot.window(xlim =c(0,307), ylim = c(0.8,1))
axis(side = 1, at = axisAt-18350, labels = format(axisAt, "%d-%b-%Y"), col = colAxis, col.axis = colMean)
axis(side = 2, col = colAxis, col.axis = colMean)
title(ylab = paste0("Mean ",qLable), col.lab = colMean)
lines(km_fit, xlab="Days", main = 'C19 incidence', col = c(colLow, NA,colMedium,NA, colHigh), lwd=2, mark.time=FALSE, ylim = c(0.8, 1)) 
#legend("bottomleft", fill = c(colLow, colMedium, colHigh), legend = paste0(c("Lowest 10% PGS of ", "Median PGS of ", "Highest 10% PGS of "), prsLabel), bty = "n")





dummy <- vragenLong[daysSeq,c(q,colnames(prs)[-1],"gender_recent","age_recent","age2_recent","household_recent","have_childs_at_home_recent","chronic_recent", "days", "days2")]
dummy$days <- daysSeq
dummy$days2 <- dummy$days * dummy$days
for(prsCol in colnames(prs)[-1]){
  dummy[,prsCol] <- mean(prs[,prsCol])
  #dummy[,prsCol] <- 0
}
dummy[,"age_recent"] <- mean(pheno3$age_recent)
#dummy[,"age_recent"] <- 0
dummy[,"age2_recent"] <- dummy[,"age_recent"] * dummy[,"age_recent"]
dummy[,"household_recent"] <- factor(levels(dummy[,"household_recent"])[1], levels = levels(dummy[,"household_recent"]))
dummy[,"have_childs_at_home_recent"] <- factor(levels(dummy[,"have_childs_at_home_recent"])[1], levels = levels(dummy[,"have_childs_at_home_recent"]))
dummy[,"gender_recent"] <- factor(levels(dummy[,"gender_recent"])[1], levels = levels(dummy[,"gender_recent"]))
dummy[,"chronic_recent"] <- factor(levels(dummy[,"chronic_recent"])[1], levels = levels(dummy[,"chronic_recent"]))

fam <- NA

if(qInfo["Type"] == "gaussian"){
  fam <- gaussian()
} else if (qInfo["Type"] == "binomial"){
  fam <- binomial(link = "logit")
} else {
  stop("error")
}

coef=metaRes[,"y"]


dummy[,effectName] <- prsRange[10]
highPrs <- predict_meta(df = dummy, coefficients = coef, family = fam)#, family = binomial(link = "logit")
dummy[,effectName] <- prsRange[6]
medianPrs <- predict_meta(df = dummy, coefficients = coef, family = fam)#, family = binomial(link = "logit")
dummy[,effectName] <- prsRange[2]
lowPrs <- predict_meta(df = dummy, coefficients = coef, family = fam)#, family = binomial(link = "logit")



plot.new()
plot.window(xlim = c(startdate,endDate), ylim = range(lowPrs, medianPrs, highPrs))
axis(side = 1, at = axisAt, labels = format(axisAt, "%d-%b-%Y"), col = colAxis, col.axis = colMean)
axis(side = 2, col = colAxis, col.axis = colMean)
title(ylab = paste0("Predicted of ", qLable, "\n by full model"), col.lab = colMean)
points(startdate +daysSeq, lowPrs, col = colLow, type = "l", lwd = 2)
points(startdate +daysSeq, medianPrs, col = colMedium, type = "l", lwd = 2)
points(startdate +daysSeq, highPrs, col = colHigh, type = "l", lwd = 2)

legend("topleft", fill = rev(c(colLow, colMedium, colHigh)), legend = paste0(rev(c("Lowest 10% PGS of ", "Median PGS of ", "Highest 10% PGS for ")), prsLabel), bty = "n", cex = 0.6, pt.cex = 1, text.col = colMean, border= NA)
  
coef[!grepl(effectName, names(coef))] <- 0

dummy[,effectName] <- prsRange[10]
highPrs <- predict_meta(df = dummy, coefficients = coef, family = fam)#, family = binomial(link = "logit")
dummy[,effectName] <- prsRange[6]
medianPrs <- predict_meta(df = dummy, coefficients = coef, family = fam)#, family = binomial(link = "logit")
dummy[,effectName] <- prsRange[2]
lowPrs <- predict_meta(df = dummy, coefficients = coef, family = fam)#, family = binomial(link = "logit")


plot.new()
plot.window(xlim = c(startdate,endDate), ylim = range(lowPrs, medianPrs, highPrs))
axis(side = 1, at = axisAt, labels = format(axisAt, "%d-%b-%Y"), col = colAxis, col.axis = colMean)
axis(side = 2, col = colAxis, col.axis = colMean)
title(ylab = paste0("Contribution of ", prsLabel, " PGS\non ", qLable, " prediction"), col.lab = colMean)
points(startdate +daysSeq, lowPrs, col = colLow, type = "l", lwd = 2)
points(startdate +daysSeq, medianPrs, col = colMedium, type = "l", lwd = 2)
points(startdate +daysSeq, highPrs, col = colHigh, type = "l", lwd = 2)



replication <- read.delim("replication_data_points.txt",row.names =1 , stringsAsFactors = F, header = F)
str(replication)

betas <- replication["COVID-19 susceptibility / Ever positive SARS-CoV-2 PCR test",!is.na(replication["COVID-19 susceptibility / Ever positive SARS-CoV-2 PCR test",])]
days <- replication[1,!is.na(replication["COVID-19 susceptibility / Ever positive SARS-CoV-2 PCR test",])]


plot.new()
plot.window(xlim = c(startdate,endDate), ylim = range(betas))
axis(side = 1, at = axisAt, labels = format(axisAt, "%d-%b-%Y"), col = colAxis, col.axis = colMean)
axis(side = 2, col = colAxis, col.axis = colMean)
title(ylab = paste0("Regression coeffcients for ", prsLabel, " PGS \n per questionnaire"), col.lab = colMean)

points(startdate + as.numeric(days), as.numeric(betas), col = colMean, type = "p", lwd = 2, pch = 16)


days2 <- as.numeric(startdate) + as.numeric(days)
abline(lm(as.numeric(betas) ~ days2), xpd = F, col = colMean)

cor.test(as.numeric(betas), as.numeric(days))


dev.off()


