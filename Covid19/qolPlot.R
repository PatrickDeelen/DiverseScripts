
range(vragenLong[,"days"], na.rm=T)
range(vragenLong[,qNameMap["responsdatum covid-vragenlijst",2]], na.rm=T)


#library("viridis")
#library(scales)

#show_col(viridis_pal(option= "C")(12))
#dev.off()

qName <- "hoe waardeert u uw kwaliteit van leven over de afgelopen 14 dagen? (include 7 days)"
q<-qNameMap[qName,2]
metaRes <- resultList[[qName]]
effect <- "Life.satisfaction:days"
metaRes<- as.matrix(metaRes)

effectName <- sub(":days", "", effect)


colHigh = "#6300A7"
colMedium = "#D5546E"
colLow = "#FCD225"
colAxis = "grey70"
colMean = "#808080"

prsRange <- quantile(prs[,effectName],probs = seq(0,1,0.1))

prsLabel = prsLabels[effectName]

qInfo <- selectedQ[q,]

qLable <- qInfo[,"English.label"]

daysSeq <- qInfo[,"firstDay"]:qInfo[,"lastDay"]


pdf("qolFigure.pdf", width = 3.5, height = 7, useDingbats = F)
#rpng()
layout(matrix(c(1,2,3,4), nrow=4, byrow = T), heights = c(1,1,1,1))

par(mar = c(3,5,1,1), xpd = NA, las = 1, cex.lab = 0.6, cex.axis = 0.6)

library(tidyverse)
qualityOfLifeQuestion <- "hoe.waardeert.u.uw.kwaliteit.van.leven.over.de.afgelopen.14.dagen...include.7.days."
qualityOfLifeTable <- vragenLong %>% 
  select(vl2, !!!qualityOfLifeQuestion, responsdatum.covid.vragenlijst) %>%
  filter(!is.na(responsdatum.covid.vragenlijst) & !is.na(get(qualityOfLifeQuestion))) %>%
  group_by(vl2) %>%
  summarise(avgResponseDate = as.Date(as.character(mean(responsdatum.covid.vragenlijst))), avgQoL = mean(get(qualityOfLifeQuestion)))



days <- as.numeric(difftime(qualityOfLifeTable[,"avgResponseDate", drop = T], startdate ,units="days"))
qol <- qualityOfLifeTable[,"avgQoL", drop = T]


endDate <- startdate + 307

axisAt <- c(startdate,startdate+100,startdate+200, endDate)


plot.new()
plot.window(xlim = c(startdate,endDate), ylim = c(6.8,7.8))
axis(side = 1, at = axisAt, labels = format(axisAt, "%d-%b-%Y"), col = colAxis, col.axis = colMean)
axis(side = 2, col = colAxis, col.axis = colMean)
title(ylab = paste0("Mean ",qLable), col.lab = colMean)
points(startdate + as.numeric(days) , qol, col = colMean, type = "p", lwd = 2, pch = 16)


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

legend("bottomleft", fill = rev(c(colLow, colMedium, colHigh)), legend = paste0(rev(c("Lowest 10% PGS for ", "Median PGS for ", "Highest 10% PGS for ")), prsLabel), bty = "n", cex = 0.6, pt.cex = 1, text.col = colMean, border= NA)
  
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

betas <- replication["Life satisfaction / Precieved quality of life include 7 days",!is.na(replication["Life satisfaction / Precieved quality of life include 7 days",])]
days <- replication[1,!is.na(replication["Life satisfaction / Precieved quality of life include 7 days",])]


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


