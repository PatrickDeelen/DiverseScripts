#remoter::server(verbose = T, port = 55556, password = "laberkak", sync = T)

remoter::client("localhost", port = 55556, password = "laberkak")

setwd("/groups/umcg-lifelines/tmp01/projects/ov20_0554/analysis/pgs_correlations/")



library("nlme")
library(heatmap3)
#library("ordinal")
#library("GLMMadaptive")
library(readr)
library(lme4)
#library(meta)
library(survival)




#Functions:
spread_factor_columns <- function(df) {
  variables <- colnames(df)
  # For all columns, check if it is a factor,
  # if so, spread levels across columns as 0/1.
  # else, copy the column as is.
  out <- do.call("cbind", lapply(variables, function(variable) {
    col <- df[,variable]
    if (is.factor(col)) {
      # For all levels, add a numeric column (0/1) indicating the
      # factor value for the row.
      uniqueLevels <- levels(col)
      colDf <- as.data.frame(lapply(uniqueLevels, function(factorLevel) {
        levelBool <- as.numeric(col == factorLevel)
      }))
      colnames(colDf) <- paste0(variable, uniqueLevels)
      return(colDf)
    } else {
      colDf <- data.frame(col)
      colnames(colDf) <- variable
      return(colDf)
    }
  }))
  return(out)
}
# Function to predict outcome values given a new dataset, 
# named coefficient list (from a meta-analysis),
# the model on which analysis is based
predict_meta <- function(df, coefficients, family) {
  # Spread all factor levels across columns for all columns that are factors.
  dfSpread <- spread_factor_columns(df)
  # Calculate the predicted values per term
  predictionPerTerm <- do.call("cbind", sapply(
    names(coefficients), 
    function(coefficientLabel) {
      # Coefficient
      coefficient <- coefficients[coefficientLabel]
      # Parse coefficient label,
      # Find out which columns should be multiplied
      variablesToMultiply <- strsplit(coefficientLabel, ":", fixed = T)[[1]]
      if (all(variablesToMultiply == "(Intercept)")) {
        return(data.frame(coefficientLabel = rep(coefficient, nrow(dfSpread))))
      } else if (length(variablesToMultiply) == 2) {
        return(Reduce(`*`, dfSpread[variablesToMultiply]) * coefficient)
      } else if (length(variablesToMultiply == 1)) {
        return(dfSpread[variablesToMultiply] * coefficient)
      } else {
        stop("Something went wrong, or we are going to multiply more than two terms. Check if this works first!")
      }
    }))
  # Calculate the predicted values on a regular linear scale
  predicted <- rowSums(predictionPerTerm)
  # Return values according to the linear inverse of the link function
  return(family$linkinv(predicted))
}


inverseVarianceMeta <- function(resultsPerArray, seCol, valueCol){
  x <- as.data.frame(resultsPerArray[[1]][,FALSE])
  x$sumYDivSe2 <- 0
  x$sum1DivSe2 <- 0
  
  for(array in names(resultsPerArray)){
    se2 <- resultsPerArray[[array]][,seCol] * resultsPerArray[[array]][,seCol]
    x$sumYDivSe2 <- x$sumYDivSe2 + (resultsPerArray[[array]][,valueCol]/ se2)
    x$sum1DivSe2 <- x$sum1DivSe2 + (1/se2)
  }
  
  metaRes <- as.data.frame(resultsPerArray[[1]][,FALSE])
  metaRes$y <- x$sumYDivSe2/x$sum1DivSe2
  metaRes$se <- sqrt(1/x$sum1DivSe2)
  metaRes$z <- metaRes$y/metaRes$se 
  metaRes$p <- 2*pnorm(-abs(metaRes$z))
  return(metaRes)
}


validationSamples <- read.delim("/groups/umcg-lifelines/tmp01/projects/ov20_0554/analysis/risky_behaviour/PRS_correlation/validationSamples.txt", header = F)[,1]

vragenLongTest <- vragenLong[ (vragenLong$PROJECT_PSEUDO_ID %in% validationSamples),]
table(vragenLongTest$vl)




table(selectedQ[,"Type"])

q=qLoop[[4]]
q<-qNameMap["hoe waardeert u uw kwaliteit van leven over de afgelopen 14 dagen?",2]
q<-qNameMap["hoe waardeert u uw kwaliteit van leven over de afgelopen 14 dagen? (include 7 days)",2]
q<-qNameMap["everC19Pos",2]
q<-qNameMap["hoeveel zorgen maakte u zich de afgelopen 14 dagen over de corona-crisis?",2]
q<-qNameMap["Positive tested cumsum",2]
q<-qNameMap["BMI",2]
q<-qNameMap["hoe waardeert u uw kwaliteit van leven over de afgelopen 14 dagen?",2]
q<-qNameMap["ik heb vertrouwen in de aanpak van de corona-crisis door de nederlandse regering",2]
q<-"Mini.combined..depressief"
resultList <- lapply(qLoop, function(q){
  #zScores = tryCatch({
    
    print(q)
    
    qInfo <- selectedQ[q,]
    usedPrs <- colnames(prs)[-1]
   # usedPrs <- usedPrs[!usedPrs %in% "Cigarettes.per.day"]
    #usedPrs <- c("BMI_gwas", "Life.satisfaction", "Neuroticism")
    #usedPrs <- "Life.satisfaction"
    #usedPrs <- "Neuroticism"
   #usedPrs <- "BMI_gwas"
    #usedPrs <- "Cigarettes.per.day"
    #usedPrs <- "COVID.19.susceptibility"
    #usedPrs <- "Anxiety.tension"
    #usedPrs <- "COVID.19.severity"
    #usedPrs <- "Worry.vulnerability"
    #
    
    fixedString <- paste(q, "~((gender_recent+age_recent+age2_recent+household_recent+have_childs_at_home_recent+chronic_recent +", paste0(usedPrs, collapse = " + ") ,")*days + I(days^2) ) ")
    randomString <- "1|PROJECT_PSEUDO_ID"
    fixedModel <- as.formula(fixedString)
    randomModel <- as.formula(paste0("~",randomString))
    fullModel <- as.formula(paste0(fixedString, "+ (", randomString, ")"))
    
    resultsPerArray <- lapply(arrayList, function(array){
      
      #& vragenLongValidation$gender_recent == "female" 
      d <- vragenLongTest[!is.na(vragenLongTest[,q])  & vragenLongTest$array == array,c("PROJECT_PSEUDO_ID", q,usedPrs,"gender_recent","age_recent","age2_recent","household_recent","have_childs_at_home_recent","chronic_recent", "days", "days2", "days3", "vl")]
      table(d[,q])
      coef <- 0
      
      if(qInfo["Type"] == "gaussian" & qInfo["Mixed"]){
        print("test1")
        res <-  lme(fixed = fixedModel, random=randomModel, data=d,na.action=na.omit, control = lmeControl(opt = "optim"))#control = lmeControl(opt = "optim")
        coef <- summary(res)$tTable
      } else if (qInfo["Type"] == "gaussian" & !qInfo["Mixed"]) {
        print("test2")
        stop("Not implement")
      } else if (qInfo["Type"] == "binomial" & qInfo["Mixed"]) {
        print("test3")
        if(max(d[,q])==2){
          d[,q] <- d[,q] -1
        }
        if(sum(range(d[,q])==0:1)!=2){
          stop("not binomal")
        }
        glmMerFit <- glmer(fullModel, data = d, family = binomial, nAGQ=0 )
        coef <- summary(glmMerFit)$coefficients
        colnames(coef)[1:2]<-c("Value", "Std.Error")
      } else if (qInfo["Type"] == "binomial" & !qInfo["Mixed"]) {
        print("test4")
        d[,q] <- as.factor(d[,q])
        glmBinomFit <- glm(fixedModel ,family=binomial(link='logit'),data=d)
        coef <- summary(glmBinomFit)$coefficients
        colnames(coef)[1:2]<-c("Value", "Std.Error")
      }

      return(coef)
    })
    
    #resultsPerArray[["Gsa"]]
    #resultsPerArray[["Cyto"]]
    
      metaRes <- inverseVarianceMeta(resultsPerArray, "Std.Error", "Value")
      metaRes
      return(as.matrix(metaRes))
  #}, error = function(e){print("error: ", q); return(NULL)})
 # return(zScores)
})

glm


pheno3$test <- is.na(pheno3[,"covt17_responsedate_adu_q_1"])


cor.test(pheno3$age_recent, pheno3$covt01_bmi)

summary(glm(test~covt04_bmi,family=binomial(link='logit'),data=pheno3))

summary(glm(test~covt04_bmi+age_recent,family=binomial(link='logit'),data=pheno3))


summary(glm(test~covt04_bmi+gender_recent+age_recent+age2_recent+household_recent+have_childs_at_home_recent+chronic_recent ,family=binomial(link='logit'),data=pheno3))


usedPrs <- colnames(prs)[-1]




summary(glm(fixedModel <- as.formula(paste("test~covt04_bmi+covt11_qualityoflife_adu_q_2+gender_recent+age_recent+age2_recent+household_recent+have_childs_at_home_recent+chronic_recent +", paste0(usedPrs, collapse = " + "))),family=binomial(link='logit'),data=pheno3))



t.test(pheno3[x,"covt04_bmi"], pheno3[!x,"covt04_bmi"], na.action=na.omit)

plot(ranef(res))
dev.off()

plot(res)
dev.off()





zScoreList <- lapply(resultList, function(x){return(x[,"z"])})
pvalueList <- lapply(resultList, function(x){return(x[,"p"])})

# combine into z-score matrix excluding intercept
zscores <- do.call("cbind", zScoreList)
pvalues <- do.call("cbind", pvalueList)
#zscores <- zscores[,colnames(zscores) != "(Intercept)"]

str(zscores)

write.table(zscores, file = "zscoreMatrix.txt", sep = "\t", quote = F, col.names = NA)
write.table(pvalues, file = "pvaluesMatrix.txt", sep = "\t", quote = F, col.names = NA)


#below is testingground

plot(resultsPerArray[["Cyto"]][,"zscore"], resultsPerArray[["Gsa"]][,"zscore"])
dev.off()


str(anova(res4,res5))


d <- vragenLong[vragenLong$array == array & !is.na(vragenLong[,q]),]

prsTrait = "BMI_gwas"

prsRange <- quantile(d[,prsTrait],probs = seq(0,1,0.1))
dayRange <- range(d[,"days"])


str(res3)

str(tTable)

coef <- as.matrix(metaRes)[,"y"]
str(coef)
length(coef)


days = seq(dayRange[1], dayRange[2], 1)
days2 = days * days

highPrs <- coef["(Intercept)"] + coef["days"] * days + coef["days2"] * days2 + prsRange[10] * days * coef[paste0(prsTrait,":days")] + coef[prsTrait] * prsRange[10]
lowPrs <- coef["(Intercept)"] + coef["days"] * days + coef["days2"] * days2 + prsRange[1] * days * coef[paste0(prsTrait,":days")] + coef[prsTrait] * prsRange[1]

plot(days, highPrs, ylim = range(lowPrs, highPrs),  col = "red", type = "l", ylab = q, xlab = "Dagen sinds 30 maart 2020", main = paste0("Red is high ", prsTrait, " PRS\nBlue is low PRS"))
points(days, lowPrs, col = "blue", type = "l")
dev.off()

highPrs <- coef["(Intercept)"] + coef["days"] * days + coef["days2"] * days2 + prsRange[10] * days * coef[paste0(prsTrait,":days")]
lowPrs <- coef["(Intercept)"] + coef["days"] * days + coef["days2"] * days2 + prsRange[1] * days * coef[paste0(prsTrait,":days")]

plot(days, highPrs, ylim = range(lowPrs, highPrs),  col = "red", type = "l", ylab = q, xlab = "Dagen sinds 30 maart 2020", main = paste0("Red is high ", prsTrait, " PRS\nBlue is low PRS"))
points(days, lowPrs, col = "blue", type = "l")
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

fixedModel <- as.formula(paste(q, "~((gender_recent+age_recent+age2_recent+household_recent+have_childs_at_home_recent+chronic_recent +", paste(colnames(prs)[-1], collapse = " + ") ,")*days + days2 ) "))

lmFit <- glm(fixedModel,data=vragenLong[vragenLong$array == array,c("PROJECT_PSEUDO_ID", q,colnames(prs)[-1],"gender_recent","age_recent","age2_recent","household_recent","have_childs_at_home_recent","chronic_recent", "days", "days2", "array" )])

#### @Robert fit for normal binom below

q<-qNameMap["everC19Pos",2]
  
#qq<-qNameMap["zijn deze zorgen bijna elke dag aanwezig in de afgelopen 7 dagen?",2]
d <- vragenLong[vragenLong$array == array,c("PROJECT_PSEUDO_ID", q,colnames(prs)[-1],"gender_recent","age_recent","age2_recent","household_recent","have_childs_at_home_recent","chronic_recent", "days", "days2", "array" )]
d[,q] <- as.factor(d[,q])
table(d[,q])
fixedModel <- as.formula(paste(q, "~((gender_recent+age_recent+age2_recent+household_recent+have_childs_at_home_recent+chronic_recent +", paste(colnames(prs)[-1], collapse = " + ") ,")*days + days2  ) "))
glmBinomFit <- glm(fixedModel ,family=binomial(link='logit'),data=d)
summary(glmBinomFit)


prsTrait = "Major.depressive.disorder.in.trauma.exposed.individuals"
prsTrait = "COVID.19.susceptibility"
prsTrait = "Anxiety.tension"
prsTrait = "Schizophrenia"
prsTrait = "EduYears"
prsTrait = "Life.satisfaction"
prsTrait = "BMI_gwas"


plotThreshold = 3.18



str(resultList)
#metaRes <- resultList[["Positive tested cumsum"]]
q <- "Positive.tested.cumsum" 

metaRes <- as.matrix(metaRes)

qName <- names(resultList)[2]



endDate <- startdate + 307

axisAt <- c(startdate,startdate+100,startdate+200, endDate)

  pdf("interactionPlots.pdf", width = 10, height = 6)
  for(qName in names(resultList)){
    q<-qNameMap[qName,2]
    plotEffectsOverTime(resultList[[qName]],q)
  }
  grDevices::dev.off()

qName <- "Positive tested cumsum"
qName <- "hoe waardeert u uw kwaliteit van leven over de afgelopen 14 dagen?"
q<-qNameMap[qName,2]
metaRes <- resultList[[qName]]
effect <- "Life.satisfaction:days"
metaRes<- as.matrix(metaRes)

metaRes <-  summary(res)$tTable
colnames(metaRes)[1] <- "y"

plotEffectsOverTime <- function(metaRes, q){
  metaResSelected <- metaRes[grepl(":days", row.names(metaRes)) & abs(metaRes[,"z"])>=plotThreshold,,drop=F]
  if(nrow(metaResSelected)> 1){
    
    for(effect in rownames(metaResSelected)){
     
      effectName <- sub(":days", "", effect)
      
      if(effectName %in% colnames(prs)){
        
        prsRange <- quantile(prs[,effectName],probs = seq(0,1,0.1))
        
        prsLabel = prsLabels[effectName]
        
        qInfo <- selectedQ[q,]
        
        qLable <- qInfo[,"English.label"]
        
        daysSeq <- qInfo[,"firstDay"]:qInfo[,"lastDay"]
        
        colHigh = "#6300A7"
        colMedium = "#D5546E"
        colLow = "#FCD225"
        colAxis = "grey70"
        colMean = "#808080"
        
        
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
                
        #rpng(width = 1000, height = 800)
        
        layout(matrix(c(1,1,2,3,4,4), nrow=3, byrow = T), heights = c(0.1,0.8,0.1))
        par(mar = c(0,0,0,0), xpd = NA)
        
        plot.new()
        plot.window(xlim = 0:1, ylim = 0:1)
        text(0.5,0.5,paste0("Model fitted on '", qLable, "' stratified by '", prsLabel, "'"), cex = 2 , font = 2)
        
        dummy[,effectName] <- prsRange[10]
        highPrs <- predict_meta(df = dummy, coefficients = coef, family = fam)#, family = binomial(link = "logit")
        dummy[,effectName] <- prsRange[6]
        medianPrs <- predict_meta(df = dummy, coefficients = coef, family = fam)#, family = binomial(link = "logit")
        dummy[,effectName] <- prsRange[2]
        lowPrs <- predict_meta(df = dummy, coefficients = coef, family = fam)#, family = binomial(link = "logit")
        
        p <- predict(res, dummy, level = 0, interval = "confidence",)
        str(p)
        y <- prsRange[c(10,6,2)]
        
        summary(res)$tTable
        
        intervals(res)
        
        
        qnorm(0.95)  * ( 2.280416e-02 + 1.530986e-04 + 4.249592e-07 + 4.717006e-05)
        
        resTest <- res
        
        str(resTest)
        
        library(ggeffects)
        (x <- ggeffect(res, terms = c("days [daysSeq]", "Life.satisfaction [y]") , type = "fixed" ) )
        plot(x)
        dev.off()
        
        gginteraction(res)
        
        ggeffect(res, type = "fixed")
        
        test <- predict(res, newdata = dummy, interval = "confidence", level = 0.95, levels=0)
        str(test)
        #plot(lowPrs, p)
        #dev.off()
        
        rpng()
        par(mar = c(3,5,1,0), xpd = NA)
        plot.new()
        plot.window(xlim = c(startdate,endDate), ylim = range(lowPrs, medianPrs, highPrs))
        axis(side = 1, at = axisAt, labels = format(axisAt, "%d-%b-%Y"), col = colAxis, col.axis = colMean)
        axis(side = 2, col = colAxis, col.axis = colMean)
        title(ylab = paste0("Predicted of ", qLable, "\n by full model"), col.lab = colMean)
        points(startdate +daysSeq, lowPrs, col = colLow, type = "l", lwd = 2)
        points(startdate +daysSeq, medianPrs, col = colMedium, type = "l", lwd = 2)
        points(startdate +daysSeq, highPrs, col = colHigh, type = "l", lwd = 2)
        
        c(startdate +daysSeq, rev(startdate +daysSeq))
        
        ciExtend <- 0.03
        polygon(c(startdate +daysSeq, rev(startdate +daysSeq)), y = c(lowPrs + ciExtend, rev(lowPrs) - ciExtend))
        dev.off()
        
        coef[!grepl(effectName, names(coef))] <- 0
        
        dummy[,effectName] <- prsRange[10]
        highPrs <- predict_meta(df = dummy, coefficients = coef, family = fam)#, family = binomial(link = "logit")
        dummy[,effectName] <- prsRange[6]
        medianPrs <- predict_meta(df = dummy, coefficients = coef, family = fam)#, family = binomial(link = "logit")
        dummy[,effectName] <- prsRange[2]
        lowPrs <- predict_meta(df = dummy, coefficients = coef, family = fam)#, family = binomial(link = "logit")
        
        par(mar = c(3,5,1,0), xpd = NA)
        plot.new()
        plot.window(xlim = c(startdate,endDate), ylim = range(lowPrs, medianPrs, highPrs))
        axis(side = 1, at = axisAt, labels = format(axisAt, "%d-%b-%Y"), col = colAxis, col.axis = colMean)
        axis(side = 2, col = colAxis, col.axis = colMean)
        title(ylab = paste0("Contribution of ", prsLabel, " PGS\non ", qLable, " prediction"), col.lab = colMean)
        points(startdate +daysSeq, lowPrs, col = colLow, type = "l", lwd = 2)
        points(startdate +daysSeq, medianPrs, col = colMedium, type = "l", lwd = 2)
        points(startdate +daysSeq, highPrs, col = colHigh, type = "l", lwd = 2)
        
        
        par(mar = c(0,0,0,0), xpd = NA)
        plot.new()
        plot.window(xlim = 0:1, ylim = 0:1)
        legend("center", fill = c(colLow, colMedium, colHigh), legend = paste0(c("Lowest 10% PGS for ", "Median PGS for ", "Highest 10% PGS for "), prsLabel), bty = "n")
      
        
        dev.off()
        
      }
    }
    
    
    
  }

}

dev.off()
prsRange <- quantile(prs[,prsTrait],probs = seq(0,1,0.1))
  
dummy <- vragenLong[1:307,c(q,colnames(prs)[-1],"gender_recent","age_recent","age2_recent","household_recent","have_childs_at_home_recent","chronic_recent", "days", "days2")]
dummy$days <- 1:307
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

#test <- predict(glmBinomFit, type = "terms", newdata = dummy)
#test[,"Neuroticism:days"]

#coef = coef[,1]
coef=as.matrix(metaRes)[,"y"]
#coef=resultsPerArray[["Gsa"]][,1]

coef["days"] <- 0
coef["days2"] <- 0

#plot(coef)
#dev.off()

#gaussian()

fam = binomial(link = "logit")

dummy[,prsTrait] <- prsRange[10]
highPrs <- predict_meta(df = dummy, coefficients = coef, formula = fixedModel, family = fam)#, family = binomial(link = "logit")
dummy[,prsTrait] <- prsRange[6]
medianPrs <- predict_meta(df = dummy, coefficients = coef, formula = fixedModel, family = fam)#, family = binomial(link = "logit")
dummy[,prsTrait] <- prsRange[2]
lowPrs <- predict_meta(df = dummy, coefficients = coef, formula = fixedModel, family = fam)#, family = binomial(link = "logit")


rpng(width = 800, height = 800)
plot.new()
plot.window(xlim = range(dummy$days), ylim = range(lowPrs, medianPrs, highPrs))
axis(side = 1)
axis(side = 2)
title(xlab = "days", ylab = q, main = prsTrait)
points(lowPrs, col = "blue", type = "l")
points(medianPrs, col = "green", type = "l")
points(highPrs, col = "red", type = "l")
dev.off()

median(pheno3[,"covt01_bmi"], na.rm=T)
median(pheno3[,"covt17_bmi"], na.rm=T)




barplot(table(vragenLong[,q]))
dev.off()

######### @Robert mixed model below

#
modelTest <- as.formula(paste(q, "~((gender_recent+age_recent+age2_recent+household_recent+have_childs_at_home_recent+chronic_recent +", paste(colnames(prs)[-1], collapse = " + ") ,")*days + days2 )+ (1|vl)  + (1|PROJECT_PSEUDO_ID) "))


model <- as.formula(paste(q, "~((gender_recent+age_recent+age2_recent+household_recent+have_childs_at_home_recent+chronic_recent +", paste(colnames(prs)[-1], collapse = " + ") ,")*days + days2 ) + (1|PROJECT_PSEUDO_ID) "))
d <- vragenLong[vragenLong$array == array & !is.na(vragenLong[,q]) & !is.na(vragenLong[,"days"]),c("PROJECT_PSEUDO_ID", q,colnames(prs),"gender_recent","age_recent","age2_recent","household_recent","have_childs_at_home_recent","chronic_recent", "days", "days2", "array", "vl" )]
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


glmMerFit <- glmer(model, data = d2, family = binomial, nAGQ=0 )#glmerControl(optimizer = "nloptwrap2") #control = glmerControl(optimizer = "bobyqa")  
summary(glmMerFit)$coefficients

glmMerFit <- glmer(modelTest, data = d2, family = binomial, nAGQ=0 )#glmerControl(optimizer = "nloptwrap2") #control = glmerControl(optimizer = "bobyqa")  
summary(glmMerFit)$coefficients

glmMerFit2 <- lmer(modelTest, data = d2 )
summary(lmMerFit)



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
ordinalFit <- clmm(modelTest, data = d,na.action=na.omit, link = "logit", nAGQ = 1)

wawares <- lme(fixed=cijfer~((COVID.19.susceptibility)*days), random=~1|PROJECT_PSEUDO_ID, data= vragenLong,na.action=na.omit)
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











library(survival)



q<-qNameMap["everC19Pos",2]

usedPrs <- "Anxiety.tension"
usedPrs <- "COVID.19.susceptibility"
survModel <- as.formula(paste(" Surv(days, ",q,") ~((gender_recent+age_recent+age2_recent+household_recent+have_childs_at_home_recent+chronic_recent +", paste0(usedPrs, collapse = " + ") ,")*days + days2  ) "))

d <- vragenLong[!is.na(vragenLong[,q]) & vragenLong$array == array,c("PROJECT_PSEUDO_ID", q,usedPrs,"gender_recent","age_recent","age2_recent","household_recent","have_childs_at_home_recent","chronic_recent", "days", "days2")]


    
km_fit <- survfit(Surv(days, everC19Pos) ~ have_childs_at_home_recent, data=d)
summary(km_fit, times = c(1,30,60,90*(1:10)))          
  
rpng(width = 800, height = 800)
plot(km_fit, xlab="Days", main = 'C19 incidence', col = 1:2, lwd=2, mark.time=FALSE) 
dev.off()


coxph(formula = survModel, data = d)

x<-vragenLong[x[,"PROJECT_PSEUDO_ID"]=="01682c11-d9a5-4434-a1f6-bf3b418dea30" | x[,"PROJECT_PSEUDO_ID"]=="1523a46-5a67-4b4e-be30-9d79210e0f10",]
str(x)
split(x, x[,"PROJECT_PSEUDO_ID"])

vragenLong[vragenLong[,"everC19Pos"] == 1,"PROJECT_PSEUDO_ID"][10]


vlDates <- read.delim("sendDates.txt", stringsAsFactors = T, row.names = 1)
vlDates[,1] <- as.Date(vlDates[,1], format = "%Y-%m-%d")
vlDates <- vlDates[order(vlDates[,1]),,drop=F]


id="01682c11-d9a5-4434-a1f6-bf3b418dea30"
id="09b31120-caf0-4e3a-9006-31294b1b4411"


lastDate <- max(vragenLong[,qNameMap["responsdatum covid-vragenlijst",2]], na.rm =T)
usedPrs <- colnames(prs)[-1]

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
      censorDate <- lastDate
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

prsTrait = "COVID.19.susceptibility"
prsTrait = "Anxiety.tension"
prsRange <- quantile(prs[,prsTrait],probs = seq(0,1,0.1))
prsRange <- quantile(prs[,prsTrait],probs = c(0,0.1,0.4,0.6,0.9,1))
covid19Events2$prsQuantile <- cut(covid19Events2[,prsTrait],breaks = prsRange, include.lowest = T)
length(table(covid19Events2$prsQuantile))
q<-qNameMap["hebt u een coronavirus/covid-19 infectie (gehad)?",2]

table(covid19Events2[,"eventDays"], useNA="always")
table(covid19Events2[,q], useNA="always")

survModel <- as.formula(paste(" Surv(eventDays, ",q,") ~prsQuantile"))
survModel <- as.formula(paste(" Surv(eventDays, ",q,") ~1"))

km_fit <- survfit(survModel, data=covid19Events2)
summary(km_fit, times = c(1,180,max(covid19Events2$eventDays)))          

colHigh = "firebrick2"
colMedium = "springgreen2"
colLow = "darkturquoise"

prsLabel = prsLabels[prsTrait]

rpng(width = 800, height = 800)
pdf("kmplot.pdf", useDingbats = F)
par(xpd = NA)
plot(km_fit, xlab="Days", main = 'C19 incidence', col = c(colLow, NA,colMedium,NA, colHigh), lwd=2, mark.time=FALSE, ylim = c(0.8, 1)) 
legend("bottomleft", fill = c(colLow, colMedium, colHigh), legend = paste0(c("Lowest 10% PGS for ", "Median PGS for ", "Highest 10% PGS for "), prsLabel), bty = "n")
dev.off()


#usedPrs <- c("COVID.19.susceptibility", "Anxiety.tension")

array="Gsa"
survModel <- as.formula(paste(" Surv(eventDays, ",q,") ~(", paste0(usedPrs, collapse = " + "), ")"))
survModel <- as.formula(paste(" Surv(eventDays, ",q,") ~", paste0(c(confounders,usedPrs), collapse = " + "), "+", paste0("tt(",c(confounders, usedPrs, "days", "days2"),")", collapse = " + ")))
a <- coxph(formula = survModel, data = covid19Events2[covid19Events2$array == "Gsa",], tt=function(x,t,...){as.numeric(x)*t})
b <- coxph(formula = survModel, data = covid19Events2[covid19Events2$array == "Cyto",], tt=function(x,t,...){as.numeric(x)*t})

a
b
summary(a)$coefficients

plot(survfit(a))


veteran
plot(veteran$time, log(veteran$time+20))
dev.off()

vfit <- coxph(Surv(time, status) ~ trt + prior + karno, veteran)
vfit


zp <- cox.zph(vfit)
plot(zp[3]) # a plot for the 3rd variable in the fit
abline(0,0, col=2)
abline(h= vfit$coef[3], col=3, lwd=2, lty=2)
dev.off()


vfit3 <- coxph(Surv(time, status) ~ trt + prior + karno + tt(karno),
              data=veteran,
               tt = function(x, t, ...) x * log(t+20))
vfit3

vfit4 <- coxph(Surv(time, status) ~ trt + prior + karno + tt(karno),
               data=veteran,
               tt = function(x, t, ...) (x * t))
vfit4

vfit5 <- coxph(Surv(time, status) ~ trt + prior + karno + tt(karno),
               data=veteran,
               tt = function(x, t, ...) cbind(x * t,x * t * t ))

vfit5

plot(zp[3])
abline(coef(vfit3)[3:4], col=2)
dev.off()




prsTrait = "COVID.19.susceptibility"
prsTrait = "Anxiety.tension"
prsRange <- quantile(prs[,prsTrait],probs = seq(0,1,0.1))

dummy <- vragenLong[1:307,c(q,"vl", qNameMap["responsdatum covid-vragenlijst",2],usedPrs,"gender_recent","age_recent","age2_recent","household_recent","have_childs_at_home_recent","chronic_recent", "array", "days","days2")]
dummy$days <- 1:307
dummy$eventDays <- 1:307
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
dummy[,q]

dummy[,prsTrait] <- prsRange[10]
dummy[,prsTrait] <- prsRange[6]
dummy[,prsTrait] <- prsRange[1]

str(a)

  csurv1 <- predict(a,newdata=dummy, type = "risk" )
plot(csurv1, xlab="Days")
dev.off()

dummy[,prsTrait] <- prsRange[10]
highPrs <-  predict(a, newdata=dummy, type = type)
dummy[,prsTrait] <- prsRange[6]
medianPrs <- predict(a, newdata=dummy, type = type)
dummy[,prsTrait] <- prsRange[2]
lowPrs <- predict(a, newdata=dummy, type = type)


rpng()
plot.new()
plot.window(xlim = range(dummy$days), ylim = range(lowPrs, medianPrs, highPrs))
axis(side = 1)
axis(side = 2)
title(xlab = "days", ylab = "C19 pos")
points(lowPrs, col = "blue", type = "l")
points(medianPrs, col = "green", type = "l")
points(highPrs, col = "red", type = "l")
dev.off()



predict.coxph

mean(pheno3[pheno3[,"array"] == "Gsa","age_recent"])
mean(pheno3[pheno3[,"array"] == "Cyto","age_recent"])

layout(matrix(1:2,nrow = 2))
hist(pheno3[pheno3[,"array"] == "Gsa","age_recent"], breaks = 30)
hist(pheno3[pheno3[,"array"] == "Cyto","age_recent"], breaks = 30)
dev.off()





table(pheno3[pheno3[,"array"] == "Gsa","gender_recent"])
table(pheno3[pheno3[,"array"] == "Cyto","gender_recent"])

table(pheno3[pheno3[,"array"] == "Gsa","household_recent"])
table(pheno3[pheno3[,"array"] == "Cyto","household_recent"])


table(pheno3[pheno3[,"array"] == "Gsa","have_childs_at_home_recent"])
table(pheno3[pheno3[,"array"] == "Cyto","have_childs_at_home_recent"])



inclusionPerVl <- read.delim("/groups/umcg-lifelines/tmp01/projects/ov20_0554/analysis/risky_behaviour/PRS_correlation/inclusionPerVl.txt", row.names =1)
inclusionPerVl <- inclusionPerVl[pheno3$PROJECT_PSEUDO_ID ,]
pdf("completedQperPerson.pdf")
par(xpd = NA)
barplot(table(apply(inclusionPerVl,1,sum)), main = "Number of completed questionnaires per participant", xlab = "Completed questionnaires", ylab = "Number of particpants", las = 2)
dev.off()


