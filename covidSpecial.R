#remoter::server(verbose = T, port = 55556, password = "laberkak", sync = T)

remoter::client("localhost", port = 55556, password = "laberkak")


setwd("/groups/umcg-lifelines/tmp01/projects/ov20_0554/")

library(data.table)
covidRes <- fread("./data/raw/covid_questionnaires/week10/covid19-week10-labels.dat", stringsAsFactors = F)

covidResWeek6 <- fread("./data/raw/covid_questionnaires/week6/rl02/covid19-week6-2-labels.dat", stringsAsFactors = F)
covidResWeek8 <- fread("./data/raw/covid_questionnaires/week8/covid19-week8-labels.dat", stringsAsFactors = F)

covidLab <- fread("./data/raw/covid_questionnaires/week10/covid19-week10_M.dat", stringsAsFactors = F)
covidLabWeek6 <- fread("./data/raw/covid_questionnaires/week6/covid19-week6_M.dat", stringsAsFactors = F)
covidLabWeek8 <- fread("./data/raw/covid_questionnaires/week8/covid19-week8_M.dat", stringsAsFactors = F)

covidResWeek6_B <- merge(covidRes[,1,drop=F],covidResWeek6, by = "PSEUDOIDEXT", all.x = T, all.y=F)
covidResWeek6_B <- covidResWeek6_B[match(covidRes$PSEUDOIDEXT, covidResWeek6_B$PSEUDOIDEXT)]

covidResWeek8_B <- merge(covidRes[,1,drop=F],covidResWeek8, by = "PSEUDOIDEXT", all.x = T, all.y=F)
covidResWeek8_B <- covidResWeek8_B[match(covidRes$PSEUDOIDEXT, covidResWeek8_B$PSEUDOIDEXT)]


demo = fread("/groups/umcg-lifelines/tmp01/projects/ov20_0554/analysis/demographics/Participant (Pa_99_G).dat",header=T,data.table=F)

covidRes$gender <- factor(demo$GESLACHT[match(covidRes$PSEUDOIDEXT, demo$PSEUDOIDEXT)], levels = c(1,2), labels = c("man", "vrouw"))

all(covidResWeek6_B$PSEUDOIDEXT == covidRes$PSEUDOIDEXT )
all(covidResWeek8_B$PSEUDOIDEXT == covidRes$PSEUDOIDEXT )

str(covidRes)

query <- function(q){
  print(covidLab[covidLab$PSEUDOIDEXT==q,2])
  print(table(covidRes[[q]]))
  
  
  
  x <- table(factor(covidRes[[q]], exclude = c("8888", "9999")))
  
  print((x *100)/sum(x))

  
}

queryW6 <- function(q){
  print(covidLabWeek6[covidLabWeek6$PSEUDOIDEXT==q,2])
  print(table(covidResWeek6_B[[q]]))
  x <- table(factor(covidResWeek6_B[[q]], exclude = c("8888", "9999")))
  print((x *100)/sum(x))
}

queryW8 <- function(q){
  print(covidLabWeek8[covidLabWeek8$PSEUDOIDEXT==q,2])
  print(table(covidResWeek8_B[[q]]))
  x <- table(factor(covidResWeek8_B[[q]], exclude = c("8888", "9999")))
  print((x *100)/sum(x))
}


  query("COVID75")

table(covidRes$COVID85)
table(covidRes$COVID75)

query("COVID85")
query("COVID87H")
x <- table(factor(covidRes[["COVID85"]][covidRes[["COVID87H"]]!="Nee"], exclude = c("8888", "9999")))
print((x *100)/sum(x))

x <- table(factor(covidRes[["COVID85"]][covidRes[["COVID87H"]]=="Ja, naar het buitenland" | covidRes[["COVID87H"]]=="Ja, naar het buitenland én in Nederland"], exclude = c("8888", "9999")))
print((x *100)/sum(x))



query("COVID86")

query("COVID87C")

x <- table(factor(covidRes[["COVID87C"]], exclude = c("8888", "9999")))
y <- (x *100)/sum(x)
y["Ja, voor 1 september"] + y["Ja, zowel voor als na 1 september"]

query("COVID87D")
sort(table(tolower(covidRes$COVID87D)))

keysNederland <- tolower(c("nederland", "holland", "eigen", "drenthe", "flevoland", "Friesland", "Gelderland", "Groningen", "Limburg", "Noord Brabant", "Noord Holland", "Overijssel", "Zuid Holland", "Utrecht", "Zeeland", "Texel", "Vlieland", "Terschelling", "Ameland", "Schiermonnikoog"))

zomervakantieGangers <- covidRes[["COVID87C"]] == "Ja, voor 1 september" | covidRes[["COVID87C"]] == "Ja, zowel voor als na 1 september"

z <- lapply(keysNederland, agrep, x = tolower(covidRes$COVID87D[zomervakantieGangers]))

thuisBlijvers <- unique(do.call(c, z))

length(thuisBlijvers) * 100 / sum(zomervakantieGangers)





sort(table(tolower(covidRes$COVID87D[zomervakantieGangers][thuisBlijvers])))

write.table(sort(table(tolower(covidRes$COVID87D[zomervakantieGangers][-thuisBlijvers]))), file = "buitenlandVakantie.txt", sep = "\t", quote = F, row.names = F)



x <- sort(table(factor(tolower(covidRes$COVID87E_O9czmy[zomervakantieGangers][-thuisBlijvers]), exclude = c("8888", "9999"))))
y <- (x *100)/sum(x)

query("COVID87E_O9czmy")


landCodes <- read.delim("country_codes_v3.txt", sep = ",", stringsAsFactors = F)

risicoLanden <- landCodes$country[landCodes$code == "rood" | landCodes$code == "oranje"]

z <- sapply(risicoLanden, agrep, x = tolower(covidRes$COVID87D[zomervakantieGangers]), max.distance =0.05)
risicoVakanties <- do.call(c, z)


z2 <- sapply(z, function(x){return(covidRes$COVID87D[zomervakantieGangers][x])})

risicoVakanties2 <- do.call(c, z2)


matches <- cbind(risicoVakanties, risicoVakanties2)

write.table(matches, file = "landMatch.txt", quote = F, sep = "\t")

risicoLanden[67]




  COVID76 <- factor(covidRes[["COVID76"]], exclude = c("8888", "9999"), levels = c("Minder dan 3 maanden", "3 maanden", "6 maanden", "1 jaar", "1,5 jaar", "2 jaar", "3 jaar", "Langer dan 3 jaar", "Weet ik niet"))

x <- table(COVID76)
y <- (x *100)/sum(x)

rpng(width = 800, height = 800)
par(mar = c(12, 4, 4, 2) + 0.1)
barplot(y, ylab = "Percentage deelnemers", las = 3, main = "Hoe lang denkt u dat de pandemie nog gaat duren?")
dev.off()

query("COVID81")
COVID81 <- factor(covidRes[["COVID81"]], exclude = c("8888", "9999"), levels = c("Minder dan 3 maanden", "3 maanden", "6 maanden", "1 jaar", "1,5 jaar", "2 jaar", "3 jaar", "Langer dan 3 jaar", "Weet ik niet"))

x <- table(COVID81)
y <- (x *100)/sum(x)

rpng(width = 800, height = 800)
par(mar = c(12, 4, 4, 2) + 0.1)
barplot(y, ylab = "Percentage deelnemers", las = 3, main = "Hoe lang denkt u dat het duurt voordat er een vaccin komt?")
dev.off()

rpng(width = 800, height = 800)
par(mar = c(12, 12, 4, 2) + 0.1, las = 2)
plot(COVID76, COVID81, xlab = "")
dev.off()


cor.test(as.numeric(COVID76), as.numeric(COVID81))
plot(as.numeric(COVID76), as.numeric(COVID81))

map <- table(COVID76, COVID81)

library(heatmap3)
library(RColorBrewer)

rpng(width = 800, height = 800)
par(mar = c(12, 2, 2, 12) + 0.1)
heatmap3(log2(map), scale= "none", Rowv = NA, Colv = NA, col = colorRampPalette( brewer.pal(9, "YlOrBr"))(100), xlab = "Hoe lang denkt u dat de pandemie nog gaat duren?
", ylab = "Hoe lang denkt u dat het duurt voordat er een vaccin komt?", margins = c(15,15))
dev.off()


query("COVID82")
x <- table(factor(covidRes[["COVID82"]], exclude = c("8888", "9999")))
y <- (x *100)/sum(x)
pie(y, main = "Laat u zich vaccineren als het\nvaccin tegen COVID19 beschikbaar is?")
dev.off()




query("COVID78")
x <- table(factor(covidRes[["COVID78"]], exclude = c("8888", "9999")))
y <- (x *100)/sum(x)
pie(y, main = "Verwacht u een tweede golf aan corona infecties?")
dev.off()



query("COVID79")
x <- table(factor(covidRes[["COVID79"]], exclude = c("8888", "9999")))
y <- (x *100)/sum(x)
rpng(width = 800, height = 800)
par(mar = c(12, 2, 12, 2) + 0.1)
pie(y, main = "Als er een tweede golf komt,\nwelk scenario heeft uw voorkeur?")
dev.off()



query("COVID72A")
query("COVID72B")

rpng(width = 800, height = 800)
layout(matrix(1:2, nrow =2))
x <- table(factor(covidRes[["COVID72A"]], exclude = c("8888", "9999")))
y <- (x *100)/sum(x)

pie(y, main = "Hebt u de intentie u te laten testen als u lichte klachten hebt")

x <- table(factor(covidRes[["COVID72B"]], exclude = c("8888", "9999")))
y <- (x *100)/sum(x)

pie(y, main = "Hebt u de intentie u te laten testen als u ernstige klachten hebt")


dev.off()


query("COVID88")
query("COVID88B")

sort(table(tolower(covidRes$COVID88)))



library(wordcloud2) 
library(webshot)
library("htmlwidgets")

misString <- gsub("\\.", "", tolower(covidRes$COVID88))
misString <- gsub(",", "", misString)
misString <- gsub("de\\s", "", misString, perl=TRUE)
misString <- gsub("het\\s", "", misString, perl=TRUE)
misString <- gsub("\\?+", "", misString)
misString <- gsub("/", "", misString)
misString <- gsub("-", "", misString)
misString <- gsub("_", "", misString)
misString <- trimws(misString)


x <- table(factor(misString, exclude = c("8888", "9999", "nee", "weet ik niet", "niet veel", "0", "000", "niets", "-", "x", "weet niet", "niet", "", "geen idee", "niks", "?", "ik ga niets missen", "ik mis niets", "nvt", "ik ga niks missen", "ik mis niks", "ik zou niet weten", "eigenlijk niks", "eigenlijk niets", "weinig")))
y <- (x *100)/sum(x)
sort(y)

write.table(x, "cloudCountMissen.txt", sep = "\t", quote = F, row.names = F)
write.table(y, "cloudPercentageMissen.txt", sep = "\t", quote = F, row.names = F)


    sort(x)
cloud <- wordcloud2(x[x>20], shape = "cardioid", minSize =3)
saveWidget(cloud,"tmp.html",selfcontained = F)
webshot("tmp.html","fig_1.pdf", delay =20)



query("COVID91B1")
x <- covidRes[["COVID91B1"]] != "Nooit" | covidRes[["COVID91B1"]] != "9999" | covidRes[["COVID91B1"]] != "8888"
table(covidRes[["COVID90C_Opanxh"]][x])

table(covidRes[["COVID90A"]][x])

query("COVID96O")
x <- table(factor(covidRes[["COVID96O"]], exclude = c("8888", "9999")))
y <- (x *100)/sum(x)

rpng(width = 800, height = 800)

layout(matrix(1:2, nrow =2))

x <- table(factor(covidRes[["COVID96O"]], exclude = c("8888", "9999")))
y <- (x *100)/sum(x)
pie(y, main = "Dragen van een mondkapje in het openbaar vervoer")

x <- table(factor(covidRes[["COVID100O"]], exclude = c("8888", "9999")))
y <- (x *100)/sum(x)
pie(y, main = "Volhouden dragen van een mondkapje in het openbaar vervoer")


dev.off()


q <- "COVID99G"
query(q)
rpng(width = 800, height = 800)
x <- table(factor(covidRes[[q]], exclude = c("8888", "9999")))
y <- (x *100)/sum(x)
pie(y, main = covidLab[covidLab$PSEUDOIDEXT==q,2])
dev.off()



q <- "COVID98_O2nxzv"
query(q)
rpng(width = 800, height = 800)
x <- table(factor(covidRes[[q]], exclude = c("8888", "9999")))
y <- (x *100)/sum(x)
pie(y, main = "Wanneer moeten de COVID19 maatregelen weer strenger worden?\n(Als het coronavirus reproductiegetal groter is dan 1")
dev.off()

covidRes$COVID98_O2nxzv_COVID98_O95r9w <- covidRes$COVID98_O2nxzv
covidRes$COVID98_O2nxzv_COVID98_O95r9w[ covidRes$COVID98_O95r9w == "Als het coronavirus reproductiegetal groter is dan 2. Het reproductie getal geeft aan hoeveel andere personen een corona" ] <- "Als het coronavirus reproductiegetal groter is dan 1. Het reproductie getal geeft aan hoeveel andere personen een corona"

q <- "COVID98_O2nxzv_COVID98_O95r9w"
rpng(width = 800, height = 800)
x <- table(factor(covidRes[[q]], exclude = c("8888", "9999")))
y <- (x *100)/sum(x)
pie(y, main = "Wanneer moeten de COVID19 maatregelen weer strenger worden?\n(Als het coronavirus reproductiegetal groter is dan 2")
dev.off()




q <- "COVID98_Oe2njr"
query(q)
rpng(width = 800, height = 800)
x <- table(factor(covidRes[[q]], exclude = c("8888", "9999")))
y <- (x *100)/sum(x)
pie(y, main = covidLab[covidLab$PSEUDOIDEXT==q,2])
dev.off()

q <- "COVID98_Orpyjv"
query(q)
rpng(width = 800, height = 800)
x <- table(factor(covidRes[[q]], exclude = c("8888", "9999")))
y <- (x *100)/sum(x)
pie(y, main = covidLab[covidLab$PSEUDOIDEXT==q,2])
dev.off()

q <- "COVID98_Oshtms"
query(q)
rpng(width = 800, height = 800)
x <- table(factor(covidRes[[q]], exclude = c("8888", "9999")))
y <- (x *100)/sum(x)
pie(y, main = covidLab[covidLab$PSEUDOIDEXT==q,2])
dev.off()

covidRes$COVID98_Oshtms_2 <- 0
covidRes$COVID98_Oshtms_2[covidRes$COVID98_Oe2njr == "Als 25% van de IC bedden vol liggen"] <- 1
covidRes$COVID98_Oshtms_2[covidRes$COVID98_Orpyjv == "Als 50% van de IC bedden vol liggen"] <- 1
covidRes$COVID98_Oshtms_2[covidRes$COVID98_Oshtms == "Als 75% van de IC bedden vol liggen"] <- 1

rpng(width = 800, height = 800)
  x <- table(factor(covidRes[["COVID98_Oshtms_2"]], exclude = c("8888", "9999")))
y <- (x *100)/sum(x)
pie(y, main = "75% IC vol")
dev.off()




query("COVID81")


covidRes$COVID81


x <- table(factor(covidRes[["COVID81"]], exclude = c("8888", "9999", "Weet ik niet")))
y <- (x *100)/sum(x)



x <- table(factor(covidRes[["COVID99H"]], exclude = c("8888", "9999", "Niet van toepassing")))
y <- (x *100)/sum(x)
y

covidRes$COVID14A1[covidRes$COVID14A1 > 30] <- NA
covidRes$COVID14B1[covidRes$COVID14B1 > 30] <- NA
covidResWeek6_B$COVID14A1[covidResWeek6_B$COVID14A1 > 30] <- NA
covidResWeek6_B$COVID14B1[covidResWeek6_B$COVID14B1 > 30] <- NA

covidRes$COVID14A1[covidRes$COVID14A1 < 0] <- NA
covidRes$COVID14B1[covidRes$COVID14B1  < 0] <- NA
covidResWeek6_B$COVID14A1[covidResWeek6_B$COVID14A1 < 0] <- NA
covidResWeek6_B$COVID14B1[covidResWeek6_B$COVID14B1 < 0] <- NA 

sum((covidRes$COVID14A1 > 0 & covidRes$COVID14A1< 25) | (!is.na(covidResWeek6_B$COVID14A1) & covidResWeek6_B$COVID14A1 > 0 & covidResWeek6_B$COVID14A1< 25))

sum((covidRes$COVID14B1 > 0 & covidRes$COVID14B1< 25) | (!is.na(covidResWeek6_B$COVID14B1) & covidResWeek6_B$COVID14B1 > 0 & covidResWeek6_B$COVID14B1< 25))

covidRes$COVID14A1_B <- pmax(covidRes$COVID14A1,covidResWeek6_B$COVID14A1, na.rm=T)
covidRes$COVID14B1_B <- pmax(covidRes$COVID14B1,covidResWeek6_B$COVID14B1, na.rm=T)

covidRes$COVID14A1_COVID14B1 <-  covidRes$COVID14A1_B + covidRes$COVID14B1_B

hist(covidRes$COVID14A1_COVID14B1)
dev.off()

dim(covidRes)
dim(covidResWeek6_B)

sum(covidResWeek6$COVID14A1 > 0 & covidResWeek6$COVID14A1< 25)


queryW6("COVID14A1")
queryW6("COVID14B1")
queryW6("COVID15")

onderwijsJongeKinderen <- covidRes$COVID99G[(covidRes$COVID14A1 > 0 & covidRes$COVID14A1< 25) | (!is.na(covidResWeek6_B$COVID14A1) & covidResWeek6_B$COVID14A1 > 0 & covidResWeek6_B$COVID14A1< 25)]
x <- table(factor(onderwijsJongeKinderen, exclude = c("8888", "9999", "Niet van toepassing")))
y <- (x *100)/sum(x)
rpng(width = 800, height = 800)
par(mar = c(12, 2, 12, 2) + 0.1)
pie(y, main = "Organiseren van thuisonderwijs voor mijn kinderen\n(ouders kinderen <= 12")
dev.off()

onderwijsOudereKinderen <- covidRes$COVID99G[(covidRes$COVID14B1 > 0 & covidRes$COVID14B1< 25) | (!is.na(covidResWeek6_B$COVID14B1) & covidResWeek6_B$COVID14B1 > 0 & covidResWeek6_B$COVID14B1< 25)]
x <- table(factor(onderwijsOudereKinderen, exclude = c("8888", "9999", "Niet van toepassing")))
y <- (x *100)/sum(x)
rpng(width = 800, height = 800)
par(mar = c(12, 2, 12, 2) + 0.1)
pie(y, main = "Organiseren van thuisonderwijs voor mijn kinderen\n(ouders kinderen > 12")
dev.off()



onderwijsMeerdereKinderen <- covidRes$COVID99G[!is.na(covidRes$COVID14A1_COVID14B1) & covidRes$COVID14A1_COVID14B1 >=3 ]
x <- table(factor(onderwijsMeerdereKinderen, exclude = c("8888", "9999", "Niet van toepassing")))
y <- (x *100)/sum(x)
rpng(width = 800, height = 800)
par(mar = c(12, 2, 12, 2) + 0.1)
pie(y, main = "Organiseren van thuisonderwijs meerdere kinderen")
dev.off()

rpng(width = 800, height = 800)
layout(matrix(1:2, nrow =2))


onderwijsMeerdereKinderenMan <- covidRes$COVID99G[!is.na(covidRes$COVID14A1_COVID14B1) & covidRes$COVID14A1_COVID14B1 > 0 & covidRes$gender == "man"]
x <- table(factor(onderwijsMeerdereKinderenMan, exclude = c("8888", "9999", "Niet van toepassing")))
y <- (x *100)/sum(x)
pie(y, main = "Organiseren van thuisonderwijs man")

onderwijsMeerdereKinderenVrouw <- covidRes$COVID99G[!is.na(covidRes$COVID14A1_COVID14B1) & covidRes$COVID14A1_COVID14B1 > 0 & covidRes$gender == "vrouw"]
x <- table(factor(onderwijsMeerdereKinderenVrouw, exclude = c("8888", "9999", "Niet van toepassing")))
y <- (x *100)/sum(x)
pie(y, main = "Organiseren van thuisonderwijs vrouw")
dev.off()

covidResWeek8_B$COVID68H4[covidResWeek8_B$COVID68H4 == "8888"] <- NA
covidResWeek8_B$COVID68H4[covidResWeek8_B$COVID68H4 == "9999"] <- NA


onderwijsMeerdereKinderenSlechteBalance <- covidRes$COVID99G[!is.na(covidResWeek8_B$COVID68H4) & (covidResWeek8_B$COVID68H4 == "Erg slecht" | covidResWeek8_B$COVID68H4 == "Slecht")]
x <- table(factor(onderwijsMeerdereKinderenSlechteBalance, exclude = c("8888", "9999", "Niet van toepassing")))
y <- (x *100)/sum(x)
rpng(width = 800, height = 800)
par(mar = c(12, 2, 12, 2) + 0.1)
pie(y, main = "Organiseren van thuisonderwijs\nslecht of erg slechte prive werk balans in april")
dev.off()

queryW8("COVID68H4")


onderwijsGeenKinderen <- covidRes$COVID99G[!((covidRes$COVID14B1 > 0 & covidRes$COVID14B1< 25) | (!is.na(covidResWeek6_B$COVID14B1) & covidResWeek6_B$COVID14B1 > 0 & covidResWeek6_B$COVID14B1< 25)) & !((covidRes$COVID14A1 > 0 & covidRes$COVID14A1< 25) | (!is.na(covidResWeek6_B$COVID14A1) & covidResWeek6_B$COVID14A1 > 0 & covidResWeek6_B$COVID14A1< 25))]
x <- table(factor(onderwijsGeenKinderen, exclude = c("8888", "9999")))
y <- (x *100)/sum(x)
rpng(width = 800, height = 800)
par(mar = c(12, 2, 12, 2) + 0.1)
pie(y, main = "Organiseren van thuisonderwijs voor mijn kinderen\n(deelnemers zonder thuiswonende kinderen)")
dev.off()


pdf("bereidheidNieuweLockdown.pdf")
par(mar = c(12, 2, 12, 2) + 0.1)
bereidheidNieuweLockdown <- sapply(paste0("COVID99", LETTERS[1:8]), function(q){
  x <- table(factor(covidRes[[q]], exclude = c("8888", "9999", "Niet van toepassing"), levels = c("Helemaal mee oneens", "Oneens", "Een beetje oneens", "Niet mee eens, niet mee oneens", "Een beetje mee eens", "Mee eens", "Helemaal mee eens")))
  y <- (x *100)/sum(x)
  pie(y, main = covidLab[covidLab$PSEUDOIDEXT==q,2])
  return(y)  
})
dev.off()
write.table(bereidheidNieuweLockdown, file = "bereidheidNieuweLockdown.txt", sep = "\t", quote = F)

q <- "G"
query(q)
rpng(width = 800, height = 800)
x <- table(factor(covidRes[[q]], exclude = c("8888", "9999", "Niet van toepassing"), levels = c("Helemaal mee oneens", "Oneens", "Een beetje oneens", "Niet mee eens, niet mee oneens", "Een beetje mee eens", "Mee eens", "Helemaal mee eens")))
y <- (x *100)/sum(x)
pie(y, main = covidLab[covidLab$PSEUDOIDEXT==q,2])
dev.off()

print(levels(factor(covidRes[[q]], exclude = c("8888", "9999", "Niet van toepassing"))), sep = "\n")
