


library(data.table)
covidRes <- fread("/groups/umcg-lifelines/tmp01/projects/ov20_0554/data/raw/covid_questionnaires/week20/covid19-week20-labels.dat", quote="", na.strings = c('8888','9999','NA'))
participantsData <- fread("/groups/umcg-lifelines/tmp01/releases/pheno_lifelines/v1/tab_separated_labels/Participant (Pa_99_G).dat")
covidRes2 <- merge(covidRes, participantsData, by = "PSEUDOIDEXT")

str(participantsData)
qMap <- read.delim("/groups/umcg-lifelines/tmp01/projects/ov20_0554/data/raw/covid_questionnaires/week20/covid19-week20_M.txt", col.names = c("qId","qName"))

str(qMap)


covidRes2$ageGroup <- cut(covidRes2$AGE_QUESTIONNAIRE,breaks = seq(10,100,20))


x <- table(covidRes2$COVID269B, covidRes2$ageGroup)
round(apply(x, 2, function(x){x/sum(x)*100}),2)[c("Helemaal niet mee eens", "Niet mee eens", "Niet mee eens, niet mee oneens",  "Mee eens" ,  "Helemaal mee eens"   ),]
print(unique(covidRes2$COVID269B), sep = ", ", quote = T)



