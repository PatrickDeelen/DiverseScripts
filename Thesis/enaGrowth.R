setwd("C:\\Users\\patri\\Dropbox\\UMCG\\Thesis")

ena <- read.delim("./Chapter 9/ena_2018_08_13.txt", stringsAsFactors = F)

str(ena)

library(zoo)
ena$DATEQ <- as.factor(as.yearqtr(as.Date(ena$first_public, format = "%Y-%m-%d")))


qCount <- table(factor(ena$DATEQ[ ena$DATEQ != "NA QNA" & ena$DATEQ != "2018 Q3"]))#

qCount2 <- qCount

for(i in 2:length(qCount2)){
  qCount2[i] <- qCount2[i] + qCount2[i-1]
}

x<- (1:length(qCount2))
#xExp <- exp(x)

fit3 <- nls(qCount2 ~ a * x + exp(a + b * x), start = list(a = 1, b = 1))
cor.test(qCount2, predict(fit3, list(x=x)))

plot(as.numeric(qCount2), predict(fit3, list(x=x)))
abline(0,1)

summary(fit3)


pdf("./Chapter 9/enaGrowth.pdf", height = 5, width = 7)
par(xpd = NA, las = 2)

barplot(qCount2)
lines(seq(from = 0.5, to = 40.3, by = 1.2), predict(fit3, list(x = x)), col = "red", lwd = 1, xpd = NA)
dev.off()


fit = lm(log(qCount2) ~ log(x))
plot(as.numeric(qCount2))
lines(x, x ^ fit$coefficients[2], col = "red")
