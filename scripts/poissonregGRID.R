### Regression Analysis of SCBI Camera Grid Data
getwd()
deerdata <- read.csv("Deer_DF.csv", header = T)
require(AICcmodavg)
source("N:/Capacity Building and Academic Programs/SI-Mason Grad & Prof Training/Individual Courses--Folders/Stats for Ecol & Cons Bio/2018/CourseMaterials/R_scripts/R_scripts/functions/pairsPannelFunctions.r")
pairs(deerdata[,c(6:12)], diag.panel = panel.hist, lower.panel = panel.smooth, upper.panel = panel.cor)
str(deerdata)

hist(deerdata$Num_Stems)
hist(log10(deerdata$Num_Stems))
hist(deerdata$EDD)
hist(deerdata$Cam_Nights, breaks = 5)

glm.po1 <- glm(Number_of_Sequences ~ Log + EDD + Cam_Height + log10(Num_Stems) + Num_Oaks +offset(log(Cam_Nights)), data = deerdata, family = poisson)
summary(glm.po1)

k<- log(nrow(deerdata))
drop1(glm.po1, k=k) #using BIC

glm.po2 <- update(glm.po1, ~. -Cam_Height)
summary(glm.po2)
drop1(glm.po2, k=k)

#test for overdispersion
sum.po2 <- summary(glm.po2)
phi.po2 <- sum.po2$deviance/sum.po2$df.residual 

#____________________________________
#Response Curves for our Covariates
#____________________________________

meanEDD <- mean(deerdata$EDD)
medoak <- median(deerdata$Num_Oaks)
medstems <- median(deerdata$Num_Stems)

#Response Curves for Number of Stems. Will use average EDD, and median # of oaks. We will draw graphs for logs present and for logs absent

seq.stems <- seq(min(deerdata$Num_Stems), max(deerdata$Num_Stems), length.out = 100)
nd.seq.stem.log <- data.frame(Cam_Nights = 61, EDD = meanEDD, Num_Oaks = medoak, Log = "YES", Num_Stems = seq.stems)
pred.stems.l <- predict(glm.po2, newdata = nd.seq.stem.log, se=TRUE, type = "link")
plot(Number_of_Sequences ~ Num_Stems, data=deerdata, pch=16, col="orange", las=1, log = "x", xlab = "log(Num_Stems)", ylab = "# Deer of Sequences")
lines(nd.seq.stem.log$Num_Stems, exp(pred.stems$fit), col="purple")
lines(nd.seq.stem.log$Num_Stems, exp(pred.stems$fit + 1.96*pred.stems$se.fit), lty=2, col="purple")
lines(nd.seq.stem.log$Num_Stems, exp(pred.stems$fit - 1.96*pred.stems$se.fit), lty=2, col="purple")

nd.seq.stem.nolog <- data.frame(Cam_Nights = 61, EDD = meanEDD, Num_Oaks = medoak, Log = "NO", Num_Stems = seq.stems)
pred.stems.nl <- predict(glm.po2, newdata = nd.seq.stem.nolog, se=TRUE, type = "link")
lines(nd.seq.stem.nolog$Num_Stems, exp(pred.stems.nl$fit), col="green")
lines(nd.seq.stem.nolog$Num_Stems, exp(pred.stems.nl$fit + 1.96*pred.stems.nl$se.fit), lty=2, col="green")
lines(nd.seq.stem.nolog$Num_Stems, exp(pred.stems.nl$fit - 1.96*pred.stems.nl$se.fit), lty=2, col="green")

#Response Curve for Effective Detection Distance. Will use median # of stems, median # of oaks, and will show with both log and without.
#Think about labeling these two outliers points so we can see what cameras they are
seq.EDD <- seq(min(deerdata$EDD), max(deerdata$EDD), length.out = 100)
nd.seq.EDD.log <- data.frame(Cam_Nights = 61, EDD = seq.EDD, Num_Oaks = medoak, Log = "YES", Num_Stems = medstems)
pred.EDD.l <- predict(glm.po2, newdata = nd.seq.EDD.log, se=TRUE, type = "link")
plot(Number_of_Sequences ~ EDD, data=deerdata, pch=16, col="orange", las=1, xlab = "Effective Detection Distance (m)", ylab = "# of Deer Sequences")
lines(nd.seq.EDD.log$EDD, exp(pred.EDD.l$fit), col="purple")
lines(nd.seq.EDD.log$EDD, exp(pred.EDD.l$fit + 1.96*pred.EDD.l$se.fit), lty=2, col="purple")
lines(nd.seq.EDD.log$EDD, exp(pred.EDD.l$fit - 1.96*pred.EDD.l$se.fit), lty=2, col="purple")

nd.seq.EDD.nolog <- data.frame(Cam_Nights = 61, EDD = seq.EDD, Num_Oaks = medoak, Log = "NO", Num_Stems = medstems)
pred.EDD.nl <- predict(glm.po2, newdata = nd.seq.EDD.nolog, se=TRUE, type = "link")
lines(nd.seq.EDD.nolog$EDD, exp(pred.EDD.nl$fit), col="green")
lines(nd.seq.EDD.nolog$EDD, exp(pred.EDD.nl$fit + 1.96*pred.EDD.nl$se.fit), lty=2, col="green")
lines(nd.seq.EDD.nolog$EDD, exp(pred.EDD.nl$fit - 1.96*pred.EDD.nl$se.fit), lty=2, col="green")
