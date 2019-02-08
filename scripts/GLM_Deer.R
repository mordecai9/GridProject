###Regression Analysis of SCBI Camera Grid Data

require(AICcmodavg)
require(MASS) #for glm.nb function
source("scripts/pairsPannelFunctions.r")


# Data Import and Exploration - WTD ---------------------------------------

load("data/deerDataSum.RData")
load("data/deerDataFall.RData")
load("data/deerDataWin.RData")
load("data/deerDataSpr.RData")

pairs(deerDataSum[,c(4,7,13,14,15,17)], diag.panel = panel.hist, lower.panel = panel.smooth, upper.panel = panel.cor)
str(deerDataSum)

hist(deerDataSum$Num_Stems)
hist(log10(deerDataSum$Num_Stems)) #because one site has so many stems, we need to log10 transform this parameter for all models
hist(deerDataSum$Summer.Fall.EDD)

#Univariate exploratory plots for each of our 5 covariates
par.default <- par(no.readonly = T)
par(mfrow = c(2,2))
plot(nSeqs ~ Summer.Fall.EDD, data = deerDataSum, main = "Summer")
plot(nSeqs ~ Summer.Fall.EDD, data = deerDataFall, main = "Fall")
plot(nSeqs ~ Winter.Spring.EDD, data = deerDataWin, main = "Winter")
plot(nSeqs ~ Winter.Spring.EDD, data = deerDataSpr, main = "Spring")

plot(nSeqs ~ OakDBH, data = deerDataSum, main = "Summer")
plot(nSeqs ~ OakDBH, data = deerDataFall, main = "Fall")
plot(nSeqs ~ OakDBH, data = deerDataWin, main = "Winter")
plot(nSeqs ~ OakDBH, data = deerDataSpr, main = "Spring")

plot(nSeqs ~ log10(Num_Stems), data = deerDataSum, main = "Summer")
plot(nSeqs ~ log10(Num_Stems), data = deerDataFall, main = "Fall")
plot(nSeqs ~ log10(Num_Stems), data = deerDataWin, main = "Winter")
plot(nSeqs ~ log10(Num_Stems), data = deerDataSpr, main = "Spring")

plot(nSeqs ~ Height_cm, data = deerDataSum, main = "Summer")
plot(nSeqs ~ Height_cm, data = deerDataFall, main = "Fall")
plot(nSeqs ~ Height_cm, data = deerDataWin, main = "Winter")
plot(nSeqs ~ Height_cm, data = deerDataSpr, main = "Spring")

plot(nSeqs ~ Log.in.View, data = deerDataSum, main = "Summer")
plot(nSeqs ~ Log.in.View, data = deerDataFall, main = "Fall")
plot(nSeqs ~ Log.in.View, data = deerDataWin, main = "Winter")
plot(nSeqs ~ Log.in.View, data = deerDataSpr, main = "Spring")
par(par.default)

hist(deerDataSum$Deploy.Duration, breaks = 5)


# Full Poisson GLM and Overdispersion -------------------------------------

#Exploratory Poisson model for Summer data to see if a negative binomial model is needed due to potential overdispersion
glm.po.fullS <- glm(nSeqs ~ Log.in.View + Summer.Fall.EDD + Height_cm + log10(Num_Stems) + OakDBH +offset(log(Deploy.Duration)), data = deerDataSum, family = poisson)
summary(glm.po.fullS)

#test for overdispersion indicates we should be using negative binomial 
sum.po.fullS <- summary(glm.po.fullS)
phi.po.fullS <- sum.po.fullS$deviance/sum.po.fullS$df.residual 

#Exploratory Poisson model for Fall data
glm.po.fullF <- glm(nSeqs ~ Log.in.View + Summer.Fall.EDD + Height_cm + log10(Num_Stems) + OakDBH +offset(log(Deploy.Duration)), data = deerDataFall, family = poisson)
summary(glm.po.fullF)

#test for overdispersion indicates we should be using negative binomial 
sum.po.fullF <- summary(glm.po.fullF)
phi.po.fullF <- sum.po.fullF$deviance/sum.po.fullF$df.residual 

#Exploratory Poisson model for Winter data
glm.po.fullW <- glm(nSeqs ~ Log.in.View + Winter.Spring.EDD + Height_cm + log10(Num_Stems) + OakDBH +offset(log(Deploy.Duration)), data = deerDataWin, family = poisson)
summary(glm.po.fullW)

#test for overdispersion indicates we should be using negative binomial 
sum.po.fullW <- summary(glm.po.fullW)
phi.po.fullW <- sum.po.fullW$deviance/sum.po.fullW$df.residual 

#Exploratory Poisson model for Spring data
glm.po.fullSp <- glm(nSeqs ~ Log.in.View + Winter.Spring.EDD + Height_cm + log10(Num_Stems) + OakDBH +offset(log(Deploy.Duration)), data = deerDataSpr, family = poisson)
summary(glm.po.fullSp)

#test for overdispersion indicates we should be using negative binomial 
sum.po.fullSp <- summary(glm.po.fullSp)
phi.po.fullSp <- sum.po.fullSp$deviance/sum.po.fullSp$df.residual 


# Negative Binomial Regression - WTD --------------------------------------
#Summer

glm.nb.fullS <-
  glm.nb(
    nSeqs ~ Log.in.View + Summer.Fall.EDD + Height_cm + log10(Num_Stems) + OakDBH +
      offset(log(Deploy.Duration)),
    data = deerDataSum
  )
summary(glm.nb.fullS)

glm.nb.SI <- glm.nb(nSeqs ~ 1 + offset(log(Deploy.Duration)), data = deerDataSum)

glm.nb.Slog <-
  glm.nb(nSeqs ~ Log.in.View + offset(log(Deploy.Duration)), data = deerDataSum)

glm.nb.Shgt <-
  glm.nb(nSeqs ~ Height_cm + offset(log(Deploy.Duration)), data = deerDataSum)

glm.nb.Soak <-
  glm.nb(nSeqs ~ OakDBH + offset(log(Deploy.Duration)), data = deerDataSum)

glm.nb.Sstems <-
  glm.nb(nSeqs ~ Num_Stems + offset(log(Deploy.Duration)), data = deerDataSum)

glm.nb.Sedd <-
  glm.nb(nSeqs ~ Summer.Fall.EDD + offset(log(Deploy.Duration)), data = deerDataSum)

glm.nb.Slog_stems <-
  glm.nb(nSeqs ~ Log.in.View + Num_Stems + offset(log(Deploy.Duration)), data = deerDataSum)

mod.namesDS <- c("Full", "log", "height", "oaks", "stems", "EDD", "Log + Stems", "Intercept")
glmNBDS <- cbind(mod.namesDS, c(AICc(glm.nb.fullS), AICc(glm.nb.Slog), AICc(glm.nb.Shgt), AICc(glm.nb.Soak), AICc(glm.nb.Sstems), AICc(glm.nb.Sedd), AICc(glm.nb.Slog_stems), AICc(glm.nb.SI)))

glmNBDS

summary(glm.nb.Slog) #Likely to be best model, deer detections in Summer much less likely when a log is in view

#Fall
glm.nb.fullF <- glm.nb(nSeqs ~ Log.in.View + Summer.Fall.EDD + Height_cm + log10(Num_Stems) + OakDBH +offset(log(Deploy.Duration)), data = deerDataFall)
summary(glm.nb.fullF)

glm.nb.FI <- glm.nb(nSeqs ~ 1 + offset(log(Deploy.Duration)), data = deerDataFall)

glm.nb.Flog <-
  glm.nb(nSeqs ~ Log.in.View + offset(log(Deploy.Duration)), data = deerDataFall)

glm.nb.Fhgt <-
  glm.nb(nSeqs ~ Height_cm + offset(log(Deploy.Duration)), data = deerDataFall)

glm.nb.FhgtLog <-
  glm.nb(nSeqs ~ Height_cm + Log.in.View + offset(log(Deploy.Duration)), data = deerDataFall)

glm.nb.Foak <-
  glm.nb(nSeqs ~ OakDBH + offset(log(Deploy.Duration)), data = deerDataFall)

glm.nb.Fstems <-
  glm.nb(nSeqs ~ Num_Stems + offset(log(Deploy.Duration)), data = deerDataFall)

glm.nb.Fedd <-
  glm.nb(nSeqs ~ Summer.Fall.EDD + offset(log(Deploy.Duration)), data = deerDataFall)

mod.namesDF <- c("Full", "log", "height", "log + height", "oaks", "stems", "EDD", "Intercept")
glmNBDF <- cbind(mod.namesDF, c(AICc(glm.nb.fullF), AICc(glm.nb.Flog), AICc(glm.nb.Fhgt), AICc(glm.nb.FhgtLog), AICc(glm.nb.Foak), AICc(glm.nb.Fstems), AICc(glm.nb.Fedd), AICc(glm.nb.FI)))
glmNBDF  #No real significant parameters it looks like in Fall for deer.

#Winter
glm.nb.fullW <- glm.nb(nSeqs ~ Log.in.View + Winter.Spring.EDD + Height_cm + log10(Num_Stems) + OakDBH +offset(log(Deploy.Duration)), data = deerDataWin)
summary(glm.nb.fullW)

glm.nb.WI <- glm.nb(nSeqs ~ 1 + offset(log(Deploy.Duration)), data = deerDataWin)

glm.nb.Wlog <-
  glm.nb(nSeqs ~ Log.in.View + offset(log(Deploy.Duration)), data = deerDataWin)

glm.nb.Whgt <-
  glm.nb(nSeqs ~ Height_cm + offset(log(Deploy.Duration)), data = deerDataWin)

glm.nb.WhgtLog <-
  glm.nb(nSeqs ~ Height_cm + Log.in.View + offset(log(Deploy.Duration)), data = deerDataWin)

glm.nb.Woak <-
  glm.nb(nSeqs ~ OakDBH + offset(log(Deploy.Duration)), data = deerDataWin)

glm.nb.Wstems <-
  glm.nb(nSeqs ~ Num_Stems + offset(log(Deploy.Duration)), data = deerDataWin)

glm.nb.Wedd <-
  glm.nb(nSeqs ~ Winter.Spring.EDD + offset(log(Deploy.Duration)), data = deerDataWin)

mod.namesDW <- c("Full", "log", "height", "log + height", "oaks", "stems", "EDD", "Intercept")
glmNBDW <- cbind(mod.namesDW, c(AICc(glm.nb.fullW), AICc(glm.nb.Wlog), AICc(glm.nb.Whgt), AICc(glm.nb.WhgtLog), AICc(glm.nb.Woak), AICc(glm.nb.Wstems), AICc(glm.nb.Wedd), AICc(glm.nb.WI)))
glmNBDW #Nothing too interesting in winter


#Spring
glm.nb.fullSp <- glm.nb(nSeqs ~ Log.in.View + Winter.Spring.EDD + Height_cm + log10(Num_Stems) + OakDBH +offset(log(Deploy.Duration)), data = deerDataSpr)
summary(glm.nb.fullSp)

glm.nb.SpI <- glm.nb(nSeqs ~ 1 + offset(log(Deploy.Duration)), data = deerDataSpr)

glm.nb.Splog <-
  glm.nb(nSeqs ~ Log.in.View + offset(log(Deploy.Duration)), data = deerDataSpr)

glm.nb.Sphgt <-
  glm.nb(nSeqs ~ Height_cm + offset(log(Deploy.Duration)), data = deerDataSpr)

glm.nb.SphgtLog <-
  glm.nb(nSeqs ~ Height_cm + Log.in.View + offset(log(Deploy.Duration)), data = deerDataSpr)

glm.nb.Spoak <-
  glm.nb(nSeqs ~ OakDBH + offset(log(Deploy.Duration)), data = deerDataSpr)

glm.nb.Spstems <-
  glm.nb(nSeqs ~ Num_Stems + offset(log(Deploy.Duration)), data = deerDataSpr)

glm.nb.Spedd <-
  glm.nb(nSeqs ~ Winter.Spring.EDD + offset(log(Deploy.Duration)), data = deerDataSpr)

mod.namesDSp <- c("Full", "log", "height", "log + height", "oaks", "stems", "EDD", "Intercept")
glmNBDSp <- cbind(mod.namesDSp, c(AICc(glm.nb.fullSp), AICc(glm.nb.Splog), AICc(glm.nb.Sphgt), AICc(glm.nb.SphgtLog), AICc(glm.nb.Spoak), AICc(glm.nb.Spstems), AICc(glm.nb.Spedd), AICc(glm.nb.SpI)))
glmNBDSp #Nothing too interesting in Spring either...

#____________________________________
#Response Curves for our Covariates
#____________________________________
# NEED TO EDIT TEXT BELOW THIS LINE
meanEDD <- mean(deerDataSum$EDD)
medoak <- median(deerDataSum$Num_Oaks)
medstems <- median(deerDataSum$Num_Stems)

#Response Curves for Number of Stems. Will use average EDD, and median # of oaks. We will draw graphs for logs present and for logs absent

seq.stems <- seq(min(deerDataSum$Num_Stems), max(deerDataSum$Num_Stems), length.out = 100)
nd.seq.stem.log <- data.frame(Cam_Nights = 61, EDD = meanEDD, Num_Oaks = medoak, Log = "YES", Num_Stems = seq.stems)
pred.stems.l <- predict(glm.po2, newdata = nd.seq.stem.log, se=TRUE, type = "link")
plot(Number_of_Sequences ~ Num_Stems, data=deerDataSum, pch=16, col="orange", las=1, log = "x", xlab = "log(Num_Stems)", ylab = "# Deer of Sequences")
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
seq.EDD <- seq(min(deerDataSum$EDD), max(deerDataSum$EDD), length.out = 100)
nd.seq.EDD.log <- data.frame(Cam_Nights = 61, EDD = seq.EDD, Num_Oaks = medoak, Log = "YES", Num_Stems = medstems)
pred.EDD.l <- predict(glm.po2, newdata = nd.seq.EDD.log, se=TRUE, type = "link")
plot(Number_of_Sequences ~ EDD, data=deerDataSum, pch=16, col="orange", las=1, xlab = "Effective Detection Distance (m)", ylab = "# of Deer Sequences")
lines(nd.seq.EDD.log$EDD, exp(pred.EDD.l$fit), col="purple")
lines(nd.seq.EDD.log$EDD, exp(pred.EDD.l$fit + 1.96*pred.EDD.l$se.fit), lty=2, col="purple")
lines(nd.seq.EDD.log$EDD, exp(pred.EDD.l$fit - 1.96*pred.EDD.l$se.fit), lty=2, col="purple")

nd.seq.EDD.nolog <- data.frame(Cam_Nights = 61, EDD = seq.EDD, Num_Oaks = medoak, Log = "NO", Num_Stems = medstems)
pred.EDD.nl <- predict(glm.po2, newdata = nd.seq.EDD.nolog, se=TRUE, type = "link")
lines(nd.seq.EDD.nolog$EDD, exp(pred.EDD.nl$fit), col="green")
lines(nd.seq.EDD.nolog$EDD, exp(pred.EDD.nl$fit + 1.96*pred.EDD.nl$se.fit), lty=2, col="green")
lines(nd.seq.EDD.nolog$EDD, exp(pred.EDD.nl$fit - 1.96*pred.EDD.nl$se.fit), lty=2, col="green")
