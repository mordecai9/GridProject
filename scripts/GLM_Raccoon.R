###Regression Analysis of SCBI Camera Grid Data for Raccoon

require(AICcmodavg)
require(MASS) #for glm.nb function
source("scripts/pairsPannelFunctions.r")


# Data Import and Exploration - Raccoon ---------------------------------------
# We won't be using winter data here since bear really weren't detected in winter
load("data/racDataSum.RData")
load("data/racDataFall.RData")
load("data/racDataWin.RData")
load("data/racDataSpr.RData")

pairs(bearDataSum[,c(4,7,13,14,15,17)], diag.panel = panel.hist, lower.panel = panel.smooth, upper.panel = panel.cor)
str(bearDataSum)

#Looking at how many cameras picked up rac in each season. Looks like enough in each, though a bit low in summer
barplot(racDataSum$nSeqs)
barplot(racDataFall$nSeqs)
barplot(racDataWin$nSeqs)
barplot(racDataSpr$nSeqs)

hist(racDataSum$Num_Stems)
hist(log10(racDataSum$Num_Stems)) #because one site has so many stems, we need to log10 transform this parameter for all models
hist(racDataSum$Raccoon.EDD)


#Univariate exploratory plots for each of our 5 covariates
par.default <- par(no.readonly = T)
par(mfrow = c(2,2))
plot(nSeqs ~ Raccoon.EDD, data = racDataSum, main = "Summer")
plot(nSeqs ~ Raccoon.EDD, data = racDataFall, main = "Fall")
plot(nSeqs ~ Raccoon.EDD, data = racDataWin, main = "Winter")
plot(nSeqs ~ Raccoon.EDD, data = racDataSpr, main = "Spring")

plot(nSeqs ~ OakDBH, data = racDataSum, main = "Summer")
plot(nSeqs ~ OakDBH, data = racDataFall, main = "Fall")
plot(nSeqs ~ OakDBH, data = racDataWin, main = "Winter")
plot(nSeqs ~ OakDBH, data = racDataSpr, main = "Spring")

plot(nSeqs ~ log10(Num_Stems), data = racDataSum, main = "Summer")
plot(nSeqs ~ log10(Num_Stems), data = racDataFall, main = "Fall")
plot(nSeqs ~ log10(Num_Stems), data = racDataWin, main = "Winter")
plot(nSeqs ~ log10(Num_Stems), data = racDataSpr, main = "Spring")

plot(nSeqs ~ Height_cm, data = racDataSum, main = "Summer")
plot(nSeqs ~ Height_cm, data = racDataFall, main = "Fall")
plot(nSeqs ~ Height_cm, data = racDataWin, main = "Winter")
plot(nSeqs ~ Height_cm, data = racDataSpr, main = "Spring")

plot(nSeqs ~ Log.in.View, data = racDataSum, main = "Summer")
plot(nSeqs ~ Log.in.View, data = racDataFall, main = "Fall")
plot(nSeqs ~ Log.in.View, data = racDataWin, main = "Winter")
plot(nSeqs ~ Log.in.View, data = racDataSpr, main = "Spring")
par(par.default)

hist(racDataSum$Deploy.Duration, breaks = 5)


# Full Poisson GLM and Overdispersion -------------------------------------

#Exploratory Poisson model for Summer data to see if a negative binomial model is needed due to potential overdispersion
glm.po.fullS <- glm(nSeqs ~ Log.in.View + Raccoon.EDD + Height_cm + log10(Num_Stems) + OakDBH +offset(log(Deploy.Duration)), data = racDataSum, family = poisson)
summary(glm.po.fullS)

#test for overdispersion is 2.25
sum.po.fullS <- summary(glm.po.fullS)
phi.po.fullS <- sum.po.fullS$deviance/sum.po.fullS$df.residual 

#Exploratory Poisson model for Fall data
glm.po.fullF <- glm(nSeqs ~ Log.in.View + Raccoon.EDD + Height_cm + log10(Num_Stems) + OakDBH +offset(log(Deploy.Duration)), data = racDataFall, family = poisson)
summary(glm.po.fullF)

#test for overdispersion (4.12) indicates we should be using negative binomial 
sum.po.fullF <- summary(glm.po.fullF)
phi.po.fullF <- sum.po.fullF$deviance/sum.po.fullF$df.residual 

#Exploratory Poisson model for Winter data
glm.po.fullW <- glm(nSeqs ~ Log.in.View + Raccoon.EDD + Height_cm + log10(Num_Stems) + OakDBH +offset(log(Deploy.Duration)), data = racDataWin, family = poisson)
summary(glm.po.fullW)

#test for overdispersion indicates we should be using negative binomial 
sum.po.fullW <- summary(glm.po.fullW)
phi.po.fullW <- sum.po.fullW$deviance/sum.po.fullW$df.residual 

#Exploratory Poisson model for Spring data
glm.po.fullSp <- glm(nSeqs ~ Log.in.View + Raccoon.EDD + Height_cm + log10(Num_Stems) + OakDBH +offset(log(Deploy.Duration)), data = racDataSpr, family = poisson)
summary(glm.po.fullSp)

#test for overdispersion is 2.66 
sum.po.fullSp <- summary(glm.po.fullSp)
phi.po.fullSp <- sum.po.fullSp$deviance/sum.po.fullSp$df.residual 


# Negative Binomial Regression - Raccoon --------------------------------------
#Summer
#START FIXING CODE HERE

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

mod.namesDS <- c("Full", "log", "height", "oaks", "stems", "EDD", "Intercept")
glmNBDS <- cbind(mod.namesDS, c(AICc(glm.nb.fullS), AICc(glm.nb.Slog), AICc(glm.nb.Shgt), AICc(glm.nb.Soak), AICc(glm.nb.Sstems), AICc(glm.nb.Sedd), AICc(glm.nb.SI)))
glmNBDS

summary(glm.nb.log)

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
glmNBDF

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
glmNBDW


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
glmNBDSp