###Regression Analysis of SCBI Camera Grid Data

require(AICcmodavg)
require(MASS) #for glm.nb function
source("scripts/pairsPannelFunctions.r")


# Data Import and Exploration - WTD ---------------------------------------

load("data/sqDataSum.RData")
load("data/sqDataFall.RData")
load("data/sqDataWin.RData")
load("data/sqDataSpr.RData")

pairs(sqDataSum[,c(4,7,8,11,12,13,15)], diag.panel = panel.hist, lower.panel = panel.smooth, upper.panel = panel.cor)
str(sqDataSum)

hist(sqDataSum$Num_Stems)
hist(log10(sqDataSum$Num_Stems)) #because one site has so many stems, we need to log10 transform this parameter for all models
hist(sqDataSum$Squirrel_EDD_4S)

#Univariate exploratory plots for each of our 5 covariates
par.default <- par(no.readonly = T)
par(mfrow = c(2,2))
plot(nSeqs ~ Squirrel_EDD_4S, data = sqDataSum, main = "Summer")
plot(nSeqs ~ Squirrel_EDD_4S, data = sqDataFall, main = "Fall")
plot(nSeqs ~ Squirrel_EDD_WSp, data = sqDataWin, main = "Winter")
plot(nSeqs ~ Squirrel_EDD_WSp, data = sqDataSpr, main = "Spring")

plot(nSeqs ~ OakDBH, data = sqDataSum, main = "Summer")
plot(nSeqs ~ OakDBH, data = sqDataFall, main = "Fall")
plot(nSeqs ~ OakDBH, data = sqDataWin, main = "Winter")
plot(nSeqs ~ OakDBH, data = sqDataSpr, main = "Spring")

plot(nSeqs ~ log10(Num_Stems), data = sqDataSum, main = "Summer")
plot(nSeqs ~ log10(Num_Stems), data = sqDataFall, main = "Fall")
plot(nSeqs ~ log10(Num_Stems), data = sqDataWin, main = "Winter")
plot(nSeqs ~ log10(Num_Stems), data = sqDataSpr, main = "Spring")

#These seems to show a quadratic relationship actually, which might make sense
plot(nSeqs ~ Height_cm, data = sqDataSum, main = "Summer")
plot(nSeqs ~ Height_cm, data = sqDataFall, main = "Fall")
plot(nSeqs ~ Height_cm, data = sqDataWin, main = "Winter")
plot(nSeqs ~ Height_cm, data = sqDataSpr, main = "Spring")

plot(nSeqs ~ Log.in.View, data = sqDataSum, main = "Summer")
plot(nSeqs ~ Log.in.View, data = sqDataFall, main = "Fall")
plot(nSeqs ~ Log.in.View, data = sqDataWin, main = "Winter")
plot(nSeqs ~ Log.in.View, data = sqDataSpr, main = "Spring")
par(par.default)

hist(sqDataSum$Deploy.Duration, breaks = 5)


# Full Poisson GLM and Overdispersion -------------------------------------

#Exploratory Poisson model for Summer data to see if a negative binomial model is needed due to potential overdispersion
glm.po.fullS <- glm(nSeqs ~ Log.in.View + Squirrel_EDD_4S + Height_cm + log10(Num_Stems) + OakDBH +offset(log(Deploy.Duration)), data = sqDataSum, family = poisson)
summary(glm.po.fullS)

#test for overdispersion indicates we should be using negative binomial 
sum.po.fullS <- summary(glm.po.fullS)
phi.po.fullS <- sum.po.fullS$deviance/sum.po.fullS$df.residual 

#Exploratory Poisson model for Fall data
glm.po.fullF <- glm(nSeqs ~ Log.in.View + Squirrel_EDD_4S + Height_cm + log10(Num_Stems) + OakDBH +offset(log(Deploy.Duration)), data = sqDataFall, family = poisson)
summary(glm.po.fullF)

#test for overdispersion indicates we should be using negative binomial 
sum.po.fullF <- summary(glm.po.fullF)
phi.po.fullF <- sum.po.fullF$deviance/sum.po.fullF$df.residual 

#Exploratory Poisson model for Winter data
glm.po.fullW <- glm(nSeqs ~ Log.in.View + Squirrel_EDD_WSp + Height_cm + log10(Num_Stems) + OakDBH +offset(log(Deploy.Duration)), data = sqDataWin, family = poisson)
summary(glm.po.fullW)

#test for overdispersion indicates we should be using negative binomial 
sum.po.fullW <- summary(glm.po.fullW)
phi.po.fullW <- sum.po.fullW$deviance/sum.po.fullW$df.residual 

#Exploratory Poisson model for Spring data
glm.po.fullSp <- glm(nSeqs ~ Log.in.View + Squirrel_EDD_WSp + Height_cm + log10(Num_Stems) + OakDBH +offset(log(Deploy.Duration)), data = sqDataSpr, family = poisson)
summary(glm.po.fullSp)

#test for overdispersion indicates we should be using negative binomial 
sum.po.fullSp <- summary(glm.po.fullSp)
phi.po.fullSp <- sum.po.fullSp$deviance/sum.po.fullSp$df.residual 


# Negative Binomial Regression - Squirrels --------------------------------------
#
#Summer - next try log + EDD + height (and maybe quadratic on height?)

glm.nb.fullS <-
  glm.nb(
    nSeqs ~ Log.in.View + Squirrel_EDD_4S + Height_cm + log10(Num_Stems) + OakDBH +
      offset(log(Deploy.Duration)),
    data = sqDataSum
  )
summary(glm.nb.fullS)

glm.nb.SI <- glm.nb(nSeqs ~ 1 + offset(log(Deploy.Duration)), data = sqDataSum)

glm.nb.Slog <-
  glm.nb(nSeqs ~ Log.in.View + offset(log(Deploy.Duration)), data = sqDataSum)

glm.nb.Shgt <-
  glm.nb(nSeqs ~ Height_cm + offset(log(Deploy.Duration)), data = sqDataSum)

glm.nb.Soak <-
  glm.nb(nSeqs ~ OakDBH + offset(log(Deploy.Duration)), data = sqDataSum)

glm.nb.Sstems <-
  glm.nb(nSeqs ~ Num_Stems + offset(log(Deploy.Duration)), data = sqDataSum)

glm.nb.Sedd <-
  glm.nb(nSeqs ~ Squirrel_EDD_4S + offset(log(Deploy.Duration)), data = sqDataSum)

mod.namesDS <- c("Full", "log", "height", "oaks", "stems", "EDD", "Intercept")
glmNBDS <- cbind(mod.namesDS, c(AICc(glm.nb.fullS), AICc(glm.nb.Slog), AICc(glm.nb.Shgt), AICc(glm.nb.Soak), AICc(glm.nb.Sstems), AICc(glm.nb.Sedd), AICc(glm.nb.SI)))
glmNBDS



#Fall - next try log + EDD
glm.nb.fullF <- glm.nb(nSeqs ~ Log.in.View + Squirrel_EDD_4S + Height_cm + log10(Num_Stems) + OakDBH +offset(log(Deploy.Duration)), data = sqDataFall)
summary(glm.nb.fullF)

glm.nb.FI <- glm.nb(nSeqs ~ 1 + offset(log(Deploy.Duration)), data = sqDataFall)

glm.nb.Flog <-
  glm.nb(nSeqs ~ Log.in.View + offset(log(Deploy.Duration)), data = sqDataFall)

glm.nb.Fhgt <-
  glm.nb(nSeqs ~ Height_cm + offset(log(Deploy.Duration)), data = sqDataFall)

glm.nb.FhgtLog <-
  glm.nb(nSeqs ~ Height_cm + Log.in.View + offset(log(Deploy.Duration)), data = sqDataFall)

glm.nb.Foak <-
  glm.nb(nSeqs ~ OakDBH + offset(log(Deploy.Duration)), data = sqDataFall)

glm.nb.Fstems <-
  glm.nb(nSeqs ~ Num_Stems + offset(log(Deploy.Duration)), data = sqDataFall)

glm.nb.Fedd <-
  glm.nb(nSeqs ~ Squirrel_EDD_4S + offset(log(Deploy.Duration)), data = sqDataFall)

mod.namesDF <- c("Full", "log", "height", "log + height", "oaks", "stems", "EDD", "Intercept")
glmNBDF <- cbind(mod.namesDF, c(AICc(glm.nb.fullF), AICc(glm.nb.Flog), AICc(glm.nb.Fhgt), AICc(glm.nb.FhgtLog), AICc(glm.nb.Foak), AICc(glm.nb.Fstems), AICc(glm.nb.Fedd), AICc(glm.nb.FI)))
glmNBDF

#Winter - Next try log + stems
glm.nb.fullW <- glm.nb(nSeqs ~ Log.in.View + Squirrel_EDD_WSp + Height_cm + log10(Num_Stems) + OakDBH +offset(log(Deploy.Duration)), data = sqDataWin)
summary(glm.nb.fullW)

glm.nb.WI <- glm.nb(nSeqs ~ 1 + offset(log(Deploy.Duration)), data = sqDataWin)

glm.nb.Wlog <-
  glm.nb(nSeqs ~ Log.in.View + offset(log(Deploy.Duration)), data = sqDataWin)

glm.nb.Whgt <-
  glm.nb(nSeqs ~ Height_cm + offset(log(Deploy.Duration)), data = sqDataWin)

glm.nb.WhgtLog <-
  glm.nb(nSeqs ~ Height_cm + Log.in.View + offset(log(Deploy.Duration)), data = sqDataWin)

glm.nb.Woak <-
  glm.nb(nSeqs ~ OakDBH + offset(log(Deploy.Duration)), data = sqDataWin)

glm.nb.Wstems <-
  glm.nb(nSeqs ~ Num_Stems + offset(log(Deploy.Duration)), data = sqDataWin)

glm.nb.Wedd <-
  glm.nb(nSeqs ~ Squirrel_EDD_WSp + offset(log(Deploy.Duration)), data = sqDataWin)

mod.namesDW <- c("Full", "log", "height", "log + height", "oaks", "stems", "EDD", "Intercept")
glmNBDW <- cbind(mod.namesDW, c(AICc(glm.nb.fullW), AICc(glm.nb.Wlog), AICc(glm.nb.Whgt), AICc(glm.nb.WhgtLog), AICc(glm.nb.Woak), AICc(glm.nb.Wstems), AICc(glm.nb.Wedd), AICc(glm.nb.WI)))
glmNBDW


#Spring
glm.nb.fullSp <- glm.nb(nSeqs ~ Log.in.View + Squirrel_EDD_WSp + Height_cm + log10(Num_Stems) + OakDBH +offset(log(Deploy.Duration)), data = sqDataSpr)
summary(glm.nb.fullSp)

glm.nb.SpI <- glm.nb(nSeqs ~ 1 + offset(log(Deploy.Duration)), data = sqDataSpr)

glm.nb.Splog <-
  glm.nb(nSeqs ~ Log.in.View + offset(log(Deploy.Duration)), data = sqDataSpr)

glm.nb.Sphgt <-
  glm.nb(nSeqs ~ Height_cm + offset(log(Deploy.Duration)), data = sqDataSpr)

glm.nb.SphgtLog <-
  glm.nb(nSeqs ~ Height_cm + Log.in.View + offset(log(Deploy.Duration)), data = sqDataSpr)

glm.nb.Spoak <-
  glm.nb(nSeqs ~ OakDBH + offset(log(Deploy.Duration)), data = sqDataSpr)

glm.nb.Spstems <-
  glm.nb(nSeqs ~ Num_Stems + offset(log(Deploy.Duration)), data = sqDataSpr)

glm.nb.Spedd <-
  glm.nb(nSeqs ~ Squirrel_EDD_WSp + offset(log(Deploy.Duration)), data = sqDataSpr)

mod.namesDSp <- c("Full", "log", "height", "log + height", "oaks", "stems", "EDD", "Intercept")
glmNBDSp <- cbind(mod.namesDSp, c(AICc(glm.nb.fullSp), AICc(glm.nb.Splog), AICc(glm.nb.Sphgt), AICc(glm.nb.SphgtLog), AICc(glm.nb.Spoak), AICc(glm.nb.Spstems), AICc(glm.nb.Spedd), AICc(glm.nb.SpI)))
glmNBDSp