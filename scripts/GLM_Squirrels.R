###Regression Analysis of SCBI Camera Grid Data

rm(list = ls())

require(AICcmodavg)
require(MASS) #for glm.nb function
require(tidyverse)
source("scripts/pairsPannelFunctions.r")

# Data Import and Exploration - All Squirrels ---------------------------------------

load("data/sqDataSum.RData")
load("data/sqDataFall.RData")
load("data/sqDataWin.RData")
load("data/sqDataSpr.RData")

pairs(sqDataSum[,c(4,9,10,12,13,15,17)], diag.panel = panel.hist, lower.panel = panel.smooth, upper.panel = panel.cor)
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


# Negative Binomial Regression - Squirrels in Summer --------------------------------

#First test to see if interaction between log and EDD improves that two variable model

glm.nb.Slog_Edd <-
  glm.nb(nSeqs ~ Log.in.View + Squirrel_EDD_4S + offset(log(Deploy.Duration)), data = sqDataSum)
AICc(glm.nb.Slog_Edd)

glm.nb.SlogXEdd <-
  glm.nb(nSeqs ~ Log.in.View * Squirrel_EDD_4S + offset(log(Deploy.Duration)), data = sqDataSum)
AICc(glm.nb.SlogXEdd) #No support for interaction 


#Compare height alone to height quadratic
glm.nb.SHgt2 <-
  glm.nb(
    nSeqs ~ poly(Height_cm, 2, raw = T) +
      offset(log(Deploy.Duration)),
    data = sqDataSum
  )
AICc(glm.nb.SHgt2)
AICc(glm.nb.Shgt) #no support for quadratic effect

#Full Model Selection for Squirrels in Summer

glm.nb.fullS <-
  glm.nb(
    nSeqs ~ Log.in.View + Squirrel_EDD_4S + Height_cm + log10(Num_Stems) + OakDBH +
      offset(log(Deploy.Duration)),
    data = sqDataSum
  )


ModListSumSq <- list(NA)


ModListSumSq[[1]] <- glm.nb.SI <- glm.nb(nSeqs ~ 1 + offset(log(Deploy.Duration)), data = sqDataSum)

ModListSumSq[[2]] <- glm.nb.Slog <-
  glm.nb(nSeqs ~ Log.in.View + offset(log(Deploy.Duration)), data = sqDataSum)

ModListSumSq[[3]] <- glm.nb.Shgt <-
  glm.nb(nSeqs ~ Height_cm + offset(log(Deploy.Duration)), data = sqDataSum)

ModListSumSq[[4]] <- glm.nb.Soak <-
  glm.nb(nSeqs ~ OakDBH + offset(log(Deploy.Duration)), data = sqDataSum)

ModListSumSq[[5]] <- glm.nb.Sstems <-
  glm.nb(nSeqs ~ log10(Num_Stems) + offset(log(Deploy.Duration)), data = sqDataSum)

ModListSumSq[[6]] <- glm.nb.Sedd <-
  glm.nb(nSeqs ~ Squirrel_EDD_4S + offset(log(Deploy.Duration)), data = sqDataSum)

ModListSumSq[[7]] <- glm.nb.Slog_stems <-
  glm.nb(nSeqs ~ Log.in.View + log10(Num_Stems) + offset(log(Deploy.Duration)), data = sqDataSum)

ModListSumSq[[8]] <-glm.nb.Slog_Edd <-
  glm.nb(nSeqs ~ Log.in.View + Squirrel_EDD_4S + offset(log(Deploy.Duration)), data = sqDataSum)

ModListSumSq[[9]] <- glm.nb.Slog_Oak <-
  glm.nb(nSeqs ~ Log.in.View + OakDBH + offset(log(Deploy.Duration)), data = sqDataSum)

ModListSumSq[[10]] <- glm.nb.Slog_Hgt <-
  glm.nb(nSeqs ~ Log.in.View + Height_cm + offset(log(Deploy.Duration)), data = sqDataSum)

ModListSumSq[[11]] <- glm.nb.SEdd_Hgt <-
  glm.nb(nSeqs ~ Squirrel_EDD_4S + Height_cm + offset(log(Deploy.Duration)), data = sqDataSum)

ModListSumSq[[12]] <-  glm.nb.SEdd_stems <-
  glm.nb(nSeqs ~ Squirrel_EDD_4S + log10(Num_Stems) + offset(log(Deploy.Duration)), data = sqDataSum)

ModListSumSq[[13]] <- glm.nb.SEdd_oak <-
  glm.nb(nSeqs ~ Squirrel_EDD_4S + OakDBH + offset(log(Deploy.Duration)), data = sqDataSum)

ModListSumSq[[14]] <- glm.nb.Sstems_hgt <-
  glm.nb(nSeqs ~ Height_cm + log10(Num_Stems) + offset(log(Deploy.Duration)), data = sqDataSum)

ModListSumSq[[15]] <- glm.nb.Sstems_oak <-
  glm.nb(nSeqs ~ OakDBH + log10(Num_Stems) + offset(log(Deploy.Duration)), data = sqDataSum)

ModListSumSq[[16]] <- glm.nb.Shgt_oak <-
  glm.nb(nSeqs ~ Height_cm + OakDBH + offset(log(Deploy.Duration)), data = sqDataSum)

ModListSumSq[[17]] <- glm.nb.SlogStemsHgt <-
  glm.nb(
    nSeqs ~ Log.in.View + Height_cm + log10(Num_Stems) +
      offset(log(Deploy.Duration)),
    data = sqDataSum
  )

ModListSumSq[[18]] <- glm.nb.SlogStemsOak <-
  glm.nb(
    nSeqs ~ Log.in.View + OakDBH + log10(Num_Stems) +
      offset(log(Deploy.Duration)),
    data = sqDataSum
  )

ModListSumSq[[19]] <- glm.nb.SlogStemsEdd <-
  glm.nb(
    nSeqs ~ Log.in.View + Squirrel_EDD_4S + log10(Num_Stems) +
      offset(log(Deploy.Duration)),
    data = sqDataSum
  )

ModListSumSq[[20]] <- glm.nb.SHgtStemsEdd <-
  glm.nb(
    nSeqs ~ Height_cm + Squirrel_EDD_4S + log10(Num_Stems) +
      offset(log(Deploy.Duration)),
    data = sqDataSum
  )

ModListSumSq[[21]] <- glm.nb.SHgtStemsOak <-
  glm.nb(
    nSeqs ~ Height_cm + OakDBH + log10(Num_Stems) +
      offset(log(Deploy.Duration)),
    data = sqDataSum
  )

ModListSumSq[[22]] <- glm.nb.SHgtEddOak <-
  glm.nb(
    nSeqs ~ Height_cm + OakDBH + Squirrel_EDD_4S +
      offset(log(Deploy.Duration)),
    data = sqDataSum
  )
ModListSumSq[[23]] <- glm.nb.SHgtEddLog <-
  glm.nb(
    nSeqs ~ Height_cm + Log.in.View + Squirrel_EDD_4S +
      offset(log(Deploy.Duration)),
    data = sqDataSum
  )
ModListSumSq[[24]] <- glm.nb.SHgtOakLog <-
  glm.nb(
    nSeqs ~ Height_cm + Log.in.View + OakDBH +
      offset(log(Deploy.Duration)),
    data = sqDataSum
  )

ModListSumSq[[25]] <- glm.nb.SOakStemsEdd <-
  glm.nb(
    nSeqs ~ log10(Num_Stems) + Squirrel_EDD_4S + OakDBH +
      offset(log(Deploy.Duration)),
    data = sqDataSum
  )

ModListSumSq[[26]] <- glm.nb.SLogOakEdd <-
  glm.nb(
    nSeqs ~ Log.in.View + Squirrel_EDD_4S + OakDBH +
      offset(log(Deploy.Duration)),
    data = sqDataSum
  )
modnamesS <- c("Intercept","Log","Height","Oak","Stems","EDD","Log + Stems","Log + EDD","Log + Oak","Log + Height","EDD + Height","EDD + Stems","EDD + Oak","Stems + Height","Stems + Oak","Height + Oak","Log + Stems + Height","Log + Stems + Oak","Log + Stems + EDD","Height + Stems + EDD","Height + Stems + Oak","Height + EDD + Oak","Height + EDD + Log","Height + Oak + Log", "Stems + EDD + Oak", "Log + EDD + Oak")
modtabSumSq <- aictab(cand.set = ModListSumSq, modnames = modnamesS)
modtabSumSq

#Best model has height, EDD and log
summary(glm.nb.SHgtEddLog)
AICc(glm.nb.SHgtEddLog)

#Summed Model Weights for Squirrels in Summer from full table. 
SumSW_edd <- sum(modtabSumSq$AICcWt[grep("EDD", modtabSumSq$Modnames)])
SumSW_OakDBH <- sum(modtabSumSq$AICcWt[grep("Oak", modtabSumSq$Modnames)])
SumSW_Stems <- sum(modtabSumSq$AICcWt[grep("Stems", modtabSumSq$Modnames)])
SumSW_Log <- sum(modtabSumSq$AICcWt[grep("Log", modtabSumSq$Modnames)])
SumSW_Height <- sum(modtabSumSq$AICcWt[grep("Height", modtabSumSq$Modnames)])

#See if this model is better with EDD and log interaction, just out of curiosity
glm.nb.SHgtEddXLog <-
  glm.nb(
    nSeqs ~ Height_cm + Log.in.View * Squirrel_EDD_4S +
      offset(log(Deploy.Duration)),
    data = sqDataSum
  )
AICc(glm.nb.SHgtEddXLog) #not worth having the interaction here
summary(glm.nb.SHgtEddXLog)


#Compare best model with height as quadratic in same model, just out of curiosity
glm.nb.SHgt2EddLog <-
  glm.nb(
    nSeqs ~ poly(Height_cm, 2, raw = T) + Log.in.View + Squirrel_EDD_4S +
      offset(log(Deploy.Duration)),
    data = sqDataSum
  )
AICc(glm.nb.SHgt2EddLog) #no support for quadratic effect

#Explained deviance of best model
1 - glm.nb.SHgtEddLog$deviance / glm.nb.SHgtEddLog$null.deviance #0.4791487
#Explained deviance of full model
1 - glm.nb.fullS$deviance / glm.nb.fullS$null.deviance #0.5029752

#Response curves here for EDD with and without logs

seq.SEdd <- seq(min(sqDataSum$Squirrel_EDD_4S), max(sqDataSum$Squirrel_EDD_4S), length.out = 100)

nd.seq.SEdd <- data.frame(Deploy.Duration = 61, Squirrel_EDD_4S = seq.SEdd, Log.in.View = "YES", Height_cm = median(sqDataSum$Height_cm))
nd.seq.SEddNo <- data.frame(Deploy.Duration = 61, Squirrel_EDD_4S = seq.SEdd, Log.in.View = "NO", Height_cm = median(sqDataSum$Height_cm))

pred.SEdd <- predict(glm.nb.SHgtEddLog, newdata = nd.seq.SEdd, se=TRUE, type = "link")
pred.SEddNo <- predict(glm.nb.SHgtEddLog, newdata = nd.seq.SEddNo, se=TRUE, type = "link")

plot(nSeqs ~ Squirrel_EDD_4S, data=sqDataSum, pch=16, col=sqDataSum$Log.in.View, las=1, xlab = "Effective Detection Distance (m)", ylab = "# Squirrel Sequences in Summer")
lines(nd.seq.SEdd$Squirrel_EDD_4S, exp(pred.SEdd$fit), col="red")
lines(nd.seq.SEdd$Squirrel_EDD_4S, exp(pred.SEdd$fit + 1.96*pred.SEdd$se.fit), lty=2, col="red")
lines(nd.seq.SEdd$Squirrel_EDD_4S, exp(pred.SEdd$fit - 1.96*pred.SEdd$se.fit), lty=2, col="red")

lines(nd.seq.SEdd$Squirrel_EDD_4S, exp(pred.SEddNo$fit), col="black")
lines(nd.seq.SEdd$Squirrel_EDD_4S, exp(pred.SEddNo$fit + 1.96*pred.SEddNo$se.fit), lty=2, col="black")
lines(nd.seq.SEdd$Squirrel_EDD_4S, exp(pred.SEddNo$fit - 1.96*pred.SEddNo$se.fit), lty=2, col="black")

#Response curves of height, with and without logs
seq.SHgt <- seq(min(sqDataSum$Height_cm), max(sqDataSum$Height_cm), length.out = 100)

nd.seq.SHgt <- data.frame(Deploy.Duration = 61, Squirrel_EDD_4S = median(sqDataSum$Squirrel_EDD_4S), Log.in.View = "YES", Height_cm = seq.SHgt)
nd.seq.SHgtNo <- data.frame(Deploy.Duration = 61, Squirrel_EDD_4S = median(sqDataSum$Squirrel_EDD_4S), Log.in.View = "NO", Height_cm = seq.SHgt)

pred.SHgtY <- predict(glm.nb.SHgtEddLog, newdata = nd.seq.SHgt, se=TRUE, type = "link")
pred.SHgtN <- predict(glm.nb.SHgtEddLog, newdata = nd.seq.SHgtNo, se=TRUE, type = "link")

plot(nSeqs ~ Height_cm, data=sqDataSum, pch=16, col=sqDataSum$Log.in.View, las=1, xlab = "Camera Height (cm)", ylab = "# Squirrel Sequences in Summer")
lines(nd.seq.SHgt$Height_cm, exp(pred.SHgtY$fit), col="red")
lines(nd.seq.SHgt$Height_cm, exp(pred.SHgtY$fit + 1.96*pred.SHgtY$se.fit), lty=2, col="red")
lines(nd.seq.SHgt$Height_cm, exp(pred.SHgtY$fit - 1.96*pred.SHgtY$se.fit), lty=2, col="red")

lines(nd.seq.SHgt$Height_cm, exp(pred.SHgtN$fit), col="black")
lines(nd.seq.SHgt$Height_cm, exp(pred.SHgtN$fit + 1.96*pred.SHgtN$se.fit), lty=2, col="black")
lines(nd.seq.SHgt$Height_cm, exp(pred.SHgtN$fit - 1.96*pred.SHgtN$se.fit), lty=2, col="black")


# Squirrel in Fall NB Full Model Selection --------------------------------

#First test to see if interaction between log and EDD improves that two variable model

glm.nb.Flog_Edd <-
  glm.nb(nSeqs ~ Log.in.View + Squirrel_EDD_4S + offset(log(Deploy.Duration)), data = sqDataFall)
AICc(glm.nb.Flog_Edd)

glm.nb.FlogXEdd <-
  glm.nb(nSeqs ~ Log.in.View * Squirrel_EDD_4S + offset(log(Deploy.Duration)), data = sqDataFall)
AICc(glm.nb.FlogXEdd) #There is support for interaction, so this has to be added to each model where these two appear together.


#Compare height alone to height quadratic
glm.nb.FHgt2 <-
  glm.nb(
    nSeqs ~ poly(Height_cm, 2, raw = T) +
      offset(log(Deploy.Duration)),
    data = sqDataFall
  )

glm.nb.Fhgt <-
  glm.nb(nSeqs ~ Height_cm + offset(log(Deploy.Duration)), data = sqDataFall)

AICc(glm.nb.FHgt2)
AICc(glm.nb.Fhgt) #no support for quadratic effect


#Full Model for Squirrels in Fall

ModListFallSq[[1]] <- glm.nb.fullF <-
  glm.nb(
    nSeqs ~ Log.in.View * Squirrel_EDD_4S + Height_cm + log10(Num_Stems) + OakDBH +
      offset(log(Deploy.Duration)),
    data = sqDataFall
  )

#Full 26 model comparison, here with log/EDD interaction

ModListFallSq <- list(NA)

ModListFallSq[[1]] <- glm.nb.FI <- glm.nb(nSeqs ~ 1 + offset(log(Deploy.Duration)), data = sqDataFall)

ModListFallSq[[2]] <- glm.nb.Flog <-
  glm.nb(nSeqs ~ Log.in.View + offset(log(Deploy.Duration)), data = sqDataFall)

ModListFallSq[[3]] <- glm.nb.Fhgt <-
  glm.nb(nSeqs ~ Height_cm + offset(log(Deploy.Duration)), data = sqDataFall)

ModListFallSq[[4]] <- glm.nb.Foak <-
  glm.nb(nSeqs ~ OakDBH + offset(log(Deploy.Duration)), data = sqDataFall)

ModListFallSq[[5]] <- glm.nb.Fstems <-
  glm.nb(nSeqs ~ log10(Num_Stems) + offset(log(Deploy.Duration)), data = sqDataFall)

ModListFallSq[[6]] <- glm.nb.Fedd <-
  glm.nb(nSeqs ~ Squirrel_EDD_4S + offset(log(Deploy.Duration)), data = sqDataFall)

ModListFallSq[[7]] <- glm.nb.Flog_stems <-
  glm.nb(nSeqs ~ Log.in.View + log10(Num_Stems) + offset(log(Deploy.Duration)), data = sqDataFall)

ModListFallSq[[8]] <- glm.nb.FlogXEdd <-
  glm.nb(nSeqs ~ Log.in.View * Squirrel_EDD_4S + offset(log(Deploy.Duration)), data = sqDataFall)

ModListFallSq[[9]] <- glm.nb.Flog_Oak <-
  glm.nb(nSeqs ~ Log.in.View + OakDBH + offset(log(Deploy.Duration)), data = sqDataFall)

ModListFallSq[[10]] <- glm.nb.Flog_Hgt <-
  glm.nb(nSeqs ~ Log.in.View + Height_cm + offset(log(Deploy.Duration)), data = sqDataFall)

ModListFallSq[[11]] <- glm.nb.FEdd_Hgt <-
  glm.nb(nSeqs ~ Squirrel_EDD_4S + Height_cm + offset(log(Deploy.Duration)), data = sqDataFall)

ModListFallSq[[12]] <-  glm.nb.FEdd_stems <-
  glm.nb(nSeqs ~ Squirrel_EDD_4S + log10(Num_Stems) + offset(log(Deploy.Duration)), data = sqDataFall)

ModListFallSq[[13]] <- glm.nb.FEdd_oak <-
  glm.nb(nSeqs ~ Squirrel_EDD_4S + OakDBH + offset(log(Deploy.Duration)), data = sqDataFall)

ModListFallSq[[14]] <- glm.nb.Fstems_hgt <-
  glm.nb(nSeqs ~ Height_cm + log10(Num_Stems) + offset(log(Deploy.Duration)), data = sqDataFall)

ModListFallSq[[15]] <- glm.nb.Fstems_oak <-
  glm.nb(nSeqs ~ OakDBH + log10(Num_Stems) + offset(log(Deploy.Duration)), data = sqDataFall)

ModListFallSq[[16]] <- glm.nb.Fhgt_oak <-
  glm.nb(nSeqs ~ Height_cm + OakDBH + offset(log(Deploy.Duration)), data = sqDataFall)

ModListFallSq[[17]] <- glm.nb.FlogStemsHgt <-
  glm.nb(
    nSeqs ~ Log.in.View + Height_cm + log10(Num_Stems) +
      offset(log(Deploy.Duration)),
    data = sqDataFall
  )

ModListFallSq[[18]] <- glm.nb.FlogStemsOak <-
  glm.nb(
    nSeqs ~ Log.in.View + OakDBH + log10(Num_Stems) +
      offset(log(Deploy.Duration)),
    data = sqDataFall
  )

ModListFallSq[[19]] <- glm.nb.FlogStemsEdd <-
  glm.nb(
    nSeqs ~ Log.in.View * Squirrel_EDD_4S + log10(Num_Stems) +
      offset(log(Deploy.Duration)),
    data = sqDataFall
  )

ModListFallSq[[20]] <- glm.nb.FHgtStemsEdd <-
  glm.nb(
    nSeqs ~ Height_cm + Squirrel_EDD_4S + log10(Num_Stems) +
      offset(log(Deploy.Duration)),
    data = sqDataFall
  )

ModListFallSq[[21]] <- glm.nb.FHgtStemsOak <-
  glm.nb(
    nSeqs ~ Height_cm + OakDBH + log10(Num_Stems) +
      offset(log(Deploy.Duration)),
    data = sqDataFall
  )

ModListFallSq[[22]] <- glm.nb.FHgtEddOak <-
  glm.nb(
    nSeqs ~ Height_cm + OakDBH + Squirrel_EDD_4S +
      offset(log(Deploy.Duration)),
    data = sqDataFall
  )
ModListFallSq[[23]] <- glm.nb.FHgtEddLog <-
  glm.nb(
    nSeqs ~ Height_cm + Log.in.View * Squirrel_EDD_4S +
      offset(log(Deploy.Duration)),
    data = sqDataFall
  )
ModListFallSq[[24]] <- glm.nb.FHgtOakLog <-
  glm.nb(
    nSeqs ~ Height_cm + Log.in.View + OakDBH +
      offset(log(Deploy.Duration)),
    data = sqDataFall
  )

ModListFallSq[[25]] <- glm.nb.FOakStemsEdd <-
  glm.nb(
    nSeqs ~ log10(Num_Stems) + Squirrel_EDD_4S + OakDBH +
      offset(log(Deploy.Duration)),
    data = sqDataFall
  )

ModListFallSq[[26]] <- glm.nb.FLogOakEdd <-
  glm.nb(
    nSeqs ~ Log.in.View * Squirrel_EDD_4S + OakDBH +
      offset(log(Deploy.Duration)),
    data = sqDataFall
  )
modnamesF <- c("Intercept","Log","Height","Oak","Stems","EDD","Log + Stems","Log * EDD","Log + Oak","Log + Height","EDD + Height","EDD + Stems","EDD + Oak","Stems + Height","Stems + Oak","Height + Oak","Log + Stems + Height","Log + Stems + Oak","Log * EDD + Stems","Height + Stems + EDD","Height + Stems + Oak","Height + EDD + Oak","Height + EDD * Log","Height + Oak + Log", "Stems + EDD + Oak", "Log * EDD + Oak")
modtabFallSq <- aictab(cand.set = ModListFallSq, modnames = modnamesF)
modtabFallSq

#Best model is clear as log * EDD
summary(glm.nb.FlogXEdd)

#Summed Model Weights Squirrels in Fall from full table. Does this still work with interaction in there? I think so, cuz it is still in the model
FallSW_edd <- sum(modtabFallSq$AICcWt[grep("EDD", modtabFallSq$Modnames)])
FallSW_OakDBH <- sum(modtabFallSq$AICcWt[grep("Oak", modtabFallSq$Modnames)])
FallSW_Stems <- sum(modtabFallSq$AICcWt[grep("Stems", modtabFallSq$Modnames)])
FallSW_Log <- sum(modtabFallSq$AICcWt[grep("Log", modtabFallSq$Modnames)])
FallSW_Height <- sum(modtabFallSq$AICcWt[grep("Height", modtabFallSq$Modnames)])


#Explained Deviance of best model
1 - glm.nb.FlogXEdd$deviance / glm.nb.FlogXEdd$null.deviance #0.49976
#Explained Deviance of best model
1 - glm.nb.fullF$deviance / glm.nb.fullF$null.deviance #0.4430968


#Response Curve for EDD log interaction----------------------------------------------------------------------- 
seq.FEdd <- seq(min(sqDataFall$Squirrel_EDD_4S), max(sqDataFall$Squirrel_EDD_4S), length.out = 100) 


nd.seq.FEdd <- data.frame(Deploy.Duration = 61, Squirrel_EDD_4S = seq.FEdd, Log.in.View = "YES")
nd.seq.FEddNo <- data.frame(Deploy.Duration = 61, Squirrel_EDD_4S = seq.FEdd, Log.in.View = "NO")

pred.FEdd <- predict(glm.nb.FlogXEdd, newdata = nd.seq.FEdd, se=TRUE, type = "link")
pred.FEddNo <- predict(glm.nb.FlogXEdd, newdata = nd.seq.FEddNo, se=TRUE, type = "link")

plot(nSeqs ~ Squirrel_EDD_4S, data=sqDataFall, pch=16, col=sqDataFall$Log.in.View, las=1, xlab = "Effective Detection Distance (m)", ylab = "# Squirrel Sequences in Fall")
lines(nd.seq.FEdd$Squirrel_EDD_4S, exp(pred.FEdd$fit), col="red")
lines(nd.seq.FEdd$Squirrel_EDD_4S, exp(pred.FEdd$fit + 1.96*pred.FEdd$se.fit), lty=2, col="red")
lines(nd.seq.FEdd$Squirrel_EDD_4S, exp(pred.FEdd$fit - 1.96*pred.FEdd$se.fit), lty=2, col="red")

lines(nd.seq.FEdd$Squirrel_EDD_4S, exp(pred.FEddNo$fit), col="black")
lines(nd.seq.FEdd$Squirrel_EDD_4S, exp(pred.FEddNo$fit + 1.96*pred.FEddNo$se.fit), lty=2, col="black")
lines(nd.seq.FEdd$Squirrel_EDD_4S, exp(pred.FEddNo$fit - 1.96*pred.FEddNo$se.fit), lty=2, col="black")

#Trying plot with ggplot instead

predVals <- rbind(as.data.frame(pred.FEdd), as.data.frame(pred.FEddNo))
allNewVals <- rbind(nd.seq.FEdd, nd.seq.FEddNo)

#I here make a dataframe combining the predictor values, with the estimated response values, then estimate lower and upper confidence limits, then convert the data back onto the real scale, from the logit scale
predPlot <- data.frame(allNewVals, predVals) %>%
  mutate(
    lcl = fit - 2*se.fit,
    ucl = fit + 2*se.fit
  ) %>%
  mutate(
    response = exp(fit),
    lcl.response = exp(lcl),
    ucl.response = exp(ucl))

#Plotting commands in ggplot. We can mess around with this further for the paper as needed. The dots represent actual data points.
SqFallEDDPlot <- ggplot(predPlot, aes(x=Squirrel_EDD_4S, y=response)) +
  # Confidence region
  geom_ribbon(aes(ymin=lcl.response, ymax=ucl.response, fill = Log.in.View), alpha=0.25, show.legend = F) +
  # Prediction Lines
  geom_line(aes(color = Log.in.View), show.legend = F) +
  scale_colour_manual("",values=c("darkolivegreen", "tomato"))+
  scale_fill_manual("",values=c("darkolivegreen", "tomato"))+
  geom_point(data = sqDataFall, aes(x=Squirrel_EDD_4S, y=nSeqs, shape=Log.in.View, color=Log.in.View), show.legend = F, size = 2.5) +
  xlab("Effective Detection Distance (m)") +   
  ylab("Count of Sequences") +
  geom_text(aes(x = 2.2, y = 220, label = "A"), size = 8)+
  theme(
    axis.line.x =  element_line(),
    axis.line.y =  element_line(),
    panel.background = element_rect(fill = "white"),
    axis.text = element_text(color = "black", size = 14),
    axis.title = element_text(size = rel(1.5))
  )

#saving this so I can use it later
myTheme <- theme(
  axis.line.x =  element_line(),
  axis.line.y =  element_line(),
  panel.background = element_rect(fill = "white"),
  axis.text = element_text(color = "black", size = 14),
  axis.title = element_text(size = rel(1.5)))
save(myTheme, file = "data/responseTheme.Rdata")

SqFallEDDPlot
#ggsave("results/SqFallEDDLog.tiff", width = 6.5, height = 4.0, units = "in" ) #saves whatever last ggplot was made

#Winter Squirrel NB Full Model Selection --------------------------------

#First test to see if interaction between log and EDD improves that two variable model

glm.nb.Wlog_Edd <-
  glm.nb(nSeqs ~ Log.in.View + Squirrel_EDD_WSp + offset(log(Deploy.Duration)), data = sqDataWin)
AICc(glm.nb.Wlog_Edd)

glm.nb.WlogXEdd <-
  glm.nb(nSeqs ~ Log.in.View * Squirrel_EDD_WSp + offset(log(Deploy.Duration)), data = sqDataWin)
AICc(glm.nb.WlogXEdd) #There is support for interaction, so this has to be added to each model where these two appear together.


#Compare height alone to height quadratic
glm.nb.WHgt2 <-
  glm.nb(
    nSeqs ~ poly(Height_cm, 2, raw = T) +
      offset(log(Deploy.Duration)),
    data = sqDataWin
  )

glm.nb.Whgt <-
  glm.nb(nSeqs ~ Height_cm + offset(log(Deploy.Duration)), data = sqDataWin)

AICc(glm.nb.WHgt2)
AICc(glm.nb.Whgt) #no support for quadratic effect

#Full Model Comparison Run

glm.nb.fullW <- glm.nb(
    nSeqs ~ Log.in.View * Squirrel_EDD_WSp + Height_cm + log10(Num_Stems) + OakDBH +
      offset(log(Deploy.Duration)),
    data = sqDataWin
  )


ModListWinSq <- list(NA)

ModListWinSq[[1]] <- glm.nb.WI <- glm.nb(nSeqs ~ 1 + offset(log(Deploy.Duration)), data = sqDataWin)

ModListWinSq[[2]] <- glm.nb.Wlog <-
  glm.nb(nSeqs ~ Log.in.View + offset(log(Deploy.Duration)), data = sqDataWin)

ModListWinSq[[3]] <- glm.nb.Whgt <-
  glm.nb(nSeqs ~ Height_cm + offset(log(Deploy.Duration)), data = sqDataWin)

ModListWinSq[[4]] <- glm.nb.Woak <-
  glm.nb(nSeqs ~ OakDBH + offset(log(Deploy.Duration)), data = sqDataWin)

ModListWinSq[[5]] <- glm.nb.Wstems <-
  glm.nb(nSeqs ~ log10(Num_Stems) + offset(log(Deploy.Duration)), data = sqDataWin)

ModListWinSq[[6]] <- glm.nb.Wedd <-
  glm.nb(nSeqs ~ Squirrel_EDD_WSp + offset(log(Deploy.Duration)), data = sqDataWin)

ModListWinSq[[7]] <- glm.nb.Wlog_stems <-
  glm.nb(nSeqs ~ Log.in.View + log10(Num_Stems) + offset(log(Deploy.Duration)), data = sqDataWin)

ModListWinSq[[8]] <-glm.nb.Wlog_Edd <-
  glm.nb(nSeqs ~ Log.in.View * Squirrel_EDD_WSp + offset(log(Deploy.Duration)), data = sqDataWin)

ModListWinSq[[9]] <- glm.nb.Wlog_Oak <-
  glm.nb(nSeqs ~ Log.in.View + OakDBH + offset(log(Deploy.Duration)), data = sqDataWin)

ModListWinSq[[10]] <- glm.nb.Wlog_Hgt <-
  glm.nb(nSeqs ~ Log.in.View + Height_cm + offset(log(Deploy.Duration)), data = sqDataWin)

ModListWinSq[[11]] <- glm.nb.WEdd_Hgt <-
  glm.nb(nSeqs ~ Squirrel_EDD_WSp + Height_cm + offset(log(Deploy.Duration)), data = sqDataWin)

ModListWinSq[[12]] <-  glm.nb.WEdd_stems <-
  glm.nb(nSeqs ~ Squirrel_EDD_WSp + log10(Num_Stems) + offset(log(Deploy.Duration)), data = sqDataWin)

ModListWinSq[[13]] <- glm.nb.WEdd_oak <-
  glm.nb(nSeqs ~ Squirrel_EDD_WSp + OakDBH + offset(log(Deploy.Duration)), data = sqDataWin)

ModListWinSq[[14]] <- glm.nb.Wstems_hgt <-
  glm.nb(nSeqs ~ Height_cm + log10(Num_Stems) + offset(log(Deploy.Duration)), data = sqDataWin)

ModListWinSq[[15]] <- glm.nb.Wstems_oak <-
  glm.nb(nSeqs ~ OakDBH + log10(Num_Stems) + offset(log(Deploy.Duration)), data = sqDataWin)

ModListWinSq[[16]] <- glm.nb.Whgt_oak <-
  glm.nb(nSeqs ~ Height_cm + OakDBH + offset(log(Deploy.Duration)), data = sqDataWin)

ModListWinSq[[17]] <- glm.nb.WlogStemsHgt <-
  glm.nb(
    nSeqs ~ Log.in.View + Height_cm + log10(Num_Stems) +
      offset(log(Deploy.Duration)),
    data = sqDataWin
  )

ModListWinSq[[18]] <- glm.nb.WlogStemsOak <-
  glm.nb(
    nSeqs ~ Log.in.View + OakDBH + log10(Num_Stems) +
      offset(log(Deploy.Duration)),
    data = sqDataWin
  )

ModListWinSq[[19]] <- glm.nb.WlogStemsEdd <-
  glm.nb(
    nSeqs ~ Log.in.View * Squirrel_EDD_WSp + log10(Num_Stems) +
      offset(log(Deploy.Duration)),
    data = sqDataWin
  )

ModListWinSq[[20]] <- glm.nb.WHgtStemsEdd <-
  glm.nb(
    nSeqs ~ Height_cm + Squirrel_EDD_WSp + log10(Num_Stems) +
      offset(log(Deploy.Duration)),
    data = sqDataWin
  )

ModListWinSq[[21]] <- glm.nb.WHgtStemsOak <-
  glm.nb(
    nSeqs ~ Height_cm + OakDBH + log10(Num_Stems) +
      offset(log(Deploy.Duration)),
    data = sqDataWin
  )

ModListWinSq[[22]] <- glm.nb.WHgtEddOak <-
  glm.nb(
    nSeqs ~ Height_cm + OakDBH + Squirrel_EDD_WSp +
      offset(log(Deploy.Duration)),
    data = sqDataWin
  )
ModListWinSq[[23]] <- glm.nb.WHgtEddLog <-
  glm.nb(
    nSeqs ~ Height_cm + Log.in.View * Squirrel_EDD_WSp +
      offset(log(Deploy.Duration)),
    data = sqDataWin
  )
ModListWinSq[[24]] <- glm.nb.WHgtOakLog <-
  glm.nb(
    nSeqs ~ Height_cm + Log.in.View + OakDBH +
      offset(log(Deploy.Duration)),
    data = sqDataWin
  )

ModListWinSq[[25]] <- glm.nb.WOakStemsEdd <-
  glm.nb(
    nSeqs ~ log10(Num_Stems) + Squirrel_EDD_WSp + OakDBH +
      offset(log(Deploy.Duration)),
    data = sqDataWin
  )

ModListWinSq[[26]] <- glm.nb.WLogOakEdd <-
  glm.nb(
    nSeqs ~ Log.in.View * Squirrel_EDD_WSp + OakDBH +
      offset(log(Deploy.Duration)),
    data = sqDataWin
  )
modnamesW <- c("Intercept","Log","Height","Oak","Stems","EDD","Log + Stems","Log * EDD","Log + Oak","Log + Height","EDD + Height","EDD + Stems","EDD + Oak","Stems + Height","Stems + Oak","Height + Oak","Log + Stems + Height","Log + Stems + Oak","Log * EDD + Stems","Height + Stems + EDD","Height + Stems + Oak","Height + EDD + Oak","Height + EDD * Log","Height + Oak + Log", "Stems + EDD + Oak", "Log * EDD + Oak")
modtabWinSq <- aictab(cand.set = ModListWinSq, modnames = modnamesW)
modtabWinSq

summary(glm.nb.Wlog) #best model only has log effects

#Summed Model Weights for Squirrels in Winter from full table. 
WinSW_edd <- sum(modtabWinSq$AICcWt[grep("EDD", modtabWinSq$Modnames)])
WinSW_OakDBH <- sum(modtabWinSq$AICcWt[grep("Oak", modtabWinSq$Modnames)])
WinSW_Stems <- sum(modtabWinSq$AICcWt[grep("Stems", modtabWinSq$Modnames)])
WinSW_Log <- sum(modtabWinSq$AICcWt[grep("Log", modtabWinSq$Modnames)])
WinSW_Height <- sum(modtabWinSq$AICcWt[grep("Height", modtabWinSq$Modnames)])



#Explained deviance of best model
1 - glm.nb.Wlog$deviance / glm.nb.Wlog$null.deviance #0.1418821
#Explained deviance of full model
1 - glm.nb.fullW$deviance / glm.nb.fullW$null.deviance #0.29


# Spring Squirrel NB Full Model Selection ---------------------------------

#First test to see if interaction between log and EDD improves that two variable model

glm.nb.SPlog_Edd <-
  glm.nb(nSeqs ~ Log.in.View + Squirrel_EDD_WSp + offset(log(Deploy.Duration)), data = sqDataSpr)
AICc(glm.nb.SPlog_Edd)

glm.nb.SPlogXEdd <-
  glm.nb(nSeqs ~ Log.in.View * Squirrel_EDD_WSp + offset(log(Deploy.Duration)), data = sqDataSpr)
AICc(glm.nb.SPlogXEdd) #There is support for interaction, so this has to be added to each model where these two appear together.


#Compare height alone to height quadratic
glm.nb.SPHgt2 <-
  glm.nb(
    nSeqs ~ poly(Height_cm, 2, raw = T) +
      offset(log(Deploy.Duration)),
    data = sqDataSpr
  )

glm.nb.SPhgt <-
  glm.nb(nSeqs ~ Height_cm + offset(log(Deploy.Duration)), data = sqDataSpr)

AICc(glm.nb.SPHgt2)
AICc(glm.nb.SPhgt) #no support for quadratic effect

#Full Spring Squirrel Model Selection

glm.nb.fullSP <-
  glm.nb(
    nSeqs ~ Log.in.View * Squirrel_EDD_WSp + Height_cm + log10(Num_Stems) + OakDBH +
      offset(log(Deploy.Duration)),
    data = sqDataSpr
  )

ModListSprSq <- list(NA)


ModListSprSq[[1]] <- glm.nb.SPI <- glm.nb(nSeqs ~ 1 + offset(log(Deploy.Duration)), data = sqDataSpr)

ModListSprSq[[2]] <- glm.nb.SPlog <-
  glm.nb(nSeqs ~ Log.in.View + offset(log(Deploy.Duration)), data = sqDataSpr)

ModListSprSq[[3]] <- glm.nb.SPhgt <-
  glm.nb(nSeqs ~ Height_cm + offset(log(Deploy.Duration)), data = sqDataSpr)

ModListSprSq[[4]] <- glm.nb.SPoak <-
  glm.nb(nSeqs ~ OakDBH + offset(log(Deploy.Duration)), data = sqDataSpr)

ModListSprSq[[5]] <- glm.nb.SPstems <-
  glm.nb(nSeqs ~ log10(Num_Stems) + offset(log(Deploy.Duration)), data = sqDataSpr)

ModListSprSq[[6]] <- glm.nb.SPedd <-
  glm.nb(nSeqs ~ Squirrel_EDD_WSp + offset(log(Deploy.Duration)), data = sqDataSpr)

ModListSprSq[[7]] <- glm.nb.SPlog_stems <-
  glm.nb(nSeqs ~ Log.in.View + log10(Num_Stems) + offset(log(Deploy.Duration)), data = sqDataSpr)

ModListSprSq[[8]] <-glm.nb.SPlog_Edd <-
  glm.nb(nSeqs ~ Log.in.View * Squirrel_EDD_WSp + offset(log(Deploy.Duration)), data = sqDataSpr)

ModListSprSq[[9]] <- glm.nb.SPlog_Oak <-
  glm.nb(nSeqs ~ Log.in.View + OakDBH + offset(log(Deploy.Duration)), data = sqDataSpr)

ModListSprSq[[10]] <- glm.nb.SPlog_Hgt <-
  glm.nb(nSeqs ~ Log.in.View + Height_cm + offset(log(Deploy.Duration)), data = sqDataSpr)

ModListSprSq[[11]] <- glm.nb.SPEdd_Hgt <-
  glm.nb(nSeqs ~ Squirrel_EDD_WSp + Height_cm + offset(log(Deploy.Duration)), data = sqDataSpr)

ModListSprSq[[12]] <-  glm.nb.SPEdd_stems <-
  glm.nb(nSeqs ~ Squirrel_EDD_WSp + log10(Num_Stems) + offset(log(Deploy.Duration)), data = sqDataSpr)

ModListSprSq[[13]] <- glm.nb.SPEdd_oak <-
  glm.nb(nSeqs ~ Squirrel_EDD_WSp + OakDBH + offset(log(Deploy.Duration)), data = sqDataSpr)

ModListSprSq[[14]] <- glm.nb.SPstems_hgt <-
  glm.nb(nSeqs ~ Height_cm + log10(Num_Stems) + offset(log(Deploy.Duration)), data = sqDataSpr)

ModListSprSq[[15]] <- glm.nb.SPstems_oak <-
  glm.nb(nSeqs ~ OakDBH + log10(Num_Stems) + offset(log(Deploy.Duration)), data = sqDataSpr)

ModListSprSq[[16]] <- glm.nb.SPhgt_oak <-
  glm.nb(nSeqs ~ Height_cm + OakDBH + offset(log(Deploy.Duration)), data = sqDataSpr)

ModListSprSq[[17]] <- glm.nb.SPlogStemsHgt <-
  glm.nb(
    nSeqs ~ Log.in.View + Height_cm + log10(Num_Stems) +
      offset(log(Deploy.Duration)),
    data = sqDataSpr
  )

ModListSprSq[[18]] <- glm.nb.SPlogStemsOak <-
  glm.nb(
    nSeqs ~ Log.in.View + OakDBH + log10(Num_Stems) +
      offset(log(Deploy.Duration)),
    data = sqDataSpr
  )

ModListSprSq[[19]] <- glm.nb.SPlogStemsEdd <-
  glm.nb(
    nSeqs ~ Log.in.View * Squirrel_EDD_WSp + log10(Num_Stems) +
      offset(log(Deploy.Duration)),
    data = sqDataSpr
  )

ModListSprSq[[20]] <- glm.nb.SPHgtStemsEdd <-
  glm.nb(
    nSeqs ~ Height_cm + Squirrel_EDD_WSp + log10(Num_Stems) +
      offset(log(Deploy.Duration)),
    data = sqDataSpr
  )

ModListSprSq[[21]] <- glm.nb.SPHgtStemsOak <-
  glm.nb(
    nSeqs ~ Height_cm + OakDBH + log10(Num_Stems) +
      offset(log(Deploy.Duration)),
    data = sqDataSpr
  )

ModListSprSq[[22]] <- glm.nb.SPHgtEddOak <-
  glm.nb(
    nSeqs ~ Height_cm + OakDBH + Squirrel_EDD_WSp +
      offset(log(Deploy.Duration)),
    data = sqDataSpr
  )
ModListSprSq[[23]] <- glm.nb.SPHgtEddLog <-
  glm.nb(
    nSeqs ~ Height_cm + Log.in.View * Squirrel_EDD_WSp +
      offset(log(Deploy.Duration)),
    data = sqDataSpr
  )
ModListSprSq[[24]] <- glm.nb.SPHgtOakLog <-
  glm.nb(
    nSeqs ~ Height_cm + Log.in.View + OakDBH +
      offset(log(Deploy.Duration)),
    data = sqDataSpr
  )
ModListSprSq[[25]] <- glm.nb.SPOakStemsEdd <-
  glm.nb(
    nSeqs ~ log10(Num_Stems) + Squirrel_EDD_WSp + OakDBH +
      offset(log(Deploy.Duration)),
    data = sqDataSpr
  )

ModListSprSq[[26]] <- glm.nb.SPLogOakEdd <-
  glm.nb(
    nSeqs ~ Log.in.View * Squirrel_EDD_WSp + OakDBH +
      offset(log(Deploy.Duration)),
    data = sqDataSpr
  )
modnamesSP <- c("Intercept","Log","Height","Oak","Stems","EDD","Log + Stems","Log * EDD","Log + Oak","Log + Height","EDD + Height","EDD + Stems","EDD + Oak","Stems + Height","Stems + Oak","Height + Oak","Log + Stems + Height","Log + Stems + Oak","Log * EDD + Stems","Height + Stems + EDD","Height + Stems + Oak","Height + EDD + Oak","Height + EDD * Log","Height + Oak + Log", "Stems + EDD + Oak", "Log * EDD + Oak")
modtabSprSq <- aictab(cand.set = ModListSprSq, modnames = modnamesSP)
modtabSprSq

summary(glm.nb.SPlog_Edd) #best model has log * EDD

#Summed Model Weights for Squirrels in Spring from full table. 
SprSW_edd <- sum(modtabSprSq$AICcWt[grep("EDD", modtabSprSq$Modnames)])
SprSW_OakDBH <- sum(modtabSprSq$AICcWt[grep("Oak", modtabSprSq$Modnames)])
SprSW_Stems <- sum(modtabSprSq$AICcWt[grep("Stems", modtabSprSq$Modnames)])
SprSW_Log <- sum(modtabSprSq$AICcWt[grep("Log", modtabSprSq$Modnames)])
SprSW_Height <- sum(modtabSprSq$AICcWt[grep("Height", modtabSprSq$Modnames)])


#Explained deviance of best model
1 - glm.nb.SPlog_Edd$deviance / glm.nb.SPlog_Edd$null.deviance #0.3287772
#Explained deviance of full model
1 - glm.nb.fullSP$deviance / glm.nb.fullSP$null.deviance #0.3807705

#Repsonse Curve of EDD with and without logs

seq.SpEdd <- seq(min(sqDataSpr$Squirrel_EDD_WSp), max(sqDataSpr$Squirrel_EDD_WSp), length.out = 100)

nd.seq.SpEdd <- data.frame(Deploy.Duration = 61, Squirrel_EDD_WSp = seq.SpEdd, Log.in.View = "YES")
nd.seq.SpEddNo <- data.frame(Deploy.Duration = 61, Squirrel_EDD_WSp = seq.SpEdd, Log.in.View = "NO")

pred.SpEdd <- predict(glm.nb.SPlog_Edd, newdata = nd.seq.SpEdd, se=TRUE, type = "link")
pred.SpEddNo <- predict(glm.nb.SPlog_Edd, newdata = nd.seq.SpEddNo, se=TRUE, type = "link")

predValsSp <- rbind(as.data.frame(pred.SpEdd), as.data.frame(pred.SpEddNo))
allNewValsSp <- rbind(nd.seq.SpEdd, nd.seq.SpEddNo)

#I here make a dataframe combining the predictor values, with the estimated response values, then estimate lower and upper confidence limits, then convert the data back onto the real scale, from the log scale
predPlotSp <- data.frame(allNewValsSp, predValsSp) %>%
  mutate(
    lcl = fit - 2*se.fit,
    ucl = fit + 2*se.fit
  ) %>%
  mutate(
    response = exp(fit),
    lcl.response = exp(lcl),
    ucl.response = exp(ucl))

#Plotting commands in ggplot. We can mess around with this further for the paper as needed. The dots represent actual data points.

SqSpEDDPlot <- ggplot(predPlotSp, aes(x=Squirrel_EDD_WSp, y=response)) +
  # Confidence region
  geom_ribbon(aes(ymin=lcl.response, ymax=ucl.response, fill = Log.in.View), alpha=0.25, show.legend = F) +
  # Prediction Lines
  geom_line(aes(color = Log.in.View), show.legend = F) +
  scale_colour_manual("",values=c("darkolivegreen", "tomato"))+
  scale_fill_manual("",values=c("darkolivegreen", "tomato"))+
  geom_point(data = sqDataSpr, aes(x=Squirrel_EDD_WSp, y=nSeqs, shape=Log.in.View, color=Log.in.View), show.legend = F, size = 2.5) +
  xlab("Effective Detection Distance (m)") +   
  ylab("Count of Sequences") +
  geom_text(aes(x = 2.3, y = 220, label = "B"), size = 8) +
  myTheme

SqSpEDDPlot
ggsave("results/SqSpringEDDLog.tiff", width = 6.5, height = 4.0, units = "in" ) #saves whatever last ggplot was made