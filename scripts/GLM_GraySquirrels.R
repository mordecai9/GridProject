###Regression Analysis of SCBI Camera Grid Data

rm(list = ls())

require(AICcmodavg)
require(MASS) #for glm.nb function
source("scripts/pairsPannelFunctions.r")


# Data Import and Exploration - Gray Squirrels ---------------------------------------

load("data/grsqDataSum.RData")
load("data/grsqDataFall.RData")
load("data/grsqDataWin.RData")
load("data/grsqDataSpr.RData")

pairs(grsqDataSum[,c(4,7,13,14,15, 17)], diag.panel = panel.hist, lower.panel = panel.smooth, upper.panel = panel.cor)
str(grsqDataSum)

hist(grsqDataSum$Num_Stems)
hist(log10(grsqDataSum$Num_Stems)) #because one site has so many stems, we need to log10 transform this parameter for all models
hist(grsqDataSum$Squirrel_EDD_4S)

#Univariate exploratory plots for each of our 5 covariates
par.default <- par(no.readonly = T)
par(mfrow = c(2,2))
plot(nSeqs ~ Squirrel_EDD_4S, data = grsqDataSum, main = "Summer")
plot(nSeqs ~ Squirrel_EDD_4S, data = grsqDataFall, main = "Fall")
plot(nSeqs ~ Squirrel_EDD_WSp, data = grsqDataWin, main = "Winter")
plot(nSeqs ~ Squirrel_EDD_WSp, data = grsqDataSpr, main = "Spring")

plot(nSeqs ~ OakDBH, data = grsqDataSum, main = "Summer")
plot(nSeqs ~ OakDBH, data = grsqDataFall, main = "Fall")
plot(nSeqs ~ OakDBH, data = grsqDataWin, main = "Winter")
plot(nSeqs ~ OakDBH, data = grsqDataSpr, main = "Spring")

plot(nSeqs ~ log10(Num_Stems), data = grsqDataSum, main = "Summer")
plot(nSeqs ~ log10(Num_Stems), data = grsqDataFall, main = "Fall")
plot(nSeqs ~ log10(Num_Stems), data = grsqDataWin, main = "Winter")
plot(nSeqs ~ log10(Num_Stems), data = grsqDataSpr, main = "Spring")

#These seems to show a quadratic relationship actually, which might make sense
plot(nSeqs ~ Height_cm, data = grsqDataSum, main = "Summer")
plot(nSeqs ~ Height_cm, data = grsqDataFall, main = "Fall")
plot(nSeqs ~ Height_cm, data = grsqDataWin, main = "Winter")
plot(nSeqs ~ Height_cm, data = grsqDataSpr, main = "Spring")

plot(nSeqs ~ Log.in.View, data = grsqDataSum, main = "Summer")
plot(nSeqs ~ Log.in.View, data = grsqDataFall, main = "Fall")
plot(nSeqs ~ Log.in.View, data = grsqDataWin, main = "Winter")
plot(nSeqs ~ Log.in.View, data = grsqDataSpr, main = "Spring")
par(par.default)

hist(grsqDataSum$Deploy.Duration, breaks = 5)

#Looking at how many cameras picked up bear in each season
barplot(grsqDataSum$nSeqs, xlab = "Summer", ylab = "# of Seqs per camera")
barplot(grsqDataFall$nSeqs, xlab = "Fall", ylab = "# of Seqs per camera")
barplot(grsqDataWin$nSeqs, xlab = "Winter", ylab = "# of Seqs per camera")
barplot(grsqDataSpr$nSeqs, xlab = "Spring", ylab = "# of Seqs per camera")


# Full Poisson GLM and Overdispersion -------------------------------------

#Exploratory Poisson model for Summer data to see if a negative binomial model is needed due to potential overdispersion
glm.po.fullS <- glm(nSeqs ~ Log.in.View + Squirrel_EDD_4S + Height_cm + log10(Num_Stems) + OakDBH +offset(log(Deploy.Duration)), data = grsqDataSum, family = poisson)
summary(glm.po.fullS)

#test for overdispersion indicates we should be using negative binomial (2.31)
sum.po.fullS <- summary(glm.po.fullS)
phi.po.fullS <- sum.po.fullS$deviance/sum.po.fullS$df.residual 

#Exploratory Poisson model for Fall data
glm.po.fullF <- glm(nSeqs ~ Log.in.View + Squirrel_EDD_4S + Height_cm + log10(Num_Stems) + OakDBH +offset(log(Deploy.Duration)), data = grsqDataFall, family = poisson)
summary(glm.po.fullF)

#test for overdispersion indicates we should be using negative binomial (5.385) 
sum.po.fullF <- summary(glm.po.fullF)
phi.po.fullF <- sum.po.fullF$deviance/sum.po.fullF$df.residual 

#Exploratory Poisson model for Winter data
glm.po.fullW <- glm(nSeqs ~ Log.in.View + Squirrel_EDD_WSp + Height_cm + log10(Num_Stems) + OakDBH +offset(log(Deploy.Duration)), data = grsqDataWin, family = poisson)
summary(glm.po.fullW)

#test for overdispersion indicates we should be using negative binomial  (11.61)
sum.po.fullW <- summary(glm.po.fullW)
phi.po.fullW <- sum.po.fullW$deviance/sum.po.fullW$df.residual 

#Exploratory Poisson model for Spring data
glm.po.fullSp <- glm(nSeqs ~ Log.in.View + Squirrel_EDD_WSp + Height_cm + log10(Num_Stems) + OakDBH +offset(log(Deploy.Duration)), data = grsqDataSpr, family = poisson)
summary(glm.po.fullSp)

#test for overdispersion indicates we should be using negative binomial (8.43)
sum.po.fullSp <- summary(glm.po.fullSp)
phi.po.fullSp <- sum.po.fullSp$deviance/sum.po.fullSp$df.residual 


# Negative Binomial Regression - Gray Squirrels in Summer --------------------------------
#First test to see if interaction between log and EDD improves that two variable model

glm.nb.Slog_Edd <-
  glm.nb(nSeqs ~ Log.in.View + Squirrel_EDD_4S + offset(log(Deploy.Duration)), data = grsqDataSum)
AICc(glm.nb.Slog_Edd)

glm.nb.SlogXEdd <-
  glm.nb(nSeqs ~ Log.in.View * Squirrel_EDD_4S + offset(log(Deploy.Duration)), data = grsqDataSum)
AICc(glm.nb.SlogXEdd) #No support for interaction 


#Compare height alone to height quadratic
glm.nb.SHgt2 <-
  glm.nb(
    nSeqs ~ poly(Height_cm, 2, raw = T) +
      offset(log(Deploy.Duration)),
    data = grsqDataSum
  )

glm.nb.Shgt <-
  glm.nb(
    nSeqs ~ Height_cm +
      offset(log(Deploy.Duration)),
    data = grsqDataSum
  )

AICc(glm.nb.SHgt2)
AICc(glm.nb.Shgt) #no support for quadratic effect

#Full Model Selection for Gray Squirrels in Summer

glm.nb.fullS <-
  glm.nb(
    nSeqs ~ Log.in.View + Squirrel_EDD_4S + Height_cm + log10(Num_Stems) + OakDBH +
      offset(log(Deploy.Duration)),
    data = grsqDataSum
  )


ModListSumGSq <- list(NA)


ModListSumGSq[[1]] <- glm.nb.SI <- glm.nb(nSeqs ~ 1 + offset(log(Deploy.Duration)), data = grsqDataSum)

ModListSumGSq[[2]] <- glm.nb.Slog <-
  glm.nb(nSeqs ~ Log.in.View + offset(log(Deploy.Duration)), data = grsqDataSum)

ModListSumGSq[[3]] <- glm.nb.Shgt <-
  glm.nb(nSeqs ~ Height_cm + offset(log(Deploy.Duration)), data = grsqDataSum)

ModListSumGSq[[4]] <- glm.nb.Soak <-
  glm.nb(nSeqs ~ OakDBH + offset(log(Deploy.Duration)), data = grsqDataSum)

ModListSumGSq[[5]] <- glm.nb.Sstems <-
  glm.nb(nSeqs ~ log10(Num_Stems) + offset(log(Deploy.Duration)), data = grsqDataSum)

ModListSumGSq[[6]] <- glm.nb.Sedd <-
  glm.nb(nSeqs ~ Squirrel_EDD_4S + offset(log(Deploy.Duration)), data = grsqDataSum)

ModListSumGSq[[7]] <- glm.nb.Slog_stems <-
  glm.nb(nSeqs ~ Log.in.View + log10(Num_Stems) + offset(log(Deploy.Duration)), data = grsqDataSum)

ModListSumGSq[[8]] <-glm.nb.Slog_Edd <-
  glm.nb(nSeqs ~ Log.in.View + Squirrel_EDD_4S + offset(log(Deploy.Duration)), data = grsqDataSum)

ModListSumGSq[[9]] <- glm.nb.Slog_Oak <-
  glm.nb(nSeqs ~ Log.in.View + OakDBH + offset(log(Deploy.Duration)), data = grsqDataSum)

ModListSumGSq[[10]] <- glm.nb.Slog_Hgt <-
  glm.nb(nSeqs ~ Log.in.View + Height_cm + offset(log(Deploy.Duration)), data = grsqDataSum)

ModListSumGSq[[11]] <- glm.nb.SEdd_Hgt <-
  glm.nb(nSeqs ~ Squirrel_EDD_4S + Height_cm + offset(log(Deploy.Duration)), data = grsqDataSum)

ModListSumGSq[[12]] <-  glm.nb.SEdd_stems <-
  glm.nb(nSeqs ~ Squirrel_EDD_4S + log10(Num_Stems) + offset(log(Deploy.Duration)), data = grsqDataSum)

ModListSumGSq[[13]] <- glm.nb.SEdd_oak <-
  glm.nb(nSeqs ~ Squirrel_EDD_4S + OakDBH + offset(log(Deploy.Duration)), data = grsqDataSum)

ModListSumGSq[[14]] <- glm.nb.Sstems_hgt <-
  glm.nb(nSeqs ~ Height_cm + log10(Num_Stems) + offset(log(Deploy.Duration)), data = grsqDataSum)

ModListSumGSq[[15]] <- glm.nb.Sstems_oak <-
  glm.nb(nSeqs ~ OakDBH + log10(Num_Stems) + offset(log(Deploy.Duration)), data = grsqDataSum)

ModListSumGSq[[16]] <- glm.nb.Shgt_oak <-
  glm.nb(nSeqs ~ Height_cm + OakDBH + offset(log(Deploy.Duration)), data = grsqDataSum)

ModListSumGSq[[17]] <- glm.nb.SlogStemsHgt <-
  glm.nb(
    nSeqs ~ Log.in.View + Height_cm + log10(Num_Stems) +
      offset(log(Deploy.Duration)),
    data = grsqDataSum
  )

ModListSumGSq[[18]] <- glm.nb.SlogStemsOak <-
  glm.nb(
    nSeqs ~ Log.in.View + OakDBH + log10(Num_Stems) +
      offset(log(Deploy.Duration)),
    data = grsqDataSum
  )

ModListSumGSq[[19]] <- glm.nb.SlogStemsEdd <-
  glm.nb(
    nSeqs ~ Log.in.View + Squirrel_EDD_4S + log10(Num_Stems) +
      offset(log(Deploy.Duration)),
    data = grsqDataSum
  )

ModListSumGSq[[20]] <- glm.nb.SHgtStemsEdd <-
  glm.nb(
    nSeqs ~ Height_cm + Squirrel_EDD_4S + log10(Num_Stems) +
      offset(log(Deploy.Duration)),
    data = grsqDataSum
  )

ModListSumGSq[[21]] <- glm.nb.SHgtStemsOak <-
  glm.nb(
    nSeqs ~ Height_cm + OakDBH + log10(Num_Stems) +
      offset(log(Deploy.Duration)),
    data = grsqDataSum
  )

ModListSumGSq[[22]] <- glm.nb.SHgtEddOak <-
  glm.nb(
    nSeqs ~ Height_cm + OakDBH + Squirrel_EDD_4S +
      offset(log(Deploy.Duration)),
    data = grsqDataSum
  )
ModListSumGSq[[23]] <- glm.nb.SHgtEddLog <-
  glm.nb(
    nSeqs ~ Height_cm + Log.in.View + Squirrel_EDD_4S +
      offset(log(Deploy.Duration)),
    data = grsqDataSum
  )
ModListSumGSq[[24]] <- glm.nb.SHgtOakLog <-
  glm.nb(
    nSeqs ~ Height_cm + Log.in.View + OakDBH +
      offset(log(Deploy.Duration)),
    data = grsqDataSum
  )

ModListSumGSq[[25]] <- glm.nb.SOakStemsEdd <-
  glm.nb(
    nSeqs ~ log10(Num_Stems) + Squirrel_EDD_4S + OakDBH +
      offset(log(Deploy.Duration)),
    data = grsqDataSum
  )

ModListSumGSq[[26]] <- glm.nb.SLogOakEdd <-
  glm.nb(
    nSeqs ~ Log.in.View + Squirrel_EDD_4S + OakDBH +
      offset(log(Deploy.Duration)),
    data = grsqDataSum
  )
modnamesS <- c("Intercept","Log","Height","Oak","Stems","EDD","Log + Stems","Log + EDD","Log + Oak","Log + Height","EDD + Height","EDD + Stems","EDD + Oak","Stems + Height","Stems + Oak","Height + Oak","Log + Stems + Height","Log + Stems + Oak","Log + Stems + EDD","Height + Stems + EDD","Height + Stems + Oak","Height + EDD + Oak","Height + EDD + Log","Height + Oak + Log", "Stems + EDD + Oak", "Log + EDD + Oak")
modtabSumGSq <- aictab(cand.set = ModListSumGSq, modnames = modnamesS)
modtabSumGSq

#Best model has EDD and log
summary(glm.nb.Slog_Edd)
AICc(glm.nb.Slog_Edd)

#Summed Model Weights for Gray Squirrels in Summer from full table. 
SumSW_edd <- sum(modtabSumGSq$AICcWt[grep("EDD", modtabSumGSq$Modnames)])
SumSW_OakDBH <- sum(modtabSumGSq$AICcWt[grep("Oak", modtabSumGSq$Modnames)])
SumSW_Stems <- sum(modtabSumGSq$AICcWt[grep("Stems", modtabSumGSq$Modnames)])
SumSW_Log <- sum(modtabSumGSq$AICcWt[grep("Log", modtabSumGSq$Modnames)])
SumSW_Height <- sum(modtabSumGSq$AICcWt[grep("Height", modtabSumGSq$Modnames)])


#See if this model is better with EDD and log interaction out of curiosity
glm.nb.SEddXLog <-
  glm.nb(
    nSeqs ~ Log.in.View * Squirrel_EDD_4S +
      offset(log(Deploy.Duration)),
    data = grsqDataSum
  )
AICc(glm.nb.SEddXLog) #not worth having the interaction here
summary(glm.nb.SEddXLog)


#Explained deviance of best model
1 - glm.nb.Slog_Edd$deviance / glm.nb.Slog_Edd$null.deviance #0.3027614
#Explained deviance of full model
1 - glm.nb.fullS$deviance / glm.nb.fullS$null.deviance #0.4038415

#Response curves here for EDD with and without logs

seq.SEdd <- seq(min(grsqDataSum$Squirrel_EDD_4S), max(grsqDataSum$Squirrel_EDD_4S), length.out = 100)

nd.seq.SEdd <- data.frame(Deploy.Duration = 61, Squirrel_EDD_4S = seq.SEdd, Log.in.View = "YES", Height_cm = median(grsqDataSum$Height_cm))
nd.seq.SEddNo <- data.frame(Deploy.Duration = 61, Squirrel_EDD_4S = seq.SEdd, Log.in.View = "NO", Height_cm = median(grsqDataSum$Height_cm))

pred.SEdd <- predict(glm.nb.SHgtEddLog, newdata = nd.seq.SEdd, se=TRUE, type = "link")
pred.SEddNo <- predict(glm.nb.SHgtEddLog, newdata = nd.seq.SEddNo, se=TRUE, type = "link")

plot(nSeqs ~ Squirrel_EDD_4S, data=grsqDataSum, pch=16, col=grsqDataSum$Log.in.View, las=1, xlab = "Effective Detection Distance (m)", ylab = "# gr Squirrel Sequences in Summer")
lines(nd.seq.SEdd$Squirrel_EDD_4S, exp(pred.SEdd$fit), col="red")
lines(nd.seq.SEdd$Squirrel_EDD_4S, exp(pred.SEdd$fit + 1.96*pred.SEdd$se.fit), lty=2, col="red")
lines(nd.seq.SEdd$Squirrel_EDD_4S, exp(pred.SEdd$fit - 1.96*pred.SEdd$se.fit), lty=2, col="red")

lines(nd.seq.SEdd$Squirrel_EDD_4S, exp(pred.SEddNo$fit), col="black")
lines(nd.seq.SEdd$Squirrel_EDD_4S, exp(pred.SEddNo$fit + 1.96*pred.SEddNo$se.fit), lty=2, col="black")
lines(nd.seq.SEdd$Squirrel_EDD_4S, exp(pred.SEddNo$fit - 1.96*pred.SEddNo$se.fit), lty=2, col="black")



# Gray Squirrel in Fall NB Full Model Selection --------------------------------

#First test to see if interaction between log and EDD improves that two variable model

glm.nb.Flog_Edd <-
  glm.nb(nSeqs ~ Log.in.View + Squirrel_EDD_4S + offset(log(Deploy.Duration)), data = grsqDataFall)
AICc(glm.nb.Flog_Edd)

glm.nb.FlogXEdd <-
  glm.nb(nSeqs ~ Log.in.View * Squirrel_EDD_4S + offset(log(Deploy.Duration)), data = grsqDataFall)
AICc(glm.nb.FlogXEdd) #There is support for interaction, so this needs to be added into all models below.


#Compare height alone to height quadratic
glm.nb.FHgt2 <-
  glm.nb(
    nSeqs ~ poly(Height_cm, 2, raw = T) +
      offset(log(Deploy.Duration)),
    data = grsqDataFall
  )

glm.nb.Fhgt <-
  glm.nb(nSeqs ~ Height_cm + offset(log(Deploy.Duration)), data = grsqDataFall)

AICc(glm.nb.FHgt2)
AICc(glm.nb.Fhgt) #no support for quadratic effect

#Full Model Selection for Fall Gray Squirrels


glm.nb.fullF <-
  glm.nb(
    nSeqs ~ Log.in.View * Squirrel_EDD_4S + Height_cm + log10(Num_Stems) + OakDBH +
      offset(log(Deploy.Duration)),
    data = grsqDataFall
  )

ModListFallGSq <- list(NA)


ModListFallGSq[[1]] <- glm.nb.FI <- glm.nb(nSeqs ~ 1 + offset(log(Deploy.Duration)), data = grsqDataFall)

ModListFallGSq[[2]] <- glm.nb.Flog <-
  glm.nb(nSeqs ~ Log.in.View + offset(log(Deploy.Duration)), data = grsqDataFall)

ModListFallGSq[[3]] <- glm.nb.Fhgt <-
  glm.nb(nSeqs ~ Height_cm + offset(log(Deploy.Duration)), data = grsqDataFall)

ModListFallGSq[[4]] <- glm.nb.Foak <-
  glm.nb(nSeqs ~ OakDBH + offset(log(Deploy.Duration)), data = grsqDataFall)

ModListFallGSq[[5]] <- glm.nb.Fstems <-
  glm.nb(nSeqs ~ log10(Num_Stems) + offset(log(Deploy.Duration)), data = grsqDataFall)

ModListFallGSq[[6]] <- glm.nb.Fedd <-
  glm.nb(nSeqs ~ Squirrel_EDD_4S + offset(log(Deploy.Duration)), data = grsqDataFall)

ModListFallGSq[[7]] <- glm.nb.Flog_stems <-
  glm.nb(nSeqs ~ Log.in.View + log10(Num_Stems) + offset(log(Deploy.Duration)), data = grsqDataFall)

ModListFallGSq[[8]] <- glm.nb.Flog_Edd <-
  glm.nb(nSeqs ~ Log.in.View * Squirrel_EDD_4S + offset(log(Deploy.Duration)), data = grsqDataFall)

ModListFallGSq[[9]] <- glm.nb.Flog_Oak <-
  glm.nb(nSeqs ~ Log.in.View + OakDBH + offset(log(Deploy.Duration)), data = grsqDataFall)

ModListFallGSq[[10]] <- glm.nb.Flog_Hgt <-
  glm.nb(nSeqs ~ Log.in.View + Height_cm + offset(log(Deploy.Duration)), data = grsqDataFall)

ModListFallGSq[[11]] <- glm.nb.FEdd_Hgt <-
  glm.nb(nSeqs ~ Squirrel_EDD_4S + Height_cm + offset(log(Deploy.Duration)), data = grsqDataFall)

ModListFallGSq[[12]] <-  glm.nb.FEdd_stems <-
  glm.nb(nSeqs ~ Squirrel_EDD_4S + log10(Num_Stems) + offset(log(Deploy.Duration)), data = grsqDataFall)

ModListFallGSq[[13]] <- glm.nb.FEdd_oak <-
  glm.nb(nSeqs ~ Squirrel_EDD_4S + OakDBH + offset(log(Deploy.Duration)), data = grsqDataFall)

ModListFallGSq[[14]] <- glm.nb.Fstems_hgt <-
  glm.nb(nSeqs ~ Height_cm + log10(Num_Stems) + offset(log(Deploy.Duration)), data = grsqDataFall)

ModListFallGSq[[15]] <- glm.nb.Fstems_oak <-
  glm.nb(nSeqs ~ OakDBH + log10(Num_Stems) + offset(log(Deploy.Duration)), data = grsqDataFall)

ModListFallGSq[[16]] <- glm.nb.Fhgt_oak <-
  glm.nb(nSeqs ~ Height_cm + OakDBH + offset(log(Deploy.Duration)), data = grsqDataFall)

ModListFallGSq[[17]] <- glm.nb.FlogStemsHgt <-
  glm.nb(
    nSeqs ~ Log.in.View + Height_cm + log10(Num_Stems) +
      offset(log(Deploy.Duration)),
    data = grsqDataFall
  )

ModListFallGSq[[18]] <- glm.nb.FlogStemsOak <-
  glm.nb(
    nSeqs ~ Log.in.View + OakDBH + log10(Num_Stems) +
      offset(log(Deploy.Duration)),
    data = grsqDataFall
  )

ModListFallGSq[[19]] <- glm.nb.FlogStemsEdd <-
  glm.nb(
    nSeqs ~ Log.in.View * Squirrel_EDD_4S + log10(Num_Stems) +
      offset(log(Deploy.Duration)),
    data = grsqDataFall
  )

ModListFallGSq[[20]] <- glm.nb.FHgtStemsEdd <-
  glm.nb(
    nSeqs ~ Height_cm + Squirrel_EDD_4S + log10(Num_Stems) +
      offset(log(Deploy.Duration)),
    data = grsqDataFall
  )

ModListFallGSq[[21]] <- glm.nb.FHgtStemsOak <-
  glm.nb(
    nSeqs ~ Height_cm + OakDBH + log10(Num_Stems) +
      offset(log(Deploy.Duration)),
    data = grsqDataFall
  )

ModListFallGSq[[22]] <- glm.nb.FHgtEddOak <-
  glm.nb(
    nSeqs ~ Height_cm + OakDBH + Squirrel_EDD_4S +
      offset(log(Deploy.Duration)),
    data = grsqDataFall
  )
ModListFallGSq[[23]] <- glm.nb.FHgtEddLog <-
  glm.nb(
    nSeqs ~ Height_cm + Log.in.View * Squirrel_EDD_4S +
      offset(log(Deploy.Duration)),
    data = grsqDataFall
  )
ModListFallGSq[[24]] <- glm.nb.FHgtOakLog <-
  glm.nb(
    nSeqs ~ Height_cm + Log.in.View + OakDBH +
      offset(log(Deploy.Duration)),
    data = grsqDataFall
  )
ModListFallGSq[[25]] <- glm.nb.FOakStemsEdd <-
  glm.nb(
    nSeqs ~ log10(Num_Stems) + Squirrel_EDD_4S + OakDBH +
      offset(log(Deploy.Duration)),
    data = grsqDataFall
  )

ModListFallGSq[[26]] <- glm.nb.FLogOakEdd <-
  glm.nb(
    nSeqs ~ Log.in.View * Squirrel_EDD_4S + OakDBH +
      offset(log(Deploy.Duration)),
    data = grsqDataFall
  )
modnamesF <- c("Intercept","Log","Height","Oak","Stems","EDD","Log + Stems","Log * EDD","Log + Oak","Log + Height","EDD + Height","EDD + Stems","EDD + Oak","Stems + Height","Stems + Oak","Height + Oak","Log + Stems + Height","Log + Stems + Oak","Log * EDD + Stems","Height + Stems + EDD","Height + Stems + Oak","Height + EDD + Oak","Height + EDD * Log","Height + Oak + Log", "Stems + EDD + Oak", "Log * EDD + Oak")
modtabFallGSq <- aictab(cand.set = ModListFallGSq, modnames = modnamesF)
modtabFallGSq

#Best model is Log * EDD 
summary(glm.nb.FlogXEdd)

#Summed Model Weights Squirrels in Fall from full table. Note interaction not included at the moment
FallSW_edd <- sum(modtabFallGSq$AICcWt[grep("EDD", modtabFallGSq$Modnames)])
FallSW_OakDBH <- sum(modtabFallGSq$AICcWt[grep("Oak", modtabFallGSq$Modnames)])
FallSW_Stems <- sum(modtabFallGSq$AICcWt[grep("Stems", modtabFallGSq$Modnames)])
FallSW_Log <- sum(modtabFallGSq$AICcWt[grep("Log", modtabFallGSq$Modnames)])
FallSW_Height <- sum(modtabFallGSq$AICcWt[grep("Height", modtabFallGSq$Modnames)])

#Explained Deviance of best model
1 - glm.nb.FlogXEdd$deviance / glm.nb.FlogXEdd$null.deviance #0.4755511
#Explained Deviance of best model
1 - glm.nb.fullF$deviance / glm.nb.fullF$null.deviance #0.2699834

#Response Curve for EDD impact, with and without log present. Possibly worth including this in the paper.

seq.FEdd <- seq(min(grsqDataFall$Squirrel_EDD_4S), max(grsqDataFall$Squirrel_EDD_4S), length.out = 100)

nd.seq.FEdd <- data.frame(Deploy.Duration = 61, Squirrel_EDD_4S = seq.FEdd, Log.in.View = "YES")
nd.seq.FEddNo <- data.frame(Deploy.Duration = 61, Squirrel_EDD_4S = seq.FEdd, Log.in.View = "NO")

pred.FEdd <- predict(glm.nb.FlogXEdd, newdata = nd.seq.FEdd, se=TRUE, type = "link")
pred.FEddNo <- predict(glm.nb.FlogXEdd, newdata = nd.seq.FEddNo, se=TRUE, type = "link")

plot(nSeqs ~ Squirrel_EDD_4S, data=grsqDataFall, pch=16, col=grsqDataFall$Log.in.View, las=1, xlab = "Effective Detection Distance (m)", ylab = "# Squirrel Sequences in Fall")
lines(nd.seq.FEdd$Squirrel_EDD_4S, exp(pred.FEdd$fit), col="red")
lines(nd.seq.FEdd$Squirrel_EDD_4S, exp(pred.FEdd$fit + 1.96*pred.FEdd$se.fit), lty=2, col="red")
lines(nd.seq.FEdd$Squirrel_EDD_4S, exp(pred.FEdd$fit - 1.96*pred.FEdd$se.fit), lty=2, col="red")

lines(nd.seq.FEdd$Squirrel_EDD_4S, exp(pred.FEddNo$fit), col="black")
lines(nd.seq.FEdd$Squirrel_EDD_4S, exp(pred.FEddNo$fit + 1.96*pred.FEddNo$se.fit), lty=2, col="black")
lines(nd.seq.FEdd$Squirrel_EDD_4S, exp(pred.FEddNo$fit - 1.96*pred.FEddNo$se.fit), lty=2, col="black")

#Winter Gray Squirrel NB Full Model Selection --------------------------------

#First test to see if interaction between log and EDD improves that two variable model

glm.nb.Wlog_Edd <-
  glm.nb(nSeqs ~ Log.in.View + Squirrel_EDD_WSp + offset(log(Deploy.Duration)), data = grsqDataWin)
AICc(glm.nb.Wlog_Edd)

glm.nb.WlogXEdd <-
  glm.nb(nSeqs ~ Log.in.View * Squirrel_EDD_WSp + offset(log(Deploy.Duration)), data = grsqDataWin)
AICc(glm.nb.WlogXEdd) #There is support for interaction, so it must go into all models below.


#Compare height alone to height quadratic
glm.nb.WHgt2 <-
  glm.nb(
    nSeqs ~ poly(Height_cm, 2, raw = T) +
      offset(log(Deploy.Duration)),
    data = grsqDataWin
  )

glm.nb.Whgt <-
  glm.nb(nSeqs ~ Height_cm + offset(log(Deploy.Duration)), data = grsqDataWin)

AICc(glm.nb.WHgt2)
AICc(glm.nb.Whgt) #there is support for quadratic effect on height

#Full Model Comparison Run Winter Gray Squirrels


ModListWinGSq[[1]] <- glm.nb.fullW <-
  glm.nb(
    nSeqs ~ Log.in.View * Squirrel_EDD_WSp + poly(Height_cm, 2, raw = T) + log10(Num_Stems) + OakDBH +
      offset(log(Deploy.Duration)),
    data = grsqDataWin
  )

ModListWinGSq <- list(NA)

ModListWinGSq[[1]] <- glm.nb.WI <- glm.nb(nSeqs ~ 1 + offset(log(Deploy.Duration)), data = grsqDataWin)

ModListWinGSq[[2]] <- glm.nb.Wlog <-
  glm.nb(nSeqs ~ Log.in.View + offset(log(Deploy.Duration)), data = grsqDataWin)

ModListWinGSq[[3]] <- glm.nb.Whgt <-
  glm.nb(nSeqs ~ poly(Height_cm, 2, raw = T) + offset(log(Deploy.Duration)), data = grsqDataWin)

ModListWinGSq[[4]] <- glm.nb.Woak <-
  glm.nb(nSeqs ~ OakDBH + offset(log(Deploy.Duration)), data = grsqDataWin)

ModListWinGSq[[5]] <- glm.nb.Wstems <-
  glm.nb(nSeqs ~ log10(Num_Stems) + offset(log(Deploy.Duration)), data = grsqDataWin)

ModListWinGSq[[6]] <- glm.nb.Wedd <-
  glm.nb(nSeqs ~ Squirrel_EDD_WSp + offset(log(Deploy.Duration)), data = grsqDataWin)

ModListWinGSq[[7]] <- glm.nb.Wlog_stems <-
  glm.nb(nSeqs ~ Log.in.View + log10(Num_Stems) + offset(log(Deploy.Duration)), data = grsqDataWin)

ModListWinGSq[[8]] <-glm.nb.Wlog_Edd <-
  glm.nb(nSeqs ~ Log.in.View * Squirrel_EDD_WSp + offset(log(Deploy.Duration)), data = grsqDataWin)

ModListWinGSq[[9]] <- glm.nb.Wlog_Oak <-
  glm.nb(nSeqs ~ Log.in.View + OakDBH + offset(log(Deploy.Duration)), data = grsqDataWin)

ModListWinGSq[[10]] <- glm.nb.Wlog_Hgt <-
  glm.nb(nSeqs ~ Log.in.View + poly(Height_cm, 2, raw = T) + offset(log(Deploy.Duration)), data = grsqDataWin)

ModListWinGSq[[11]] <- glm.nb.WEdd_Hgt <-
  glm.nb(nSeqs ~ Squirrel_EDD_WSp + poly(Height_cm, 2, raw = T) + offset(log(Deploy.Duration)), data = grsqDataWin)

ModListWinGSq[[12]] <-  glm.nb.WEdd_stems <-
  glm.nb(nSeqs ~ Squirrel_EDD_WSp + log10(Num_Stems) + offset(log(Deploy.Duration)), data = grsqDataWin)

ModListWinGSq[[13]] <- glm.nb.WEdd_oak <-
  glm.nb(nSeqs ~ Squirrel_EDD_WSp + OakDBH + offset(log(Deploy.Duration)), data = grsqDataWin)

ModListWinGSq[[14]] <- glm.nb.Wstems_hgt <-
  glm.nb(nSeqs ~ poly(Height_cm, 2, raw = T) + log10(Num_Stems) + offset(log(Deploy.Duration)), data = grsqDataWin)

ModListWinGSq[[15]] <- glm.nb.Wstems_oak <-
  glm.nb(nSeqs ~ OakDBH + log10(Num_Stems) + offset(log(Deploy.Duration)), data = grsqDataWin)

ModListWinGSq[[16]] <- glm.nb.Whgt_oak <-
  glm.nb(nSeqs ~ poly(Height_cm, 2, raw = T) + OakDBH + offset(log(Deploy.Duration)), data = grsqDataWin)

ModListWinGSq[[17]] <- glm.nb.WlogStemsHgt <-
  glm.nb(
    nSeqs ~ Log.in.View + poly(Height_cm, 2, raw = T) + log10(Num_Stems) +
      offset(log(Deploy.Duration)),
    data = grsqDataWin
  )

ModListWinGSq[[18]] <- glm.nb.WlogStemsOak <-
  glm.nb(
    nSeqs ~ Log.in.View + OakDBH + log10(Num_Stems) +
      offset(log(Deploy.Duration)),
    data = grsqDataWin
  )

ModListWinGSq[[19]] <- glm.nb.WlogStemsEdd <-
  glm.nb(
    nSeqs ~ Log.in.View * Squirrel_EDD_WSp + log10(Num_Stems) +
      offset(log(Deploy.Duration)),
    data = grsqDataWin
  )

ModListWinGSq[[20]] <- glm.nb.WHgtStemsEdd <-
  glm.nb(
    nSeqs ~ poly(Height_cm, 2, raw = T) + Squirrel_EDD_WSp + log10(Num_Stems) +
      offset(log(Deploy.Duration)),
    data = grsqDataWin
  )

ModListWinGSq[[21]] <- glm.nb.WHgtStemsOak <-
  glm.nb(
    nSeqs ~ poly(Height_cm, 2, raw = T) + OakDBH + log10(Num_Stems) +
      offset(log(Deploy.Duration)),
    data = grsqDataWin
  )

ModListWinGSq[[22]] <- glm.nb.WHgtEddOak <-
  glm.nb(
    nSeqs ~ poly(Height_cm, 2, raw = T) + OakDBH + Squirrel_EDD_WSp +
      offset(log(Deploy.Duration)),
    data = grsqDataWin
  )
ModListWinGSq[[23]] <- glm.nb.WHgtEddLog <-
  glm.nb(
    nSeqs ~ poly(Height_cm, 2, raw = T) + Log.in.View * Squirrel_EDD_WSp +
      offset(log(Deploy.Duration)),
    data = grsqDataWin
  )
ModListWinGSq[[24]] <- glm.nb.WHgtOakLog <-
  glm.nb(
    nSeqs ~ poly(Height_cm, 2, raw = T) + Log.in.View + OakDBH +
      offset(log(Deploy.Duration)),
    data = grsqDataWin
  )
ModListWinGSq[[25]] <- glm.nb.WOakStemsEdd <-
  glm.nb(
    nSeqs ~ log10(Num_Stems) + Squirrel_EDD_WSp + OakDBH +
      offset(log(Deploy.Duration)),
    data = grsqDataWin
  )

ModListWinGSq[[26]] <- glm.nb.WLogOakEdd <-
  glm.nb(
    nSeqs ~ Log.in.View * Squirrel_EDD_WSp + OakDBH +
      offset(log(Deploy.Duration)),
    data = grsqDataWin
  )
modnamesW <- c("Intercept","Log","Height2","Oak","Stems","EDD","Log + Stems","Log * EDD","Log + Oak","Log + Height2","EDD + Height2","EDD + Stems","EDD + Oak","Stems + Height2","Stems + Oak","Height2 + Oak","Log + Stems + Height2","Log + Stems + Oak","Log * EDD + Stems","Height2 + Stems + EDD","Height2 + Stems + Oak","Height2 + EDD + Oak","Height2 + EDD * Log","Height2 + Oak + Log", "Stems + EDD + Oak", "Log * EDD + Oak")
modtabWinGSq <- aictab(cand.set = ModListWinGSq, modnames = modnamesW)
modtabWinGSq

summary(glm.nb.Whgt) #best model is height only model

#Summed Model Weights for Squirrels in Winter from full table. 
WinSW_edd <- sum(modtabWinGSq$AICcWt[grep("EDD", modtabWinGSq$Modnames)])
WinSW_OakDBH <- sum(modtabWinGSq$AICcWt[grep("Oak", modtabWinGSq$Modnames)])
WinSW_Stems <- sum(modtabWinGSq$AICcWt[grep("Stems", modtabWinGSq$Modnames)])
WinSW_Log <- sum(modtabWinGSq$AICcWt[grep("Log", modtabWinGSq$Modnames)])
WinSW_Height <- sum(modtabWinGSq$AICcWt[grep("Height", modtabWinGSq$Modnames)])



#Explained deviance of best model
1 - glm.nb.WHgt2$deviance / glm.nb.WHgt2$null.deviance #0.2348909
#Explained deviance of full model
1 - glm.nb.fullW$deviance / glm.nb.fullW$null.deviance #0.3688571


#Response curve of height. Could be interesting figure for paper
seq.WHgt <- seq(min(grsqDataWin$Height_cm), max(grsqDataWin$Height_cm), length.out = 100)

nd.seq.WHgt <- data.frame(Deploy.Duration = 61, Height_cm = seq.WHgt)

pred.WHgt <- predict(glm.nb.Whgt, newdata = nd.seq.WHgt, se=TRUE, type = "link")


plot(nSeqs ~ Height_cm, data=grsqDataWin, pch=16, las=1, xlab = "Camera Height (cm)", ylab = "# Gray Squirrel Sequences in Winter")
lines(nd.seq.WHgt$Height_cm, exp(pred.WHgt$fit), col="red")
lines(nd.seq.WHgt$Height_cm, exp(pred.WHgt$fit + 1.96*pred.WHgt$se.fit), lty=2, col="red")
lines(nd.seq.WHgt$Height_cm, exp(pred.WHgt$fit - 1.96*pred.WHgt$se.fit), lty=2, col="red")

#Making plot in ggplot instead

predPlotGrayWinHgt <- data.frame(nd.seq.WHgt, pred.WHgt) %>%
  mutate(
    lcl = fit - 2*se.fit,
    ucl = fit + 2*se.fit
  ) %>%
  mutate(
    response = exp(fit),
    lcl.response = exp(lcl),
    ucl.response = exp(ucl))

#Plotting commands in ggplot. We can mess around with this further for the paper as needed. The dots represent actual data points.
load(file = "data/responseTheme.Rdata")

GrSqWinHgtPlot <- ggplot(predPlotGrayWinHgt, aes(x=Height_cm, y=response)) +
  # Confidence region
  geom_ribbon(aes(ymin=lcl.response, ymax=ucl.response), alpha=0.25, show.legend = F) +
  # Prediction Lines
  geom_line() +
  geom_point(data = grsqDataWin, aes(x=Height_cm, y=nSeqs), show.legend = F, size = 2.5) +
  xlab("Camera height (cm)") +   
  ylab("Count of sequences") +
  geom_text(aes(x = 34, y = 50, label = "A"), size = 8)+
  myTheme

#ggsave("results/GrSqWinHgt.tiff", width = 6.5, height = 4.0, units = "in" )



# Spring Squirrel NB Full Model Selection ---------------------------------
#First test to see if interaction between log and EDD improves that two variable model

glm.nb.SPlog_Edd <-
  glm.nb(nSeqs ~ Log.in.View + Squirrel_EDD_WSp + offset(log(Deploy.Duration)), data = grsqDataSpr)
AICc(glm.nb.SPlog_Edd)

glm.nb.SPlogXEdd <-
  glm.nb(nSeqs ~ Log.in.View * Squirrel_EDD_WSp + offset(log(Deploy.Duration)), data = grsqDataSpr)
AICc(glm.nb.SPlogXEdd) #There is support for interaction, so this has to be added to each model where these two appear together.


#Compare height alone to height quadratic
glm.nb.SPHgt2 <-
  glm.nb(
    nSeqs ~ poly(Height_cm, 2, raw = T) +
      offset(log(Deploy.Duration)),
    data = grsqDataSpr
  )

glm.nb.SPhgt <-
  glm.nb(nSeqs ~ Height_cm + offset(log(Deploy.Duration)), data = grsqDataSpr)

AICc(glm.nb.SPHgt2)
AICc(glm.nb.SPhgt) #there is no support for a quadratic effect

#Full Spring Gray Squirrel Model Selection

glm.nb.fullSP <-
  glm.nb(
    nSeqs ~ Log.in.View * Squirrel_EDD_WSp + Height_cm + log10(Num_Stems) + OakDBH +
      offset(log(Deploy.Duration)),
    data = grsqDataSpr
  )

ModListSprGSq <- list(NA)


ModListSprGSq[[1]] <- glm.nb.SPI <- glm.nb(nSeqs ~ 1 + offset(log(Deploy.Duration)), data = grsqDataSpr)

ModListSprGSq[[2]] <- glm.nb.SPlog <-
  glm.nb(nSeqs ~ Log.in.View + offset(log(Deploy.Duration)), data = grsqDataSpr)

ModListSprGSq[[3]] <- glm.nb.SPhgt <-
  glm.nb(nSeqs ~ Height_cm + offset(log(Deploy.Duration)), data = grsqDataSpr)

ModListSprGSq[[4]] <- glm.nb.SPoak <-
  glm.nb(nSeqs ~ OakDBH + offset(log(Deploy.Duration)), data = grsqDataSpr)

ModListSprGSq[[5]] <- glm.nb.SPstems <-
  glm.nb(nSeqs ~ log10(Num_Stems) + offset(log(Deploy.Duration)), data = grsqDataSpr)

ModListSprGSq[[6]] <- glm.nb.SPedd <-
  glm.nb(nSeqs ~ Squirrel_EDD_WSp + offset(log(Deploy.Duration)), data = grsqDataSpr)

ModListSprGSq[[7]] <- glm.nb.SPlog_stems <-
  glm.nb(nSeqs ~ Log.in.View + log10(Num_Stems) + offset(log(Deploy.Duration)), data = grsqDataSpr)

ModListSprGSq[[8]] <-glm.nb.SPlog_Edd <-
  glm.nb(nSeqs ~ Log.in.View * Squirrel_EDD_WSp + offset(log(Deploy.Duration)), data = grsqDataSpr)

ModListSprGSq[[9]] <- glm.nb.SPlog_Oak <-
  glm.nb(nSeqs ~ Log.in.View + OakDBH + offset(log(Deploy.Duration)), data = grsqDataSpr)

ModListSprGSq[[10]] <- glm.nb.SPlog_Hgt <-
  glm.nb(nSeqs ~ Log.in.View + Height_cm + offset(log(Deploy.Duration)), data = grsqDataSpr)

ModListSprGSq[[11]] <- glm.nb.SPEdd_Hgt <-
  glm.nb(nSeqs ~ Squirrel_EDD_WSp + Height_cm + offset(log(Deploy.Duration)), data = grsqDataSpr)

ModListSprGSq[[12]] <-  glm.nb.SPEdd_stems <-
  glm.nb(nSeqs ~ Squirrel_EDD_WSp + log10(Num_Stems) + offset(log(Deploy.Duration)), data = grsqDataSpr)

ModListSprGSq[[13]] <- glm.nb.SPEdd_oak <-
  glm.nb(nSeqs ~ Squirrel_EDD_WSp + OakDBH + offset(log(Deploy.Duration)), data = grsqDataSpr)

ModListSprGSq[[14]] <- glm.nb.SPstems_hgt <-
  glm.nb(nSeqs ~ Height_cm + log10(Num_Stems) + offset(log(Deploy.Duration)), data = grsqDataSpr)

ModListSprGSq[[15]] <- glm.nb.SPstems_oak <-
  glm.nb(nSeqs ~ OakDBH + log10(Num_Stems) + offset(log(Deploy.Duration)), data = grsqDataSpr)

ModListSprGSq[[16]] <- glm.nb.SPhgt_oak <-
  glm.nb(nSeqs ~ Height_cm + OakDBH + offset(log(Deploy.Duration)), data = grsqDataSpr)

ModListSprGSq[[17]] <- glm.nb.SPlogStemsHgt <-
  glm.nb(
    nSeqs ~ Log.in.View + Height_cm + log10(Num_Stems) +
      offset(log(Deploy.Duration)),
    data = grsqDataSpr
  )

ModListSprGSq[[18]] <- glm.nb.SPlogStemsOak <-
  glm.nb(
    nSeqs ~ Log.in.View + OakDBH + log10(Num_Stems) +
      offset(log(Deploy.Duration)),
    data = grsqDataSpr
  )

ModListSprGSq[[19]] <- glm.nb.SPlogStemsEdd <-
  glm.nb(
    nSeqs ~ Log.in.View * Squirrel_EDD_WSp + log10(Num_Stems) +
      offset(log(Deploy.Duration)),
    data = grsqDataSpr
  )

ModListSprGSq[[20]] <- glm.nb.SPHgtStemsEdd <-
  glm.nb(
    nSeqs ~ Height_cm + Squirrel_EDD_WSp + log10(Num_Stems) +
      offset(log(Deploy.Duration)),
    data = grsqDataSpr
  )

ModListSprGSq[[21]] <- glm.nb.SPHgtStemsOak <-
  glm.nb(
    nSeqs ~ Height_cm + OakDBH + log10(Num_Stems) +
      offset(log(Deploy.Duration)),
    data = grsqDataSpr
  )

ModListSprGSq[[22]] <- glm.nb.SPHgtEddOak <-
  glm.nb(
    nSeqs ~ Height_cm + OakDBH + Squirrel_EDD_WSp +
      offset(log(Deploy.Duration)),
    data = grsqDataSpr
  )
ModListSprGSq[[23]] <- glm.nb.SPHgtEddLog <-
  glm.nb(
    nSeqs ~ Height_cm + Log.in.View * Squirrel_EDD_WSp +
      offset(log(Deploy.Duration)),
    data = grsqDataSpr
  )
ModListSprGSq[[24]] <- glm.nb.SPHgtOakLog <-
  glm.nb(
    nSeqs ~ Height_cm + Log.in.View + OakDBH +
      offset(log(Deploy.Duration)),
    data = grsqDataSpr
  )
ModListSprGSq[[25]] <- glm.nb.SPOakStemsEdd <-
  glm.nb(
    nSeqs ~ log10(Num_Stems) + Squirrel_EDD_WSp + OakDBH +
      offset(log(Deploy.Duration)),
    data = grsqDataSpr
  )

ModListSprGSq[[26]] <- glm.nb.SPLogOakEdd <-
  glm.nb(
    nSeqs ~ Log.in.View * Squirrel_EDD_WSp + OakDBH +
      offset(log(Deploy.Duration)),
    data = grsqDataSpr
  )
modnamesSP <- c("Intercept","Log","Height","Oak","Stems","EDD","Log + Stems","Log * EDD","Log + Oak","Log + Height","EDD + Height","EDD + Stems","EDD + Oak","Stems + Height","Stems + Oak","Height + Oak","Log + Stems + Height","Log + Stems + Oak","Log * EDD + Stems","Height + Stems + EDD","Height + Stems + Oak","Height + EDD + Oak","Height + EDD * Log","Height + Oak + Log", "Stems + EDD + Oak", "Log * EDD + Oak")
modtabSprGSq <- aictab(cand.set = ModListSprGSq, modnames = modnamesSP)
modtabSprGSq

summary(glm.nb.SPlogXEdd) #best model has log * EDD

#Summed Model Weights for Gray Squirrels in Spring from full table. 
SprSW_edd <- sum(modtabSprGSq$AICcWt[grep("EDD", modtabSprGSq$Modnames)])
SprSW_OakDBH <- sum(modtabSprGSq$AICcWt[grep("Oak", modtabSprGSq$Modnames)])
SprSW_Stems <- sum(modtabSprGSq$AICcWt[grep("Stems", modtabSprGSq$Modnames)])
SprSW_Log <- sum(modtabSprGSq$AICcWt[grep("Log", modtabSprGSq$Modnames)])
SprSW_Height <- sum(modtabSprGSq$AICcWt[grep("Height", modtabSprGSq$Modnames)])



#Explained deviance of best model
1 - glm.nb.SPlog_Edd$deviance / glm.nb.SPlog_Edd$null.deviance #0.3786397
#Explained deviance of full model
1 - glm.nb.fullSP$deviance / glm.nb.fullSP$null.deviance #0.4173307

#Response curve of EDD with and without logs. Could be a good one for the paper.

seq.SPEdd <- seq(min(grsqDataSpr$Squirrel_EDD_WSp), max(grsqDataSpr$Squirrel_EDD_WSp), length.out = 100)

nd.seq.SPEdd <- data.frame(Deploy.Duration = 61, Squirrel_EDD_WSp = seq.SPEdd, Log.in.View = "YES")
nd.seq.SPEddNo <- data.frame(Deploy.Duration = 61, Squirrel_EDD_WSp = seq.SPEdd, Log.in.View = "NO")

pred.SPEdd <- predict(glm.nb.SPlog_Edd, newdata = nd.seq.SPEdd, se=TRUE, type = "link")
pred.SPEddNo <- predict(glm.nb.SPlog_Edd, newdata = nd.seq.SPEddNo, se=TRUE, type = "link")

plot(nSeqs ~ Squirrel_EDD_WSp, data=grsqDataSpr, pch=16, col=grsqDataSpr$Log.in.View, las=1, xlab = "Effective Detection Distance (m)", ylab = "# Gray Squirrel Sequences in Spring")

lines(nd.seq.SPEdd$Squirrel_EDD_WSp, exp(pred.SPEdd$fit), col="red")
lines(nd.seq.SPEdd$Squirrel_EDD_WSp, exp(pred.SPEdd$fit + 1.96*pred.SPEdd$se.fit), lty=2, col="red")
lines(nd.seq.SPEdd$Squirrel_EDD_WSp, exp(pred.SPEdd$fit - 1.96*pred.SPEdd$se.fit), lty=2, col="red")

lines(nd.seq.SPEdd$Squirrel_EDD_WSp, exp(pred.SPEddNo$fit), col="black")
lines(nd.seq.SPEdd$Squirrel_EDD_WSp, exp(pred.SPEddNo$fit + 1.96*pred.SPEddNo$se.fit), lty=2, col="black")
lines(nd.seq.SPEdd$Squirrel_EDD_WSp, exp(pred.SPEddNo$fit - 1.96*pred.SPEddNo$se.fit), lty=2, col="black")

