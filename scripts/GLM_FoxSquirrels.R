###Regression Analysis of SCBI Camera Grid Data

rm(list = ls())

require(AICcmodavg)
require(MASS) #for glm.nb function
source("scripts/pairsPannelFunctions.r")


# Data Import and Exploration - Fox Squirrels ---------------------------------------

load("data/foxsqDataSum.RData")
load("data/foxsqDataFall.RData")
load("data/foxsqDataWin.RData")
load("data/foxsqDataSpr.RData")

pairs(foxsqDataSum[,c(4,7,13,14,15, 17)], diag.panel = panel.hist, lower.panel = panel.smooth, upper.panel = panel.cor)
str(foxsqDataSum)

hist(sqDataSum$Num_Stems)
hist(log10(sqDataSum$Num_Stems)) #because one site has so many stems, we need to log10 transform this parameter for all models
hist(sqDataSum$Squirrel_EDD_4S)

#Univariate exploratory plots for each of our 5 covariates
par.default <- par(no.readonly = T)
par(mfrow = c(2,2))
plot(nSeqs ~ Squirrel_EDD_4S, data = foxsqDataSum, main = "Summer")
plot(nSeqs ~ Squirrel_EDD_4S, data = foxsqDataFall, main = "Fall")
plot(nSeqs ~ Squirrel_EDD_WSp, data = foxsqDataWin, main = "Winter")
plot(nSeqs ~ Squirrel_EDD_WSp, data = foxsqDataSpr, main = "Spring")

plot(nSeqs ~ OakDBH, data = foxsqDataSum, main = "Summer")
plot(nSeqs ~ OakDBH, data = foxsqDataFall, main = "Fall")
plot(nSeqs ~ OakDBH, data = foxsqDataWin, main = "Winter")
plot(nSeqs ~ OakDBH, data = foxsqDataSpr, main = "Spring")

plot(nSeqs ~ log10(Num_Stems), data = foxsqDataSum, main = "Summer")
plot(nSeqs ~ log10(Num_Stems), data = foxsqDataFall, main = "Fall")
plot(nSeqs ~ log10(Num_Stems), data = foxsqDataWin, main = "Winter")
plot(nSeqs ~ log10(Num_Stems), data = foxsqDataSpr, main = "Spring")

#These seems to show a quadratic relationship actually, which might make sense
plot(nSeqs ~ Height_cm, data = foxsqDataSum, main = "Summer")
plot(nSeqs ~ Height_cm, data = foxsqDataFall, main = "Fall")
plot(nSeqs ~ Height_cm, data = foxsqDataWin, main = "Winter")
plot(nSeqs ~ Height_cm, data = foxsqDataSpr, main = "Spring")

plot(nSeqs ~ Log.in.View, data = foxsqDataSum, main = "Summer")
plot(nSeqs ~ Log.in.View, data = foxsqDataFall, main = "Fall")
plot(nSeqs ~ Log.in.View, data = foxsqDataWin, main = "Winter")
plot(nSeqs ~ Log.in.View, data = foxsqDataSpr, main = "Spring")
par(par.default)

hist(foxsqDataSum$Deploy.Duration, breaks = 5)

#Looking at how many cameras picked up bear in each season
barplot(foxsqDataSum$nSeqs, xlab = "Summer", ylab = "# of Seqs per camera")
barplot(foxsqDataFall$nSeqs, xlab = "Fall", ylab = "# of Seqs per camera")
barplot(foxsqDataWin$nSeqs, xlab = "Winter", ylab = "# of Seqs per camera")
barplot(foxsqDataSpr$nSeqs, xlab = "Spring", ylab = "# of Seqs per camera")


# Full Poisson GLM and Overdispersion -------------------------------------

#Exploratory Poisson model for Summer data to see if a negative binomial model is needed due to potential overdispersion
glm.po.fullS <- glm(nSeqs ~ Log.in.View + Squirrel_EDD_4S + Height_cm + log10(Num_Stems) + OakDBH +offset(log(Deploy.Duration)), data = foxsqDataSum, family = poisson)
summary(glm.po.fullS)

#test for overdispersion indicates we should be using negative binomial (2.7)
sum.po.fullS <- summary(glm.po.fullS)
phi.po.fullS <- sum.po.fullS$deviance/sum.po.fullS$df.residual 

#Exploratory Poisson model for Fall data
glm.po.fullF <- glm(nSeqs ~ Log.in.View + Squirrel_EDD_4S + Height_cm + log10(Num_Stems) + OakDBH +offset(log(Deploy.Duration)), data = foxsqDataFall, family = poisson)
summary(glm.po.fullF)

#test for overdispersion indicates we should be using negative binomial (3.7) 
sum.po.fullF <- summary(glm.po.fullF)
phi.po.fullF <- sum.po.fullF$deviance/sum.po.fullF$df.residual 

#Exploratory Poisson model for Winter data
glm.po.fullW <- glm(nSeqs ~ Log.in.View + Squirrel_EDD_WSp + Height_cm + log10(Num_Stems) + OakDBH +offset(log(Deploy.Duration)), data = foxsqDataWin, family = poisson)
summary(glm.po.fullW)

#test for overdispersion indicates we should be using negative binomial  (12.1)
sum.po.fullW <- summary(glm.po.fullW)
phi.po.fullW <- sum.po.fullW$deviance/sum.po.fullW$df.residual 

#Exploratory Poisson model for Spring data
glm.po.fullSp <- glm(nSeqs ~ Log.in.View + Squirrel_EDD_WSp + Height_cm + log10(Num_Stems) + OakDBH +offset(log(Deploy.Duration)), data = foxsqDataSpr, family = poisson)
summary(glm.po.fullSp)

#test for overdispersion indicates we should be using negative binomial (8.5)
sum.po.fullSp <- summary(glm.po.fullSp)
phi.po.fullSp <- sum.po.fullSp$deviance/sum.po.fullSp$df.residual 


# Negative Binomial Regression - Fox Squirrels in Summer --------------------------------

#First test to see if interaction between log and EDD improves that two variable model

glm.nb.Slog_Edd <-
  glm.nb(nSeqs ~ Log.in.View + Squirrel_EDD_4S + offset(log(Deploy.Duration)), data = foxsqDataSum)
AICc(glm.nb.Slog_Edd)

glm.nb.SlogXEdd <-
  glm.nb(nSeqs ~ Log.in.View * Squirrel_EDD_4S + offset(log(Deploy.Duration)), data = foxsqDataSum)
AICc(glm.nb.SlogXEdd) #No support for interaction 


#Compare height alone to height quadratic
glm.nb.SHgt2 <-
  glm.nb(
    nSeqs ~ poly(Height_cm, 2, raw = T) +
      offset(log(Deploy.Duration)),
    data = foxsqDataSum
  )

glm.nb.Shgt <-
  glm.nb(
    nSeqs ~ Height_cm +
      offset(log(Deploy.Duration)),
    data = foxsqDataSum
  )

AICc(glm.nb.SHgt2)
AICc(glm.nb.Shgt) #no support for quadratic effect

#Full Model Selection for Fox Squirrels in Summer

ModListSumFSq <- list(NA)

ModListSumFSq[[1]] <- glm.nb.fullS <-
  glm.nb(
    nSeqs ~ Log.in.View + Squirrel_EDD_4S + Height_cm + log10(Num_Stems) + OakDBH +
      offset(log(Deploy.Duration)),
    data = foxsqDataSum
  )

ModListSumFSq[[2]] <- glm.nb.SI <- glm.nb(nSeqs ~ 1 + offset(log(Deploy.Duration)), data = foxsqDataSum)

ModListSumFSq[[3]] <- glm.nb.Slog <-
  glm.nb(nSeqs ~ Log.in.View + offset(log(Deploy.Duration)), data = foxsqDataSum)

ModListSumFSq[[4]] <- glm.nb.Shgt <-
  glm.nb(nSeqs ~ Height_cm + offset(log(Deploy.Duration)), data = foxsqDataSum)

ModListSumFSq[[5]] <- glm.nb.Soak <-
  glm.nb(nSeqs ~ OakDBH + offset(log(Deploy.Duration)), data = foxsqDataSum)

ModListSumFSq[[6]] <- glm.nb.Sstems <-
  glm.nb(nSeqs ~ log10(Num_Stems) + offset(log(Deploy.Duration)), data = foxsqDataSum)

ModListSumFSq[[7]] <- glm.nb.Sedd <-
  glm.nb(nSeqs ~ Squirrel_EDD_4S + offset(log(Deploy.Duration)), data = foxsqDataSum)

ModListSumFSq[[8]] <- glm.nb.Slog_stems <-
  glm.nb(nSeqs ~ Log.in.View + log10(Num_Stems) + offset(log(Deploy.Duration)), data = foxsqDataSum)

ModListSumFSq[[9]] <-glm.nb.Slog_Edd <-
  glm.nb(nSeqs ~ Log.in.View + Squirrel_EDD_4S + offset(log(Deploy.Duration)), data = foxsqDataSum)

ModListSumFSq[[10]] <- glm.nb.Slog_Oak <-
  glm.nb(nSeqs ~ Log.in.View + OakDBH + offset(log(Deploy.Duration)), data = foxsqDataSum)

ModListSumFSq[[11]] <- glm.nb.Slog_Hgt <-
  glm.nb(nSeqs ~ Log.in.View + Height_cm + offset(log(Deploy.Duration)), data = foxsqDataSum)

ModListSumFSq[[12]] <- glm.nb.SEdd_Hgt <-
  glm.nb(nSeqs ~ Squirrel_EDD_4S + Height_cm + offset(log(Deploy.Duration)), data = foxsqDataSum)

ModListSumFSq[[13]] <-  glm.nb.SEdd_stems <-
  glm.nb(nSeqs ~ Squirrel_EDD_4S + log10(Num_Stems) + offset(log(Deploy.Duration)), data = foxsqDataSum)

ModListSumFSq[[14]] <- glm.nb.SEdd_oak <-
  glm.nb(nSeqs ~ Squirrel_EDD_4S + OakDBH + offset(log(Deploy.Duration)), data = foxsqDataSum)

ModListSumFSq[[15]] <- glm.nb.Sstems_hgt <-
  glm.nb(nSeqs ~ Height_cm + log10(Num_Stems) + offset(log(Deploy.Duration)), data = foxsqDataSum)

ModListSumFSq[[16]] <- glm.nb.Sstems_oak <-
  glm.nb(nSeqs ~ OakDBH + log10(Num_Stems) + offset(log(Deploy.Duration)), data = foxsqDataSum)

ModListSumFSq[[17]] <- glm.nb.Shgt_oak <-
  glm.nb(nSeqs ~ Height_cm + OakDBH + offset(log(Deploy.Duration)), data = foxsqDataSum)

ModListSumFSq[[18]] <- glm.nb.SlogStemsHgt <-
  glm.nb(
    nSeqs ~ Log.in.View + Height_cm + log10(Num_Stems) +
      offset(log(Deploy.Duration)),
    data = foxsqDataSum
  )

ModListSumFSq[[19]] <- glm.nb.SlogStemsOak <-
  glm.nb(
    nSeqs ~ Log.in.View + OakDBH + log10(Num_Stems) +
      offset(log(Deploy.Duration)),
    data = foxsqDataSum
  )

ModListSumFSq[[20]] <- glm.nb.SlogStemsEdd <-
  glm.nb(
    nSeqs ~ Log.in.View + Squirrel_EDD_4S + log10(Num_Stems) +
      offset(log(Deploy.Duration)),
    data = foxsqDataSum
  )

ModListSumFSq[[21]] <- glm.nb.SHgtStemsEdd <-
  glm.nb(
    nSeqs ~ Height_cm + Squirrel_EDD_4S + log10(Num_Stems) +
      offset(log(Deploy.Duration)),
    data = foxsqDataSum
  )

ModListSumFSq[[22]] <- glm.nb.SHgtStemsOak <-
  glm.nb(
    nSeqs ~ Height_cm + OakDBH + log10(Num_Stems) +
      offset(log(Deploy.Duration)),
    data = foxsqDataSum
  )

ModListSumFSq[[23]] <- glm.nb.SHgtEddOak <-
  glm.nb(
    nSeqs ~ Height_cm + OakDBH + Squirrel_EDD_4S +
      offset(log(Deploy.Duration)),
    data = foxsqDataSum
  )
ModListSumFSq[[24]] <- glm.nb.SHgtEddLog <-
  glm.nb(
    nSeqs ~ Height_cm + Log.in.View + Squirrel_EDD_4S +
      offset(log(Deploy.Duration)),
    data = foxsqDataSum
  )
ModListSumFSq[[25]] <- glm.nb.SHgtOakLog <-
  glm.nb(
    nSeqs ~ Height_cm + Log.in.View + OakDBH +
      offset(log(Deploy.Duration)),
    data = foxsqDataSum
  )
modnamesS <- c("Full","Intercept","Log","Height","Oak","Stems","EDD","Log + Stems","Log + EDD","Log + Oak","Log + Height","EDD + Height","EDD + Stems","EDD + Oak","Stems + Height","Stems + Oak","Height + Oak","Log + Stems + Height","Log + Stems + Oak","Log + Stems + EDD","Height + Stems + EDD","Height + Stems + Oak","Height + EDD + Oak","Height + EDD + Log","Height + Oak + Log")
modtabSumFSq <- aictab(cand.set = ModListSumFSq, modnames = modnamesS)
modtabSumFSq

#Best model has height, EDD and log
summary(glm.nb.SHgtEddLog)
AICc(glm.nb.SHgtEddLog)

#Summed Model Weights for Squirrels in Summer from full table. 
SumSW_edd <- sum(modtabSumFSq$AICcWt[grep("EDD", modtabSumFSq$Modnames)])
SumSW_OakDBH <- sum(modtabSumFSq$AICcWt[grep("Oak", modtabSumFSq$Modnames)])
SumSW_Stems <- sum(modtabSumFSq$AICcWt[grep("Stems", modtabSumFSq$Modnames)])
SumSW_Log <- sum(modtabSumFSq$AICcWt[grep("Log", modtabSumFSq$Modnames)])
SumSW_Height <- sum(modtabSumFSq$AICcWt[grep("Height", modtabSumFSq$Modnames)])


#Explained deviance of best model
1 - glm.nb.SHgtEddLog$deviance / glm.nb.SHgtEddLog$null.deviance #0.4888731
#Explained deviance of full model
1 - glm.nb.fullS$deviance / glm.nb.fullS$null.deviance #0.5444852

#Response curves here for EDD with and without logs

seq.SEdd <- seq(min(foxsqDataSum$Squirrel_EDD_4S), max(foxsqDataSum$Squirrel_EDD_4S), length.out = 100)

nd.seq.SEdd <- data.frame(Deploy.Duration = 61, Squirrel_EDD_4S = seq.SEdd, Log.in.View = "YES", Height_cm = median(foxsqDataSum$Height_cm))
nd.seq.SEddNo <- data.frame(Deploy.Duration = 61, Squirrel_EDD_4S = seq.SEdd, Log.in.View = "NO", Height_cm = median(foxsqDataSum$Height_cm))

pred.SEdd <- predict(glm.nb.SHgtEddLog, newdata = nd.seq.SEdd, se=TRUE, type = "link")
pred.SEddNo <- predict(glm.nb.SHgtEddLog, newdata = nd.seq.SEddNo, se=TRUE, type = "link")

plot(nSeqs ~ Squirrel_EDD_4S, data=foxsqDataSum, pch=16, col=foxsqDataSum$Log.in.View, las=1, xlab = "Effective Detection Distance (m)", ylab = "# Fox Squirrel Sequences in Summer")
lines(nd.seq.SEdd$Squirrel_EDD_4S, exp(pred.SEdd$fit), col="red")
lines(nd.seq.SEdd$Squirrel_EDD_4S, exp(pred.SEdd$fit + 1.96*pred.SEdd$se.fit), lty=2, col="red")
lines(nd.seq.SEdd$Squirrel_EDD_4S, exp(pred.SEdd$fit - 1.96*pred.SEdd$se.fit), lty=2, col="red")

lines(nd.seq.SEdd$Squirrel_EDD_4S, exp(pred.SEddNo$fit), col="black")
lines(nd.seq.SEdd$Squirrel_EDD_4S, exp(pred.SEddNo$fit + 1.96*pred.SEddNo$se.fit), lty=2, col="black")
lines(nd.seq.SEdd$Squirrel_EDD_4S, exp(pred.SEddNo$fit - 1.96*pred.SEddNo$se.fit), lty=2, col="black")

#Response curves of height, with and without logs
seq.SHgt <- seq(min(foxsqDataSum$Height_cm), max(foxsqDataSum$Height_cm), length.out = 100)

nd.seq.SHgt <- data.frame(Deploy.Duration = 61, Squirrel_EDD_4S = median(foxsqDataSum$Squirrel_EDD_4S), Log.in.View = "YES", Height_cm = seq.SHgt)
nd.seq.SHgtNo <- data.frame(Deploy.Duration = 61, Squirrel_EDD_4S = median(foxsqDataSum$Squirrel_EDD_4S), Log.in.View = "NO", Height_cm = seq.SHgt)

pred.SHgtY <- predict(glm.nb.SHgtEddLog, newdata = nd.seq.SHgt, se=TRUE, type = "link")
pred.SHgtN <- predict(glm.nb.SHgtEddLog, newdata = nd.seq.SHgtNo, se=TRUE, type = "link")

plot(nSeqs ~ Height_cm, data=foxsqDataSum, pch=16, col=foxsqDataSum$Log.in.View, las=1, xlab = "Camera Height (cm)", ylab = "# Squirrel Sequences in Summer")
lines(nd.seq.SHgt$Height_cm, exp(pred.SHgtY$fit), col="red")
lines(nd.seq.SHgt$Height_cm, exp(pred.SHgtY$fit + 1.96*pred.SHgtY$se.fit), lty=2, col="red")
lines(nd.seq.SHgt$Height_cm, exp(pred.SHgtY$fit - 1.96*pred.SHgtY$se.fit), lty=2, col="red")

lines(nd.seq.SHgt$Height_cm, exp(pred.SHgtN$fit), col="black")
lines(nd.seq.SHgt$Height_cm, exp(pred.SHgtN$fit + 1.96*pred.SHgtN$se.fit), lty=2, col="black")
lines(nd.seq.SHgt$Height_cm, exp(pred.SHgtN$fit - 1.96*pred.SHgtN$se.fit), lty=2, col="black")


# Fox Squirrel in Fall NB Full Model Selection --------------------------------

#First test to see if interaction between log and EDD improves that two variable model

glm.nb.Flog_Edd <-
  glm.nb(nSeqs ~ Log.in.View + Squirrel_EDD_4S + offset(log(Deploy.Duration)), data = foxsqDataFall)
AICc(glm.nb.Flog_Edd)

glm.nb.FlogXEdd <-
  glm.nb(nSeqs ~ Log.in.View * Squirrel_EDD_4S + offset(log(Deploy.Duration)), data = foxsqDataFall)
AICc(glm.nb.FlogXEdd) #There is no support for interaction.


#Compare height alone to height quadratic
glm.nb.FHgt2 <-
  glm.nb(
    nSeqs ~ poly(Height_cm, 2, raw = T) +
      offset(log(Deploy.Duration)),
    data = foxsqDataFall
  )

glm.nb.Fhgt <-
  glm.nb(nSeqs ~ Height_cm + offset(log(Deploy.Duration)), data = foxsqDataFall)

AICc(glm.nb.FHgt2)
AICc(glm.nb.Fhgt) #no support for quadratic effect

#Full Model Selection for Fall Fox Squirrels

ModListFallFSq <- list(NA)

ModListFallFSq[[1]] <- glm.nb.fullF <-
  glm.nb(
    nSeqs ~ Log.in.View + Squirrel_EDD_4S + Height_cm + log10(Num_Stems) + OakDBH +
      offset(log(Deploy.Duration)),
    data = foxsqDataFall
  )

ModListFallFSq[[2]] <- glm.nb.FI <- glm.nb(nSeqs ~ 1 + offset(log(Deploy.Duration)), data = foxsqDataFall)

ModListFallFSq[[3]] <- glm.nb.Flog <-
  glm.nb(nSeqs ~ Log.in.View + offset(log(Deploy.Duration)), data = foxsqDataFall)

ModListFallFSq[[4]] <- glm.nb.Fhgt <-
  glm.nb(nSeqs ~ Height_cm + offset(log(Deploy.Duration)), data = foxsqDataFall)

ModListFallFSq[[5]] <- glm.nb.Foak <-
  glm.nb(nSeqs ~ OakDBH + offset(log(Deploy.Duration)), data = foxsqDataFall)

ModListFallFSq[[6]] <- glm.nb.Fstems <-
  glm.nb(nSeqs ~ log10(Num_Stems) + offset(log(Deploy.Duration)), data = foxsqDataFall)

ModListFallFSq[[7]] <- glm.nb.Fedd <-
  glm.nb(nSeqs ~ Squirrel_EDD_4S + offset(log(Deploy.Duration)), data = foxsqDataFall)

ModListFallFSq[[8]] <- glm.nb.Flog_stems <-
  glm.nb(nSeqs ~ Log.in.View + log10(Num_Stems) + offset(log(Deploy.Duration)), data = foxsqDataFall)

ModListFallFSq[[9]] <- glm.nb.Flog_Edd <-
  glm.nb(nSeqs ~ Log.in.View + Squirrel_EDD_4S + offset(log(Deploy.Duration)), data = foxsqDataFall)

ModListFallFSq[[10]] <- glm.nb.Flog_Oak <-
  glm.nb(nSeqs ~ Log.in.View + OakDBH + offset(log(Deploy.Duration)), data = foxsqDataFall)

ModListFallFSq[[11]] <- glm.nb.Flog_Hgt <-
  glm.nb(nSeqs ~ Log.in.View + Height_cm + offset(log(Deploy.Duration)), data = foxsqDataFall)

ModListFallFSq[[12]] <- glm.nb.FEdd_Hgt <-
  glm.nb(nSeqs ~ Squirrel_EDD_4S + Height_cm + offset(log(Deploy.Duration)), data = foxsqDataFall)

ModListFallFSq[[13]] <-  glm.nb.FEdd_stems <-
  glm.nb(nSeqs ~ Squirrel_EDD_4S + log10(Num_Stems) + offset(log(Deploy.Duration)), data = foxsqDataFall)

ModListFallFSq[[14]] <- glm.nb.FEdd_oak <-
  glm.nb(nSeqs ~ Squirrel_EDD_4S + OakDBH + offset(log(Deploy.Duration)), data = foxsqDataFall)

ModListFallFSq[[15]] <- glm.nb.Fstems_hgt <-
  glm.nb(nSeqs ~ Height_cm + log10(Num_Stems) + offset(log(Deploy.Duration)), data = foxsqDataFall)

ModListFallFSq[[16]] <- glm.nb.Fstems_oak <-
  glm.nb(nSeqs ~ OakDBH + log10(Num_Stems) + offset(log(Deploy.Duration)), data = foxsqDataFall)

ModListFallFSq[[17]] <- glm.nb.Fhgt_oak <-
  glm.nb(nSeqs ~ Height_cm + OakDBH + offset(log(Deploy.Duration)), data = foxsqDataFall)

ModListFallFSq[[18]] <- glm.nb.FlogStemsHgt <-
  glm.nb(
    nSeqs ~ Log.in.View + Height_cm + log10(Num_Stems) +
      offset(log(Deploy.Duration)),
    data = foxsqDataFall
  )

ModListFallFSq[[19]] <- glm.nb.FlogStemsOak <-
  glm.nb(
    nSeqs ~ Log.in.View + OakDBH + log10(Num_Stems) +
      offset(log(Deploy.Duration)),
    data = foxsqDataFall
  )

ModListFallFSq[[20]] <- glm.nb.FlogStemsEdd <-
  glm.nb(
    nSeqs ~ Log.in.View + Squirrel_EDD_4S + log10(Num_Stems) +
      offset(log(Deploy.Duration)),
    data = foxsqDataFall
  )

ModListFallFSq[[21]] <- glm.nb.FHgtStemsEdd <-
  glm.nb(
    nSeqs ~ Height_cm + Squirrel_EDD_4S + log10(Num_Stems) +
      offset(log(Deploy.Duration)),
    data = foxsqDataFall
  )

ModListFallFSq[[22]] <- glm.nb.FHgtStemsOak <-
  glm.nb(
    nSeqs ~ Height_cm + OakDBH + log10(Num_Stems) +
      offset(log(Deploy.Duration)),
    data = foxsqDataFall
  )

ModListFallFSq[[23]] <- glm.nb.FHgtEddOak <-
  glm.nb(
    nSeqs ~ Height_cm + OakDBH + Squirrel_EDD_4S +
      offset(log(Deploy.Duration)),
    data = foxsqDataFall
  )
ModListFallFSq[[24]] <- glm.nb.FHgtEddLog <-
  glm.nb(
    nSeqs ~ Height_cm + Log.in.View + Squirrel_EDD_4S +
      offset(log(Deploy.Duration)),
    data = foxsqDataFall
  )
ModListFallFSq[[25]] <- glm.nb.FHgtOakLog <-
  glm.nb(
    nSeqs ~ Height_cm + Log.in.View + OakDBH +
      offset(log(Deploy.Duration)),
    data = foxsqDataFall
  )
modnamesF <- c("Full","Intercept","Log","Height","Oak","Stems","EDD","Log + Stems","Log + EDD","Log + Oak","Log + Height","EDD + Height","EDD + Stems","EDD + Oak","Stems + Height","Stems + Oak","Height + Oak","Log + Stems + Height","Log + Stems + Oak","Log + Stems + EDD","Height + Stems + EDD","Height + Stems + Oak","Height + EDD + Oak","Height + EDD + Log","Height + Oak + Log")
modtabFallFSq <- aictab(cand.set = ModListFallFSq, modnames = modnamesF)
modtabFallFSq

#Best model is clear as log and Stems and EDD
summary(glm.nb.FlogStemsEdd)

#Summed Model Weights Squirrels in Fall from full table. Note interaction not included at the moment
FallSW_edd <- sum(modtabFallFSq$AICcWt[grep("EDD", modtabFallFSq$Modnames)])
FallSW_OakDBH <- sum(modtabFallFSq$AICcWt[grep("Oak", modtabFallFSq$Modnames)])
FallSW_Stems <- sum(modtabFallFSq$AICcWt[grep("Stems", modtabFallFSq$Modnames)])
FallSW_Log <- sum(modtabFallFSq$AICcWt[grep("Log", modtabFallFSq$Modnames)])
FallSW_Height <- sum(modtabFallFSq$AICcWt[grep("Height", modtabFallFSq$Modnames)])


#Explained Deviance of best model
1 - glm.nb.FlogStemsEdd$deviance / glm.nb.FlogStemsEdd$null.deviance #0.4446317
#Explained Deviance of best model
1 - glm.nb.fullF$deviance / glm.nb.fullF$null.deviance #0.4629069

#Response Curve for EDD impact, with and without log present
#NEEDS TO BE REDONE. Could be interesting to do one for Stems, since this may be the only model where it is important.
#seq.FEdd <- seq(min(sqDataFall$Squirrel_EDD_4S), max(sqDataFall$Squirrel_EDD_4S), length.out = 100)

#nd.seq.FEdd <- data.frame(Deploy.Duration = 61, Squirrel_EDD_4S = seq.FEdd, Log.in.View = "YES")
#nd.seq.FEddNo <- data.frame(Deploy.Duration = 61, Squirrel_EDD_4S = seq.FEdd, Log.in.View = "NO")

#pred.FEdd <- predict(glm.nb.FlogXEdd, newdata = nd.seq.FEdd, se=TRUE, type = "link")
#pred.FEddNo <- predict(glm.nb.FlogXEdd, newdata = nd.seq.FEddNo, se=TRUE, type = "link")

#plot(nSeqs ~ Squirrel_EDD_4S, data=sqDataFall, pch=16, col=sqDataFall$Log.in.View, las=1, xlab = "Effective Detection Distance (m)", ylab = "# Squirrel Sequences in Fall")
#lines(nd.seq.FEdd$Squirrel_EDD_4S, exp(pred.FEdd$fit), col="red")
#lines(nd.seq.FEdd$Squirrel_EDD_4S, exp(pred.FEdd$fit + 1.96*pred.FEdd$se.fit), lty=2, col="red")
#lines(nd.seq.FEdd$Squirrel_EDD_4S, exp(pred.FEdd$fit - 1.96*pred.FEdd$se.fit), lty=2, col="red")

#lines(nd.seq.FEdd$Squirrel_EDD_4S, exp(pred.FEddNo$fit), col="black")
#lines(nd.seq.FEdd$Squirrel_EDD_4S, exp(pred.FEddNo$fit + 1.96*pred.FEddNo$se.fit), lty=2, col="black")
#lines(nd.seq.FEdd$Squirrel_EDD_4S, exp(pred.FEddNo$fit - 1.96*pred.FEddNo$se.fit), lty=2, col="black")


#Winter Fox Squirrel NB Full Model Selection --------------------------------

#First test to see if interaction between log and EDD improves that two variable model

glm.nb.Wlog_Edd <-
  glm.nb(nSeqs ~ Log.in.View + Squirrel_EDD_WSp + offset(log(Deploy.Duration)), data = foxsqDataWin)
AICc(glm.nb.Wlog_Edd)

glm.nb.WlogXEdd <-
  glm.nb(nSeqs ~ Log.in.View * Squirrel_EDD_WSp + offset(log(Deploy.Duration)), data = foxsqDataWin)
AICc(glm.nb.WlogXEdd) #There is no support for interaction


#Compare height alone to height quadratic
glm.nb.WHgt2 <-
  glm.nb(
    nSeqs ~ poly(Height_cm, 2, raw = T) +
      offset(log(Deploy.Duration)),
    data = foxsqDataWin
  )

glm.nb.Whgt <-
  glm.nb(nSeqs ~ Height_cm + offset(log(Deploy.Duration)), data = foxsqDataWin)

AICc(glm.nb.WHgt2)
AICc(glm.nb.Whgt) #no support for quadratic effect

#Full Model Comparison Run Winter Fox Squirrels

ModListWinFSq <- list(NA)
ModListWinFSq[[1]] <- glm.nb.fullW <-
  glm.nb(
    nSeqs ~ Log.in.View + Squirrel_EDD_WSp + Height_cm + log10(Num_Stems) + OakDBH +
      offset(log(Deploy.Duration)),
    data = foxsqDataWin
  )

ModListWinFSq[[2]] <- glm.nb.WI <- glm.nb(nSeqs ~ 1 + offset(log(Deploy.Duration)), data = foxsqDataWin)

ModListWinFSq[[3]] <- glm.nb.Wlog <-
  glm.nb(nSeqs ~ Log.in.View + offset(log(Deploy.Duration)), data = foxsqDataWin)

ModListWinFSq[[4]] <- glm.nb.Whgt <-
  glm.nb(nSeqs ~ Height_cm + offset(log(Deploy.Duration)), data = foxsqDataWin)

ModListWinFSq[[5]] <- glm.nb.Woak <-
  glm.nb(nSeqs ~ OakDBH + offset(log(Deploy.Duration)), data = foxsqDataWin)

ModListWinFSq[[6]] <- glm.nb.Wstems <-
  glm.nb(nSeqs ~ log10(Num_Stems) + offset(log(Deploy.Duration)), data = foxsqDataWin)

ModListWinFSq[[7]] <- glm.nb.Wedd <-
  glm.nb(nSeqs ~ Squirrel_EDD_WSp + offset(log(Deploy.Duration)), data = foxsqDataWin)

ModListWinFSq[[8]] <- glm.nb.Wlog_stems <-
  glm.nb(nSeqs ~ Log.in.View + log10(Num_Stems) + offset(log(Deploy.Duration)), data = foxsqDataWin)

ModListWinFSq[[9]] <-glm.nb.Wlog_Edd <-
  glm.nb(nSeqs ~ Log.in.View + Squirrel_EDD_WSp + offset(log(Deploy.Duration)), data = foxsqDataWin)

ModListWinFSq[[10]] <- glm.nb.Wlog_Oak <-
  glm.nb(nSeqs ~ Log.in.View + OakDBH + offset(log(Deploy.Duration)), data = foxsqDataWin)

ModListWinFSq[[11]] <- glm.nb.Wlog_Hgt <-
  glm.nb(nSeqs ~ Log.in.View + Height_cm + offset(log(Deploy.Duration)), data = foxsqDataWin)

ModListWinFSq[[12]] <- glm.nb.WEdd_Hgt <-
  glm.nb(nSeqs ~ Squirrel_EDD_WSp + Height_cm + offset(log(Deploy.Duration)), data = foxsqDataWin)

ModListWinFSq[[13]] <-  glm.nb.WEdd_stems <-
  glm.nb(nSeqs ~ Squirrel_EDD_WSp + log10(Num_Stems) + offset(log(Deploy.Duration)), data = foxsqDataWin)

ModListWinFSq[[14]] <- glm.nb.WEdd_oak <-
  glm.nb(nSeqs ~ Squirrel_EDD_WSp + OakDBH + offset(log(Deploy.Duration)), data = foxsqDataWin)

ModListWinFSq[[15]] <- glm.nb.Wstems_hgt <-
  glm.nb(nSeqs ~ Height_cm + log10(Num_Stems) + offset(log(Deploy.Duration)), data = foxsqDataWin)

ModListWinFSq[[16]] <- glm.nb.Wstems_oak <-
  glm.nb(nSeqs ~ OakDBH + log10(Num_Stems) + offset(log(Deploy.Duration)), data = foxsqDataWin)

ModListWinFSq[[17]] <- glm.nb.Whgt_oak <-
  glm.nb(nSeqs ~ Height_cm + OakDBH + offset(log(Deploy.Duration)), data = foxsqDataWin)

ModListWinFSq[[18]] <- glm.nb.WlogStemsHgt <-
  glm.nb(
    nSeqs ~ Log.in.View + Height_cm + log10(Num_Stems) +
      offset(log(Deploy.Duration)),
    data = foxsqDataWin
  )

ModListWinFSq[[19]] <- glm.nb.WlogStemsOak <-
  glm.nb(
    nSeqs ~ Log.in.View + OakDBH + log10(Num_Stems) +
      offset(log(Deploy.Duration)),
    data = foxsqDataWin
  )

ModListWinFSq[[20]] <- glm.nb.WlogStemsEdd <-
  glm.nb(
    nSeqs ~ Log.in.View + Squirrel_EDD_WSp + log10(Num_Stems) +
      offset(log(Deploy.Duration)),
    data = foxsqDataWin
  )

ModListWinFSq[[21]] <- glm.nb.WHgtStemsEdd <-
  glm.nb(
    nSeqs ~ Height_cm + Squirrel_EDD_WSp + log10(Num_Stems) +
      offset(log(Deploy.Duration)),
    data = foxsqDataWin
  )

ModListWinFSq[[22]] <- glm.nb.WHgtStemsOak <-
  glm.nb(
    nSeqs ~ Height_cm + OakDBH + log10(Num_Stems) +
      offset(log(Deploy.Duration)),
    data = foxsqDataWin
  )

ModListWinFSq[[23]] <- glm.nb.WHgtEddOak <-
  glm.nb(
    nSeqs ~ Height_cm + OakDBH + Squirrel_EDD_WSp +
      offset(log(Deploy.Duration)),
    data = foxsqDataWin
  )
ModListWinFSq[[24]] <- glm.nb.WHgtEddLog <-
  glm.nb(
    nSeqs ~ Height_cm + Log.in.View + Squirrel_EDD_WSp +
      offset(log(Deploy.Duration)),
    data = foxsqDataWin
  )
ModListWinFSq[[25]] <- glm.nb.WHgtOakLog <-
  glm.nb(
    nSeqs ~ Height_cm + Log.in.View + OakDBH +
      offset(log(Deploy.Duration)),
    data = foxsqDataWin
  )
modnamesW <- c("Full","Intercept","Log","Height","Oak","Stems","EDD","Log + Stems","Log + EDD","Log + Oak","Log + Height","EDD + Height","EDD + Stems","EDD + Oak","Stems + Height","Stems + Oak","Height + Oak","Log + Stems + Height","Log + Stems + Oak","Log + Stems + EDD","Height + Stems + EDD","Height + Stems + Oak","Height + EDD + Oak","Height + EDD + Log","Height + Oak + Log")
modtabWinFSq <- aictab(cand.set = ModListWinFSq, modnames = modnamesW)
modtabWinFSq

summary(glm.nb.Wlog) #best model is intercept only model

#Summed Model Weights for Squirrels in Winter from full table. 
WinSW_edd <- sum(modtabWinFSq$AICcWt[grep("EDD", modtabWinFSq$Modnames)])
WinSW_OakDBH <- sum(modtabWinFSq$AICcWt[grep("Oak", modtabWinFSq$Modnames)])
WinSW_Stems <- sum(modtabWinFSq$AICcWt[grep("Stems", modtabWinFSq$Modnames)])
WinSW_Log <- sum(modtabWinFSq$AICcWt[grep("Log", modtabWinFSq$Modnames)])
WinSW_Height <- sum(modtabWinFSq$AICcWt[grep("Height", modtabWinFSq$Modnames)])

#Explained deviance of best model
1 - glm.nb.WI$deviance / glm.nb.WI$null.deviance #0.0
#Explained deviance of full model
1 - glm.nb.fullW$deviance / glm.nb.fullW$null.deviance #0.05966105


# Spring Fox Squirrel NB Full Model Selection ---------------------------------

#First test to see if interaction between log and EDD improves that two variable model

glm.nb.SPlog_Edd <-
  glm.nb(nSeqs ~ Log.in.View + Squirrel_EDD_WSp + offset(log(Deploy.Duration)), data = foxsqDataSpr)
AICc(glm.nb.SPlog_Edd)

glm.nb.SPlogXEdd <-
  glm.nb(nSeqs ~ Log.in.View * Squirrel_EDD_WSp + offset(log(Deploy.Duration)), data = foxsqDataSpr)
AICc(glm.nb.SPlogXEdd) #There is support for interaction, so this has to be added to each model where these two appear together.


#Compare height alone to height quadratic
glm.nb.SPHgt2 <-
  glm.nb(
    nSeqs ~ poly(Height_cm, 2, raw = T) +
      offset(log(Deploy.Duration)),
    data = foxsqDataSpr
  )

glm.nb.SPhgt <-
  glm.nb(nSeqs ~ Height_cm + offset(log(Deploy.Duration)), data = foxsqDataSpr)

AICc(glm.nb.SPHgt2)
AICc(glm.nb.SPhgt) #there is support for a quadratic effect, so it needs to be added to all models

#Full Spring Fox Squirrel Model Selection

ModListSprFSq <- list(NA)
ModListSprFSq[[1]] <- glm.nb.fullSP <-
  glm.nb(
    nSeqs ~ Log.in.View + Squirrel_EDD_WSp + poly(Height_cm, 2, raw = T) + log10(Num_Stems) + OakDBH +
      offset(log(Deploy.Duration)),
    data = foxsqDataSpr
  )

ModListSprFSq[[2]] <- glm.nb.SPI <- glm.nb(nSeqs ~ 1 + offset(log(Deploy.Duration)), data = foxsqDataSpr)

ModListSprFSq[[3]] <- glm.nb.SPlog <-
  glm.nb(nSeqs ~ Log.in.View + offset(log(Deploy.Duration)), data = foxsqDataSpr)

ModListSprFSq[[4]] <- glm.nb.SPhgt <-
  glm.nb(nSeqs ~ poly(Height_cm, 2, raw = T) + offset(log(Deploy.Duration)), data = foxsqDataSpr)

ModListSprFSq[[5]] <- glm.nb.SPoak <-
  glm.nb(nSeqs ~ OakDBH + offset(log(Deploy.Duration)), data = foxsqDataSpr)

ModListSprFSq[[6]] <- glm.nb.SPstems <-
  glm.nb(nSeqs ~ log10(Num_Stems) + offset(log(Deploy.Duration)), data = foxsqDataSpr)

ModListSprFSq[[7]] <- glm.nb.SPedd <-
  glm.nb(nSeqs ~ Squirrel_EDD_WSp + offset(log(Deploy.Duration)), data = foxsqDataSpr)

ModListSprFSq[[8]] <- glm.nb.SPlog_stems <-
  glm.nb(nSeqs ~ Log.in.View + log10(Num_Stems) + offset(log(Deploy.Duration)), data = foxsqDataSpr)

ModListSprFSq[[9]] <-glm.nb.SPlog_Edd <-
  glm.nb(nSeqs ~ Log.in.View + Squirrel_EDD_WSp + offset(log(Deploy.Duration)), data = foxsqDataSpr)

ModListSprFSq[[10]] <- glm.nb.SPlog_Oak <-
  glm.nb(nSeqs ~ Log.in.View + OakDBH + offset(log(Deploy.Duration)), data = foxsqDataSpr)

ModListSprFSq[[11]] <- glm.nb.SPlog_Hgt <-
  glm.nb(nSeqs ~ Log.in.View + poly(Height_cm, 2, raw = T) + offset(log(Deploy.Duration)), data = foxsqDataSpr)

ModListSprFSq[[12]] <- glm.nb.SPEdd_Hgt <-
  glm.nb(nSeqs ~ Squirrel_EDD_WSp + poly(Height_cm, 2, raw = T) + offset(log(Deploy.Duration)), data = foxsqDataSpr)

ModListSprFSq[[13]] <-  glm.nb.SPEdd_stems <-
  glm.nb(nSeqs ~ Squirrel_EDD_WSp + log10(Num_Stems) + offset(log(Deploy.Duration)), data = foxsqDataSpr)

ModListSprFSq[[14]] <- glm.nb.SPEdd_oak <-
  glm.nb(nSeqs ~ Squirrel_EDD_WSp + OakDBH + offset(log(Deploy.Duration)), data = foxsqDataSpr)

ModListSprFSq[[15]] <- glm.nb.SPstems_hgt <-
  glm.nb(nSeqs ~ poly(Height_cm, 2, raw = T) + log10(Num_Stems) + offset(log(Deploy.Duration)), data = foxsqDataSpr)

ModListSprFSq[[16]] <- glm.nb.SPstems_oak <-
  glm.nb(nSeqs ~ OakDBH + log10(Num_Stems) + offset(log(Deploy.Duration)), data = foxsqDataSpr)

ModListSprFSq[[17]] <- glm.nb.SPhgt_oak <-
  glm.nb(nSeqs ~ poly(Height_cm, 2, raw = T) + OakDBH + offset(log(Deploy.Duration)), data = foxsqDataSpr)

ModListSprFSq[[18]] <- glm.nb.SPlogStemsHgt <-
  glm.nb(
    nSeqs ~ Log.in.View + poly(Height_cm, 2, raw = T) + log10(Num_Stems) +
      offset(log(Deploy.Duration)),
    data = foxsqDataSpr
  )

ModListSprFSq[[19]] <- glm.nb.SPlogStemsOak <-
  glm.nb(
    nSeqs ~ Log.in.View + OakDBH + log10(Num_Stems) +
      offset(log(Deploy.Duration)),
    data = foxsqDataSpr
  )

ModListSprFSq[[20]] <- glm.nb.SPlogStemsEdd <-
  glm.nb(
    nSeqs ~ Log.in.View + Squirrel_EDD_WSp + log10(Num_Stems) +
      offset(log(Deploy.Duration)),
    data = foxsqDataSpr
  )

ModListSprFSq[[21]] <- glm.nb.SPHgtStemsEdd <-
  glm.nb(
    nSeqs ~ poly(Height_cm, 2, raw = T) + Squirrel_EDD_WSp + log10(Num_Stems) +
      offset(log(Deploy.Duration)),
    data = foxsqDataSpr
  )

ModListSprFSq[[22]] <- glm.nb.SPHgtStemsOak <-
  glm.nb(
    nSeqs ~ poly(Height_cm, 2, raw = T) + OakDBH + log10(Num_Stems) +
      offset(log(Deploy.Duration)),
    data = foxsqDataSpr
  )

ModListSprFSq[[23]] <- glm.nb.SPHgtEddOak <-
  glm.nb(
    nSeqs ~ poly(Height_cm, 2, raw = T) + OakDBH + Squirrel_EDD_WSp +
      offset(log(Deploy.Duration)),
    data = foxsqDataSpr
  )
ModListSprFSq[[24]] <- glm.nb.SPHgtEddLog <-
  glm.nb(
    nSeqs ~ poly(Height_cm, 2, raw = T) + Log.in.View + Squirrel_EDD_WSp +
      offset(log(Deploy.Duration)),
    data = foxsqDataSpr
  )
ModListSprFSq[[25]] <- glm.nb.SPHgtOakLog <-
  glm.nb(
    nSeqs ~ poly(Height_cm, 2, raw = T) + Log.in.View + OakDBH +
      offset(log(Deploy.Duration)),
    data = foxsqDataSpr
  )
modnamesSP <- c("Full","Intercept","Log","Height2","Oak","Stems","EDD","Log + Stems","Log + EDD","Log + Oak","Log + Height2","EDD + Height2","EDD + Stems","EDD + Oak","Stems + Height2","Stems + Oak","Height2 + Oak","Log + Stems + Height2","Log + Stems + Oak","Log + Stems + EDD","Height2 + Stems + EDD","Height2 + Stems + Oak","Height2 + EDD + Oak","Height2 + EDD + Log","Height2 + Oak + Log")
modtabSprFSq <- aictab(cand.set = ModListSprFSq, modnames = modnamesSP)
modtabSprFSq

summary(glm.nb.SPlog) #best model has intercept only

#Summed Model Weights for Fox Squirrels in Spring from full table. 
SprSW_edd <- sum(modtabSprFSq$AICcWt[grep("EDD", modtabSprFSq$Modnames)])
SprSW_OakDBH <- sum(modtabSprFSq$AICcWt[grep("Oak", modtabSprFSq$Modnames)])
SprSW_Stems <- sum(modtabSprFSq$AICcWt[grep("Stems", modtabSprFSq$Modnames)])
SprSW_Log <- sum(modtabSprFSq$AICcWt[grep("Log", modtabSprFSq$Modnames)])
SprSW_Height <- sum(modtabSprFSq$AICcWt[grep("Height", modtabSprFSq$Modnames)])

#Explained deviance of best model
1 - glm.nb.SPI$deviance / glm.nb.SPI$null.deviance #0.00
#Explained deviance of full model
1 - glm.nb.fullSP$deviance / glm.nb.fullSP$null.deviance #0.2632075
