###Regression Analysis of SCBI Camera Grid Data for Raccoon

require(AICcmodavg)
require(MASS) #for glm.nb function
source("scripts/pairsPannelFunctions.r")


# Data Import and Exploration - Raccoon ---------------------------------------
load("data/racDataSum.RData")
load("data/racDataFall.RData")
load("data/racDataWin.RData")
load("data/racDataSpr.RData")

pairs(racDataSum[,c(4,7,13,14,15,17)], diag.panel = panel.hist, lower.panel = panel.smooth, upper.panel = panel.cor)
str(racDataSum)

#Looking at how many cameras picked up rac in each season. Looks like enough in each, though a bit low in summer

par.default <- par(no.readonly = T)
par (mfrow = c(2,2))
barplot(racDataSum$nSeqs,
        ylab = "# of sequences",
        xlab = "Individual Cameras",
        main = "Summer")
barplot(racDataFall$nSeqs,
        ylab = "# of sequences",
        xlab = "Individual Cameras",
        main = "Fall")
barplot(racDataWin$nSeqs,
        ylab = "# of sequences",
        xlab = "Individual Cameras",
        main = "Winter")
barplot(racDataSpr$nSeqs,
        ylab = "# of sequences",
        xlab = "Individual Cameras",
        main = "Spring")
par(par.default)

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


# Negative Binomial Regression - Raccoon Summer--------------------------------------
#Summer

ModListSumRac <- list(NA)

ModListSumRac[[1]] <- glm.nb.fullS <-
  glm.nb(
    nSeqs ~ Log.in.View + Raccoon.EDD + Height_cm + log10(Num_Stems) + OakDBH +
      offset(log(Deploy.Duration)),
    data = racDataSum
  )

ModListSumRac[[2]] <- glm.nb.SI <- glm.nb(nSeqs ~ 1 + offset(log(Deploy.Duration)), data = racDataSum)

ModListSumRac[[3]] <- glm.nb.Slog <-
  glm.nb(nSeqs ~ Log.in.View + offset(log(Deploy.Duration)), data = racDataSum)

ModListSumRac[[4]] <- glm.nb.Shgt <-
  glm.nb(nSeqs ~ Height_cm + offset(log(Deploy.Duration)), data = racDataSum)

ModListSumRac[[5]] <- glm.nb.Soak <-
  glm.nb(nSeqs ~ OakDBH + offset(log(Deploy.Duration)), data = racDataSum)

ModListSumRac[[6]] <- glm.nb.Sstems <-
  glm.nb(nSeqs ~ log10(Num_Stems) + offset(log(Deploy.Duration)), data = racDataSum)

ModListSumRac[[7]] <- glm.nb.Sedd <-
  glm.nb(nSeqs ~ Raccoon.EDD + offset(log(Deploy.Duration)), data = racDataSum)

ModListSumRac[[8]] <- glm.nb.Slog_stems <-
  glm.nb(nSeqs ~ Log.in.View + log10(Num_Stems) + offset(log(Deploy.Duration)), data = racDataSum)

ModListSumRac[[9]] <-glm.nb.Slog_Edd <-
  glm.nb(nSeqs ~ Log.in.View + Raccoon.EDD + offset(log(Deploy.Duration)), data = racDataSum)

ModListSumRac[[10]] <- glm.nb.Slog_Oak <-
  glm.nb(nSeqs ~ Log.in.View + OakDBH + offset(log(Deploy.Duration)), data = racDataSum)

ModListSumRac[[11]] <- glm.nb.Slog_Hgt <-
  glm.nb(nSeqs ~ Log.in.View + Height_cm + offset(log(Deploy.Duration)), data = racDataSum)

ModListSumRac[[12]] <- glm.nb.SEdd_Hgt <-
  glm.nb(nSeqs ~ Raccoon.EDD + Height_cm + offset(log(Deploy.Duration)), data = racDataSum)

ModListSumRac[[13]] <-  glm.nb.SEdd_stems <-
  glm.nb(nSeqs ~ Raccoon.EDD + log10(Num_Stems) + offset(log(Deploy.Duration)), data = racDataSum)

ModListSumRac[[14]] <- glm.nb.SEdd_oak <-
  glm.nb(nSeqs ~ Raccoon.EDD + OakDBH + offset(log(Deploy.Duration)), data = racDataSum)

ModListSumRac[[15]] <- glm.nb.Sstems_hgt <-
  glm.nb(nSeqs ~ Height_cm + log10(Num_Stems) + offset(log(Deploy.Duration)), data = racDataSum)

ModListSumRac[[16]] <- glm.nb.Sstems_oak <-
  glm.nb(nSeqs ~ OakDBH + log10(Num_Stems) + offset(log(Deploy.Duration)), data = racDataSum)

ModListSumRac[[17]] <- glm.nb.Shgt_oak <-
  glm.nb(nSeqs ~ Height_cm + OakDBH + offset(log(Deploy.Duration)), data = racDataSum)

ModListSumRac[[18]] <- glm.nb.SlogStemsHgt <-
  glm.nb(
    nSeqs ~ Log.in.View + Height_cm + log10(Num_Stems) +
      offset(log(Deploy.Duration)),
    data = racDataSum
  )

ModListSumRac[[19]] <- glm.nb.SlogStemsOak <-
  glm.nb(
    nSeqs ~ Log.in.View + OakDBH + log10(Num_Stems) +
      offset(log(Deploy.Duration)),
    data = racDataSum
  )

ModListSumRac[[20]] <- glm.nb.SlogStemsEdd <-
  glm.nb(
    nSeqs ~ Log.in.View + Raccoon.EDD + log10(Num_Stems) +
      offset(log(Deploy.Duration)),
    data = racDataSum
  )

ModListSumRac[[21]] <- glm.nb.SHgtStemsEdd <-
  glm.nb(
    nSeqs ~ Height_cm + Raccoon.EDD + log10(Num_Stems) +
      offset(log(Deploy.Duration)),
    data = racDataSum
  )

ModListSumRac[[22]] <- glm.nb.SHgtStemsOak <-
  glm.nb(
    nSeqs ~ Height_cm + OakDBH + log10(Num_Stems) +
      offset(log(Deploy.Duration)),
    data = racDataSum
  )

ModListSumRac[[23]] <- glm.nb.SHgtEddOak <-
  glm.nb(
    nSeqs ~ Height_cm + OakDBH + Raccoon.EDD +
      offset(log(Deploy.Duration)),
    data = racDataSum
  )
ModListSumRac[[24]] <- glm.nb.SHgtEddLog <-
  glm.nb(
    nSeqs ~ Height_cm + Log.in.View + Raccoon.EDD +
      offset(log(Deploy.Duration)),
    data = racDataSum
  )
ModListSumRac[[25]] <- glm.nb.SHgtOakLog <-
  glm.nb(
    nSeqs ~ Height_cm + Log.in.View + OakDBH +
      offset(log(Deploy.Duration)),
    data = racDataSum
  )
modnamesS <- c("Full","Intercept","Log","Height","Oak","Stems","EDD","Log + Stems","Log + EDD","Log + Oak","Log + Height","EDD + Height","EDD + Stems","EDD + Oak","Stems + Height","Stems + Oak","Height + Oak","Log + Stems + Height","Log + Stems + Oak","Log + Stems + EDD","Height + Stems + EDD","Height + Stems + Oak","Height + EDD + Oak","Height + EDD + Log","Height + Oak + Log")
modtabSumRac <- aictab(cand.set = ModListSumRac, modnames = modnamesS)
modtabSumRac

#Explained Deviance of best model
1 - glm.nb.Sstems$deviance / glm.nb.Sstems$null.deviance #0.1365
#Explained Deviance of best 3 variable model
1 - glm.nb.SlogStemsEdd$deviance / glm.nb.SlogStemsEdd$null.deviance #0.213
#Explained Deviance of Full model
1 - glm.nb.fullS$deviance / glm.nb.fullS$null.deviance #0.215

#We would conclude only stems is important. For summer, raccoon captures are more common as the # of stems in front of the camera increases. This is hard to explain, but this is the strongest effect in summer. Looking at the plot though, this is largely determined by the single site with high stem density.

#Response Curves for Num_Stems. No other parameters were important.
summary(racDataSum$Num_Stems)
seq.stems <- 10^seq(min(log10(racDataSum$Num_Stems)), max(log10(racDataSum$Num_Stems)), length.out = 100)
nd.seq.stems <- data.frame(Deploy.Duration = 61, Num_Stems = seq.stems)
pred.stems <- predict(glm.nb.Sstems, newdata = nd.seq.stems, se=TRUE, type = "link")
plot(nSeqs ~ Num_Stems, data=racDataSum, pch=16, col=rgb(.75,.25,0,.5), las=1, xlab = "Total stems", ylab = "# Raccoon Sequences in Summer")
lines(nd.seq.stems$Num_Stems, exp(pred.stems$fit), col="dark olive green")
lines(nd.seq.stems$Num_Stems, exp(pred.stems$fit + 1.96*pred.stems$se.fit), lty=2, col="dark olive green")
lines(nd.seq.stems$Num_Stems, exp(pred.stems$fit - 1.96*pred.stems$se.fit), lty=2, col="dark olive green")


# Raccoon in Fall NB Full Model Selection ---------------------------------


ModListFallRac <- list(NA)

ModListFallRac[[1]] <- glm.nb.fullF <-
  glm.nb(
    nSeqs ~ Log.in.View + Raccoon.EDD + Height_cm + log10(Num_Stems) + OakDBH +
      offset(log(Deploy.Duration)),
    data = racDataFall
  )

ModListFallRac[[2]] <- glm.nb.FI <- glm.nb(nSeqs ~ 1 + offset(log(Deploy.Duration)), data = racDataFall)

ModListFallRac[[3]] <- glm.nb.Flog <-
  glm.nb(nSeqs ~ Log.in.View + offset(log(Deploy.Duration)), data = racDataFall)

ModListFallRac[[4]] <- glm.nb.Fhgt <-
  glm.nb(nSeqs ~ Height_cm + offset(log(Deploy.Duration)), data = racDataFall)

ModListFallRac[[5]] <- glm.nb.Foak <-
  glm.nb(nSeqs ~ OakDBH + offset(log(Deploy.Duration)), data = racDataFall)

ModListFallRac[[6]] <- glm.nb.Fstems <-
  glm.nb(nSeqs ~ log10(Num_Stems) + offset(log(Deploy.Duration)), data = racDataFall)

ModListFallRac[[7]] <- glm.nb.Fedd <-
  glm.nb(nSeqs ~ Raccoon.EDD + offset(log(Deploy.Duration)), data = racDataFall)

ModListFallRac[[8]] <- glm.nb.Flog_stems <-
  glm.nb(nSeqs ~ Log.in.View + log10(Num_Stems) + offset(log(Deploy.Duration)), data = racDataFall)

ModListFallRac[[9]] <-glm.nb.Flog_Edd <-
  glm.nb(nSeqs ~ Log.in.View + Raccoon.EDD + offset(log(Deploy.Duration)), data = racDataFall)

ModListFallRac[[10]] <- glm.nb.Flog_Oak <-
  glm.nb(nSeqs ~ Log.in.View + OakDBH + offset(log(Deploy.Duration)), data = racDataFall)

ModListFallRac[[11]] <- glm.nb.Flog_Hgt <-
  glm.nb(nSeqs ~ Log.in.View + Height_cm + offset(log(Deploy.Duration)), data = racDataFall)

ModListFallRac[[12]] <- glm.nb.FEdd_Hgt <-
  glm.nb(nSeqs ~ Raccoon.EDD + Height_cm + offset(log(Deploy.Duration)), data = racDataFall)

ModListFallRac[[13]] <-  glm.nb.FEdd_stems <-
  glm.nb(nSeqs ~ Raccoon.EDD + log10(Num_Stems) + offset(log(Deploy.Duration)), data = racDataFall)

ModListFallRac[[14]] <- glm.nb.FEdd_oak <-
  glm.nb(nSeqs ~ Raccoon.EDD + OakDBH + offset(log(Deploy.Duration)), data = racDataFall)

ModListFallRac[[15]] <- glm.nb.Fstems_hgt <-
  glm.nb(nSeqs ~ Height_cm + log10(Num_Stems) + offset(log(Deploy.Duration)), data = racDataFall)

ModListFallRac[[16]] <- glm.nb.Fstems_oak <-
  glm.nb(nSeqs ~ OakDBH + log10(Num_Stems) + offset(log(Deploy.Duration)), data = racDataFall)

ModListFallRac[[17]] <- glm.nb.Fhgt_oak <-
  glm.nb(nSeqs ~ Height_cm + OakDBH + offset(log(Deploy.Duration)), data = racDataFall)

ModListFallRac[[18]] <- glm.nb.FlogStemsHgt <-
  glm.nb(
    nSeqs ~ Log.in.View + Height_cm + log10(Num_Stems) +
      offset(log(Deploy.Duration)),
    data = racDataFall
  )

ModListFallRac[[19]] <- glm.nb.FlogStemsOak <-
  glm.nb(
    nSeqs ~ Log.in.View + OakDBH + log10(Num_Stems) +
      offset(log(Deploy.Duration)),
    data = racDataFall
  )

ModListFallRac[[20]] <- glm.nb.FlogStemsEdd <-
  glm.nb(
    nSeqs ~ Log.in.View + Raccoon.EDD + log10(Num_Stems) +
      offset(log(Deploy.Duration)),
    data = racDataFall
  )

ModListFallRac[[21]] <- glm.nb.FHgtStemsEdd <-
  glm.nb(
    nSeqs ~ Height_cm + Raccoon.EDD + log10(Num_Stems) +
      offset(log(Deploy.Duration)),
    data = racDataFall
  )

ModListFallRac[[22]] <- glm.nb.FHgtStemsOak <-
  glm.nb(
    nSeqs ~ Height_cm + OakDBH + log10(Num_Stems) +
      offset(log(Deploy.Duration)),
    data = racDataFall
  )

ModListFallRac[[23]] <- glm.nb.FHgtEddOak <-
  glm.nb(
    nSeqs ~ Height_cm + OakDBH + Raccoon.EDD +
      offset(log(Deploy.Duration)),
    data = racDataFall
  )
ModListFallRac[[24]] <- glm.nb.FHgtEddLog <-
  glm.nb(
    nSeqs ~ Height_cm + Log.in.View + Raccoon.EDD +
      offset(log(Deploy.Duration)),
    data = racDataFall
  )
ModListFallRac[[25]] <- glm.nb.FHgtOakLog <-
  glm.nb(
    nSeqs ~ Height_cm + Log.in.View + OakDBH +
      offset(log(Deploy.Duration)),
    data = racDataFall
  )
modnamesF <- c("Full","Intercept","Log","Height","Oak","Stems","EDD","Log + Stems","Log + EDD","Log + Oak","Log + Height","EDD + Height","EDD + Stems","EDD + Oak","Stems + Height","Stems + Oak","Height + Oak","Log + Stems + Height","Log + Stems + Oak","Log + Stems + EDD","Height + Stems + EDD","Height + Stems + Oak","Height + EDD + Oak","Height + EDD + Log","Height + Oak + Log")
modtabFallRac <- aictab(cand.set = ModListFallRac, modnames = modnamesF)
modtabFallRac #The log model is best, not much support for EDD, but the Stems + Oak model also has reasonable support. Should discuss all three.
summary(glm.nb.Fstems_oak)


#Explained Deviance of best model A
1 - glm.nb.Flog$deviance / glm.nb.Flog$null.deviance #0.140
#Explained Deviance of best model B
1 - glm.nb.Fstems_oak$deviance / glm.nb.Fstems_oak$null.deviance #0.219
#Explained Deviance of best 3 variable model
1 - glm.nb.FlogStemsOak$deviance / glm.nb.FlogStemsOak$null.deviance #0.267
#Explained Deviance of Full model
1 - glm.nb.fullF$deviance / glm.nb.fullF$null.deviance #0.3296


#Plots for Raccoon in Fall results
plot(nSeqs ~ Log.in.View, data = racDataFall, main = "Fall")

#Oak effects in Fall
summary(racDataSum$OakDBH)
seq.oak <- seq(min(racDataFall$OakDBH), max(racDataFall$OakDBH), length.out = 100)
nd.seq.oak <- data.frame(Deploy.Duration = 61, OakDBH = seq.oak, Num_Stems = median(racDataFall$Num_Stems))
pred.oak <- predict(glm.nb.Fstems_oak, newdata = nd.seq.oak, se=TRUE, type = "link")
plot(nSeqs ~ OakDBH, data=racDataFall, pch=16, col=rgb(.75,.25,0,.5), las=1, xlab = "Total Oak DBH", ylab = "# Raccoon Sequences in Fall")
lines(nd.seq.oak$OakDBH, exp(pred.oak$fit), col="dark olive green")
lines(nd.seq.oak$OakDBH, exp(pred.oak$fit + 1.96*pred.oak$se.fit), lty=2, col="dark olive green")
lines(nd.seq.oak$OakDBH, exp(pred.oak$fit - 1.96*pred.oak$se.fit), lty=2, col="dark olive green")

#Response Curve for Num_Stems in Fall.
summary(racDataSum$Num_Stems)
seq.stems <- 10^seq(min(log10(racDataSum$Num_Stems)), max(log10(racDataSum$Num_Stems)), length.out = 100)
nd.seq.stemsF <- data.frame(Deploy.Duration = 61, Num_Stems = seq.stems, OakDBH = median(racDataFall$OakDBH))
pred.stemsF <- predict(glm.nb.Fstems_oak, newdata = nd.seq.stemsF, se=TRUE, type = "link")
plot(nSeqs ~ Num_Stems, data=racDataFall, pch=16, col=rgb(.75,.25,0,.5), las=1, xlab = "Total stems", ylab = "# Raccoon Sequences in Fall")
lines(nd.seq.stemsF$Num_Stems, exp(pred.stemsF$fit), col="dark olive green")
lines(nd.seq.stemsF$Num_Stems, exp(pred.stemsF$fit + 1.96*pred.stemsF$se.fit), lty=2, col="dark olive green")
lines(nd.seq.stemsF$Num_Stems, exp(pred.stemsF$fit - 1.96*pred.stemsF$se.fit), lty=2, col="dark olive green")


# Raccoons in Winter NB Full Model Selection ------------------------------

ModListWinRac <- list(NA)
ModListWinRac[[1]] <- glm.nb.fullW <-
  glm.nb(
    nSeqs ~ Log.in.View + Raccoon.EDD + Height_cm + log10(Num_Stems) + OakDBH +
      offset(log(Deploy.Duration)),
    data = racDataWin
  )

ModListWinRac[[2]] <- glm.nb.WI <- glm.nb(nSeqs ~ 1 + offset(log(Deploy.Duration)), data = racDataWin)

ModListWinRac[[3]] <- glm.nb.Wlog <-
  glm.nb(nSeqs ~ Log.in.View + offset(log(Deploy.Duration)), data = racDataWin)

ModListWinRac[[4]] <- glm.nb.Whgt <-
  glm.nb(nSeqs ~ Height_cm + offset(log(Deploy.Duration)), data = racDataWin)

ModListWinRac[[5]] <- glm.nb.Woak <-
  glm.nb(nSeqs ~ OakDBH + offset(log(Deploy.Duration)), data = racDataWin)

ModListWinRac[[6]] <- glm.nb.Wstems <-
  glm.nb(nSeqs ~ log10(Num_Stems) + offset(log(Deploy.Duration)), data = racDataWin)

ModListWinRac[[7]] <- glm.nb.Wedd <-
  glm.nb(nSeqs ~ Raccoon.EDD + offset(log(Deploy.Duration)), data = racDataWin)

ModListWinRac[[8]] <- glm.nb.Wlog_stems <-
  glm.nb(nSeqs ~ Log.in.View + log10(Num_Stems) + offset(log(Deploy.Duration)), data = racDataWin)

ModListWinRac[[9]] <-glm.nb.Wlog_Edd <-
  glm.nb(nSeqs ~ Log.in.View + Raccoon.EDD + offset(log(Deploy.Duration)), data = racDataWin)

ModListWinRac[[10]] <- glm.nb.Wlog_Oak <-
  glm.nb(nSeqs ~ Log.in.View + OakDBH + offset(log(Deploy.Duration)), data = racDataWin)

ModListWinRac[[11]] <- glm.nb.Wlog_Hgt <-
  glm.nb(nSeqs ~ Log.in.View + Height_cm + offset(log(Deploy.Duration)), data = racDataWin)

ModListWinRac[[12]] <- glm.nb.WEdd_Hgt <-
  glm.nb(nSeqs ~ Raccoon.EDD + Height_cm + offset(log(Deploy.Duration)), data = racDataWin)

ModListWinRac[[13]] <-  glm.nb.WEdd_stems <-
  glm.nb(nSeqs ~ Raccoon.EDD + log10(Num_Stems) + offset(log(Deploy.Duration)), data = racDataWin)

ModListWinRac[[14]] <- glm.nb.WEdd_oak <-
  glm.nb(nSeqs ~ Raccoon.EDD + OakDBH + offset(log(Deploy.Duration)), data = racDataWin)

ModListWinRac[[15]] <- glm.nb.Wstems_hgt <-
  glm.nb(nSeqs ~ Height_cm + log10(Num_Stems) + offset(log(Deploy.Duration)), data = racDataWin)

ModListWinRac[[16]] <- glm.nb.Wstems_oak <-
  glm.nb(nSeqs ~ OakDBH + log10(Num_Stems) + offset(log(Deploy.Duration)), data = racDataWin)

ModListWinRac[[17]] <- glm.nb.Whgt_oak <-
  glm.nb(nSeqs ~ Height_cm + OakDBH + offset(log(Deploy.Duration)), data = racDataWin)

ModListWinRac[[18]] <- glm.nb.WlogStemsHgt <-
  glm.nb(
    nSeqs ~ Log.in.View + Height_cm + log10(Num_Stems) +
      offset(log(Deploy.Duration)),
    data = racDataWin
  )

ModListWinRac[[19]] <- glm.nb.WlogStemsOak <-
  glm.nb(
    nSeqs ~ Log.in.View + OakDBH + log10(Num_Stems) +
      offset(log(Deploy.Duration)),
    data = racDataWin
  )

ModListWinRac[[20]] <- glm.nb.WlogStemsEdd <-
  glm.nb(
    nSeqs ~ Log.in.View + Raccoon.EDD + log10(Num_Stems) +
      offset(log(Deploy.Duration)),
    data = racDataWin
  )

ModListWinRac[[21]] <- glm.nb.WHgtStemsEdd <-
  glm.nb(
    nSeqs ~ Height_cm + Raccoon.EDD + log10(Num_Stems) +
      offset(log(Deploy.Duration)),
    data = racDataWin
  )

ModListWinRac[[22]] <- glm.nb.WHgtStemsOak <-
  glm.nb(
    nSeqs ~ Height_cm + OakDBH + log10(Num_Stems) +
      offset(log(Deploy.Duration)),
    data = racDataWin
  )

ModListWinRac[[23]] <- glm.nb.WHgtEddOak <-
  glm.nb(
    nSeqs ~ Height_cm + OakDBH + Raccoon.EDD +
      offset(log(Deploy.Duration)),
    data = racDataWin
  )
ModListWinRac[[24]] <- glm.nb.WHgtEddLog <-
  glm.nb(
    nSeqs ~ Height_cm + Log.in.View + Raccoon.EDD +
      offset(log(Deploy.Duration)),
    data = racDataWin
  )
ModListWinRac[[25]] <- glm.nb.WHgtOakLog <-
  glm.nb(
    nSeqs ~ Height_cm + Log.in.View + OakDBH +
      offset(log(Deploy.Duration)),
    data = racDataWin
  )
modnamesW <- c("Full","Intercept","Log","Height","Oak","Stems","EDD","Log + Stems","Log + EDD","Log + Oak","Log + Height","EDD + Height","EDD + Stems","EDD + Oak","Stems + Height","Stems + Oak","Height + Oak","Log + Stems + Height","Log + Stems + Oak","Log + Stems + EDD","Height + Stems + EDD","Height + Stems + Oak","Height + EDD + Oak","Height + EDD + Log","Height + Oak + Log")
modtabWinRac <- aictab(cand.set = ModListWinRac, modnames = modnamesW)
modtabWinRac

summary(glm.nb.Wlog_Hgt)
summary(glm.nb.Wstems_hgt)

#Explained Deviance of Best Model
1 - glm.nb.Wlog_Hgt$deviance / glm.nb.Wlog_Hgt$null.deviance #0.284
#Explained Deviance of best 3 variable model
1 - glm.nb.WlogStemsHgt$deviance / glm.nb.WlogStemsHgt$null.deviance #0.323796
#Explained Deviance of Full Model
1 - glm.nb.fullW$deviance / glm.nb.fullW$null.deviance #0.361776

#Plotting Height Effects for Raccoon in Winter
summary(glm.nb.Wlog_stems)
summary(glm.nb.Whgt_stems) #In Winter, raccoon captures are more frequent at cameras that are set lower and which have more tree stems in front of them. Height seems the most important.

#In Winter, raccoon captures are more frequent at cameras that are set lower and which have more tree stems in front of them

summary(glm.nb.Whgt_stems) 

seq.hgt <- seq(min(racDataWin$Height_cm), max(racDataWin$Height_cm), length.out = 100)
med.stems <- median(racDataWin$Num_Stems)

nd.seq.hgt <- data.frame(Deploy.Duration = 61, Height_cm = seq.hgt, Num_Stems = med.stems)

pred.hgt <- predict(glm.nb.Whgt_stems, newdata = nd.seq.hgt, se=TRUE, type = "link")

plot(nSeqs ~ Height_cm, data=racDataWin, pch=16, col=rgb(.75,.25,0,.5), las=1, xlab = "Camera Height (cm)", ylab = "# Raccoon Sequences in Winter")
lines(nd.seq.hgt$Height_cm, exp(pred.hgt$fit), col="dark olive green")
lines(nd.seq.hgt$Height_cm, exp(pred.hgt$fit + 1.96*pred.hgt$se.fit), lty=2, col="dark olive green")
lines(nd.seq.hgt$Height_cm, exp(pred.hgt$fit - 1.96*pred.hgt$se.fit), lty=2, col="dark olive green")




# Raccoon in Spring NB Full Model Selection -------------------------------

ModListSprRac <- list(NA)
ModListSprRac[[1]] <- glm.nb.fullSP <-
  glm.nb(
    nSeqs ~ Log.in.View + Raccoon.EDD + Height_cm + log10(Num_Stems) + OakDBH +
      offset(log(Deploy.Duration)),
    data = racDataSpr
  )

ModListSprRac[[2]] <- glm.nb.SPI <- glm.nb(nSeqs ~ 1 + offset(log(Deploy.Duration)), data = racDataSpr)

ModListSprRac[[3]] <- glm.nb.SPlog <-
  glm.nb(nSeqs ~ Log.in.View + offset(log(Deploy.Duration)), data = racDataSpr)

ModListSprRac[[4]] <- glm.nb.SPhgt <-
  glm.nb(nSeqs ~ Height_cm + offset(log(Deploy.Duration)), data = racDataSpr)

ModListSprRac[[5]] <- glm.nb.SPoak <-
  glm.nb(nSeqs ~ OakDBH + offset(log(Deploy.Duration)), data = racDataSpr)

ModListSprRac[[6]] <- glm.nb.SPstems <-
  glm.nb(nSeqs ~ log10(Num_Stems) + offset(log(Deploy.Duration)), data = racDataSpr)

ModListSprRac[[7]] <- glm.nb.SPedd <-
  glm.nb(nSeqs ~ Raccoon.EDD + offset(log(Deploy.Duration)), data = racDataSpr)

ModListSprRac[[8]] <- glm.nb.SPlog_stems <-
  glm.nb(nSeqs ~ Log.in.View + log10(Num_Stems) + offset(log(Deploy.Duration)), data = racDataSpr)

ModListSprRac[[9]] <-glm.nb.SPlog_Edd <-
  glm.nb(nSeqs ~ Log.in.View + Raccoon.EDD + offset(log(Deploy.Duration)), data = racDataSpr)

ModListSprRac[[10]] <- glm.nb.SPlog_Oak <-
  glm.nb(nSeqs ~ Log.in.View + OakDBH + offset(log(Deploy.Duration)), data = racDataSpr)

ModListSprRac[[11]] <- glm.nb.SPlog_Hgt <-
  glm.nb(nSeqs ~ Log.in.View + Height_cm + offset(log(Deploy.Duration)), data = racDataSpr)

ModListSprRac[[12]] <- glm.nb.SPEdd_Hgt <-
  glm.nb(nSeqs ~ Raccoon.EDD + Height_cm + offset(log(Deploy.Duration)), data = racDataSpr)

ModListSprRac[[13]] <-  glm.nb.SPEdd_stems <-
  glm.nb(nSeqs ~ Raccoon.EDD + log10(Num_Stems) + offset(log(Deploy.Duration)), data = racDataSpr)

ModListSprRac[[14]] <- glm.nb.SPEdd_oak <-
  glm.nb(nSeqs ~ Raccoon.EDD + OakDBH + offset(log(Deploy.Duration)), data = racDataSpr)

ModListSprRac[[15]] <- glm.nb.SPstems_hgt <-
  glm.nb(nSeqs ~ Height_cm + log10(Num_Stems) + offset(log(Deploy.Duration)), data = racDataSpr)

ModListSprRac[[16]] <- glm.nb.SPstems_oak <-
  glm.nb(nSeqs ~ OakDBH + log10(Num_Stems) + offset(log(Deploy.Duration)), data = racDataSpr)

ModListSprRac[[17]] <- glm.nb.SPhgt_oak <-
  glm.nb(nSeqs ~ Height_cm + OakDBH + offset(log(Deploy.Duration)), data = racDataSpr)

ModListSprRac[[18]] <- glm.nb.SPlogStemsHgt <-
  glm.nb(
    nSeqs ~ Log.in.View + Height_cm + log10(Num_Stems) +
      offset(log(Deploy.Duration)),
    data = racDataSpr
  )

ModListSprRac[[19]] <- glm.nb.SPlogStemsOak <-
  glm.nb(
    nSeqs ~ Log.in.View + OakDBH + log10(Num_Stems) +
      offset(log(Deploy.Duration)),
    data = racDataSpr
  )

ModListSprRac[[20]] <- glm.nb.SPlogStemsEdd <-
  glm.nb(
    nSeqs ~ Log.in.View + Raccoon.EDD + log10(Num_Stems) +
      offset(log(Deploy.Duration)),
    data = racDataSpr
  )

ModListSprRac[[21]] <- glm.nb.SPHgtStemsEdd <-
  glm.nb(
    nSeqs ~ Height_cm + Raccoon.EDD + log10(Num_Stems) +
      offset(log(Deploy.Duration)),
    data = racDataSpr
  )

ModListSprRac[[22]] <- glm.nb.SPHgtStemsOak <-
  glm.nb(
    nSeqs ~ Height_cm + OakDBH + log10(Num_Stems) +
      offset(log(Deploy.Duration)),
    data = racDataSpr
  )

ModListSprRac[[23]] <- glm.nb.SPHgtEddOak <-
  glm.nb(
    nSeqs ~ Height_cm + OakDBH + Raccoon.EDD +
      offset(log(Deploy.Duration)),
    data = racDataSpr
  )
ModListSprRac[[24]] <- glm.nb.SPHgtEddLog <-
  glm.nb(
    nSeqs ~ Height_cm + Log.in.View + Raccoon.EDD +
      offset(log(Deploy.Duration)),
    data = racDataSpr
  )
ModListSprRac[[25]] <- glm.nb.SPHgtOakLog <-
  glm.nb(
    nSeqs ~ Height_cm + Log.in.View + OakDBH +
      offset(log(Deploy.Duration)),
    data = racDataSpr
  )
modnamesSP <- c("Full","Intercept","Log","Height","Oak","Stems","EDD","Log + Stems","Log + EDD","Log + Oak","Log + Height","EDD + Height","EDD + Stems","EDD + Oak","Stems + Height","Stems + Oak","Height + Oak","Log + Stems + Height","Log + Stems + Oak","Log + Stems + EDD","Height + Stems + EDD","Height + Stems + Oak","Height + EDD + Oak","Height + EDD + Log","Height + Oak + Log")
modtabSprRac <- aictab(cand.set = ModListSprRac, modnames = modnamesSP)
modtabSprRac

summary(glm.nb.SPoak) #In Spring there is some evidence that raccoon are more frequently captured in front of cameras that have smaller total dbh of oaks.

#Explained Deviance of Best Model
1 - glm.nb.SPoak$deviance / glm.nb.SPoak$null.deviance #0.095165
#Explained Deviance of best 3 variable model
1 - glm.nb.SPHgtEddOak$deviance / glm.nb.SPHgtEddOak$null.deviance #0.16375
#Explained Deviance of Full Model
1 - glm.nb.fullSP$deviance / glm.nb.fullSP$null.deviance #0.194


#Response Curve - Oak effects in Spring
summary(racDataSpr$OakDBH)
seq.oak <- seq(min(racDataSpr$OakDBH), max(racDataSpr$OakDBH), length.out = 100)
nd.seq.oak <- data.frame(Deploy.Duration = 61, OakDBH = seq.oak)
pred.oakSP <- predict(glm.nb.SPoak, newdata = nd.seq.oak, se=TRUE, type = "link")
plot(nSeqs ~ OakDBH, data=racDataSpr, pch=16, col=rgb(.75,.25,0,.5), las=1, xlab = "Total Oak DBH", ylab = "# Raccoon Sequences in Spring")
lines(nd.seq.oak$OakDBH, exp(pred.oakSP$fit), col="dark olive green")
lines(nd.seq.oak$OakDBH, exp(pred.oakSP$fit + 1.96*pred.oakSP$se.fit), lty=2, col="dark olive green")
lines(nd.seq.oak$OakDBH, exp(pred.oakSP$fit - 1.96*pred.oakSP$se.fit), lty=2, col="dark olive green")


