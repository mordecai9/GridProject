###Regression Analysis of SCBI Camera Grid Data for Black Bear

require(AICcmodavg)
require(MASS) #for glm.nb function
source("scripts/pairsPannelFunctions.r")


# Data Import and Exploration - Black Bear ---------------------------------------
# We won't be using winter data here since bear really weren't detected in winter
load("data/bearDataSum.RData")
load("data/bearDataFall.RData")
load("data/bearDataWin.RData")
load("data/bearDataSpr.RData")

pairs(bearDataSum[,c(4,7,13,14,15,17)], diag.panel = panel.hist, lower.panel = panel.smooth, upper.panel = panel.cor)
str(bearDataSum)

#Looking at how many cameras picked up bear in each season
barplot(bearDataSum$nSeqs, xlab = "Summer", ylab = "# of Seqs per camera")
barplot(bearDataFall$nSeqs, xlab = "Fall", ylab = "# of Seqs per camera")
barplot(bearDataWin$nSeqs, xlab = "Winter", ylab = "# of Seqs per camera")
barplot(bearDataSpr$nSeqs, xlab = "Spring", ylab = "# of Seqs per camera")

hist(bearDataSum$Num_Stems)
hist(log10(bearDataSum$Num_Stems)) #because one site has so many stems, we need to log10 transform this parameter for all models
hist(bearDataSum$Summer.Fall.EDD)
hist(bearDataSpr$Winter.Spring.EDD)

#Univariate exploratory plots for each of our 5 covariates
par.default <- par(no.readonly = T)
par(mfrow = c(2,2))
plot(nSeqs ~ Summer.Fall.EDD, data = bearDataSum, main = "Summer")
plot(nSeqs ~ Summer.Fall.EDD, data = bearDataFall, main = "Fall")
plot(nSeqs ~ Winter.Spring.EDD, data = bearDataSpr, main = "Spring")

par(mfrow = c(2,2))
plot(nSeqs ~ OakDBH, data = bearDataSum, main = "Summer")
plot(nSeqs ~ OakDBH, data = bearDataFall, main = "Fall")
plot(nSeqs ~ OakDBH, data = bearDataSpr, main = "Spring")

par(mfrow = c(2,2))
plot(nSeqs ~ log10(Num_Stems), data = bearDataSum, main = "Summer")
plot(nSeqs ~ log10(Num_Stems), data = bearDataFall, main = "Fall")
plot(nSeqs ~ log10(Num_Stems), data = bearDataSpr, main = "Spring")

par(mfrow = c(2,2))
plot(nSeqs ~ Height_cm, data = bearDataSum, main = "Summer")
plot(nSeqs ~ Height_cm, data = bearDataFall, main = "Fall")
plot(nSeqs ~ Height_cm, data = bearDataSpr, main = "Spring")

par(mfrow = c(2,2))
plot(nSeqs ~ Log.in.View, data = bearDataSum, main = "Summer")
plot(nSeqs ~ Log.in.View, data = bearDataFall, main = "Fall")
plot(nSeqs ~ Log.in.View, data = bearDataSpr, main = "Spring")
par(par.default)

hist(bearDataSum$Deploy.Duration, breaks = 5)


# Full Poisson GLM and Overdispersion -------------------------------------

#Exploratory Poisson model for Summer data to see if a negative binomial model is needed due to potential overdispersion
glm.po.fullS <- glm(nSeqs ~ Log.in.View + Summer.Fall.EDD + Height_cm + log10(Num_Stems) + OakDBH +offset(log(Deploy.Duration)), data = bearDataSum, family = poisson)
summary(glm.po.fullS)

#test for overdispersion indicates not much overdispersion here at 1.56
sum.po.fullS <- summary(glm.po.fullS)
phi.po.fullS <- sum.po.fullS$deviance/sum.po.fullS$df.residual 

#Exploratory Poisson model for Fall data
glm.po.fullF <- glm(nSeqs ~ Log.in.View + Summer.Fall.EDD + Height_cm + log10(Num_Stems) + OakDBH +offset(log(Deploy.Duration)), data = bearDataFall, family = poisson)
summary(glm.po.fullF)

#test for overdispersion indicates 2.15, seems a bit borderline 
sum.po.fullF <- summary(glm.po.fullF)
phi.po.fullF <- sum.po.fullF$deviance/sum.po.fullF$df.residual 

#Exploratory Poisson model for Spring data
glm.po.fullSp <- glm(nSeqs ~ Log.in.View + Winter.Spring.EDD + Height_cm + log10(Num_Stems) + OakDBH +offset(log(Deploy.Duration)), data = bearDataSpr, family = poisson)
summary(glm.po.fullSp)

#test for overdispersion shows 1.48, not very high.
sum.po.fullSp <- summary(glm.po.fullSp)
phi.po.fullSp <- sum.po.fullSp$deviance/sum.po.fullSp$df.residual 


# Negative Binomial Regression - Black Bear in Summer---------------------------------

ModListSumBear <- list(NA)

ModListSumBear[[1]] <- glm.nb.fullS <-
  glm.nb(
    nSeqs ~ Log.in.View + Summer.Fall.EDD + Height_cm + log10(Num_Stems) + OakDBH +
      offset(log(Deploy.Duration)),
    data = bearDataSum
  )

ModListSumBear[[2]] <- glm.nb.SI <- glm.nb(nSeqs ~ 1 + offset(log(Deploy.Duration)), data = bearDataSum)

ModListSumBear[[3]] <- glm.nb.Slog <-
  glm.nb(nSeqs ~ Log.in.View + offset(log(Deploy.Duration)), data = bearDataSum)

ModListSumBear[[4]] <- glm.nb.Shgt <-
  glm.nb(nSeqs ~ Height_cm + offset(log(Deploy.Duration)), data = bearDataSum)

ModListSumBear[[5]] <- glm.nb.Soak <-
  glm.nb(nSeqs ~ OakDBH + offset(log(Deploy.Duration)), data = bearDataSum)

ModListSumBear[[6]] <- glm.nb.Sstems <-
  glm.nb(nSeqs ~ log10(Num_Stems) + offset(log(Deploy.Duration)), data = bearDataSum)

ModListSumBear[[7]] <- glm.nb.Sedd <-
  glm.nb(nSeqs ~ Summer.Fall.EDD + offset(log(Deploy.Duration)), data = bearDataSum)

ModListSumBear[[8]] <- glm.nb.Slog_stems <-
  glm.nb(nSeqs ~ Log.in.View + log10(Num_Stems) + offset(log(Deploy.Duration)), data = bearDataSum)

ModListSumBear[[9]] <-glm.nb.Slog_Edd <-
  glm.nb(nSeqs ~ Log.in.View + Summer.Fall.EDD + offset(log(Deploy.Duration)), data = bearDataSum)

ModListSumBear[[10]] <- glm.nb.Slog_Oak <-
  glm.nb(nSeqs ~ Log.in.View + OakDBH + offset(log(Deploy.Duration)), data = bearDataSum)

ModListSumBear[[11]] <- glm.nb.Slog_Hgt <-
  glm.nb(nSeqs ~ Log.in.View + Height_cm + offset(log(Deploy.Duration)), data = bearDataSum)

ModListSumBear[[12]] <- glm.nb.SEdd_Hgt <-
  glm.nb(nSeqs ~ Summer.Fall.EDD + Height_cm + offset(log(Deploy.Duration)), data = bearDataSum)

ModListSumBear[[13]] <-  glm.nb.SEdd_stems <-
  glm.nb(nSeqs ~ Summer.Fall.EDD + log10(Num_Stems) + offset(log(Deploy.Duration)), data = bearDataSum)

ModListSumBear[[14]] <- glm.nb.SEdd_oak <-
  glm.nb(nSeqs ~ Summer.Fall.EDD + OakDBH + offset(log(Deploy.Duration)), data = bearDataSum)

ModListSumBear[[15]] <- glm.nb.Sstems_hgt <-
  glm.nb(nSeqs ~ Height_cm + log10(Num_Stems) + offset(log(Deploy.Duration)), data = bearDataSum)

ModListSumBear[[16]] <- glm.nb.Sstems_oak <-
  glm.nb(nSeqs ~ OakDBH + log10(Num_Stems) + offset(log(Deploy.Duration)), data = bearDataSum)

ModListSumBear[[17]] <- glm.nb.Shgt_oak <-
  glm.nb(nSeqs ~ Height_cm + OakDBH + offset(log(Deploy.Duration)), data = bearDataSum)

ModListSumBear[[18]] <- glm.nb.SlogStemsHgt <-
  glm.nb(
    nSeqs ~ Log.in.View + Height_cm + log10(Num_Stems) +
      offset(log(Deploy.Duration)),
    data = bearDataSum
  )

ModListSumBear[[19]] <- glm.nb.SlogStemsOak <-
  glm.nb(
    nSeqs ~ Log.in.View + OakDBH + log10(Num_Stems) +
      offset(log(Deploy.Duration)),
    data = bearDataSum
  )

ModListSumBear[[20]] <- glm.nb.SlogStemsEdd <-
  glm.nb(
    nSeqs ~ Log.in.View + Summer.Fall.EDD + log10(Num_Stems) +
      offset(log(Deploy.Duration)),
    data = bearDataSum
  )

ModListSumBear[[21]] <- glm.nb.SHgtStemsEdd <-
  glm.nb(
    nSeqs ~ Height_cm + Summer.Fall.EDD + log10(Num_Stems) +
      offset(log(Deploy.Duration)),
    data = bearDataSum
  )

ModListSumBear[[22]] <- glm.nb.SHgtStemsOak <-
  glm.nb(
    nSeqs ~ Height_cm + OakDBH + log10(Num_Stems) +
      offset(log(Deploy.Duration)),
    data = bearDataSum
  )

ModListSumBear[[23]] <- glm.nb.SHgtEddOak <-
  glm.nb(
    nSeqs ~ Height_cm + OakDBH + Summer.Fall.EDD +
      offset(log(Deploy.Duration)),
    data = bearDataSum
  )
ModListSumBear[[24]] <- glm.nb.SHgtEddLog <-
  glm.nb(
    nSeqs ~ Height_cm + Log.in.View + Summer.Fall.EDD +
      offset(log(Deploy.Duration)),
    data = bearDataSum
  )
ModListSumBear[[25]] <- glm.nb.SHgtOakLog <-
  glm.nb(
    nSeqs ~ Height_cm + Log.in.View + OakDBH +
      offset(log(Deploy.Duration)),
    data = bearDataSum
  )
modnamesS <- c("Full","Intercept","Log","Height","Oak","Stems","EDD","Log + Stems","Log + EDD","Log + Oak","Log + Height","EDD + Height","EDD + Stems","EDD + Oak","Stems + Height","Stems + Oak","Height + Oak","Log + Stems + Height","Log + Stems + Oak","Log + Stems + EDD","Height + Stems + EDD","Height + Stems + Oak","Height + EDD + Oak","Height + EDD + Log","Height + Oak + Log")
modtabSumBear <- aictab(cand.set = ModListSumBear, modnames = modnamesS)
modtabSumBear

summary(glm.nb.Slog) #This is the best model. Indicates bear are captured less frequently when there is a log in view, in summer.

#Summed Model Weights for Squirrels in Summer from full table. 
SumSW_edd <- sum(modtabSumBear$AICcWt[grep("EDD", modtabSumBear$Modnames)])
SumSW_OakDBH <- sum(modtabSumBear$AICcWt[grep("Oak", modtabSumBear$Modnames)])
SumSW_Stems <- sum(modtabSumBear$AICcWt[grep("Stems", modtabSumBear$Modnames)])
SumSW_Log <- sum(modtabSumBear$AICcWt[grep("Log", modtabSumBear$Modnames)])
SumSW_Height <- sum(modtabSumBear$AICcWt[grep("Height", modtabSumBear$Modnames)])


#Explained deviance of best model
1 - glm.nb.Slog$deviance / glm.nb.Slog$null.deviance #0.1688
#Explained deviance of full model
1 - glm.nb.fullS$deviance / glm.nb.fullS$null.deviance #0.247


# Bears in Fall NB Full Model Selection -----------------------------------

ModListFallBear <- list(NA)

ModListFallBear[[1]] <- glm.nb.fullF <-
  glm.nb(
    nSeqs ~ Log.in.View + Summer.Fall.EDD + Height_cm + log10(Num_Stems) + OakDBH +
      offset(log(Deploy.Duration)),
    data = bearDataFall
  )

ModListFallBear[[2]] <- glm.nb.FI <- glm.nb(nSeqs ~ 1 + offset(log(Deploy.Duration)), data = bearDataFall)

ModListFallBear[[3]] <- glm.nb.Flog <-
  glm.nb(nSeqs ~ Log.in.View + offset(log(Deploy.Duration)), data = bearDataFall)

ModListFallBear[[4]] <- glm.nb.Fhgt <-
  glm.nb(nSeqs ~ Height_cm + offset(log(Deploy.Duration)), data = bearDataFall)

ModListFallBear[[5]] <- glm.nb.Foak <-
  glm.nb(nSeqs ~ OakDBH + offset(log(Deploy.Duration)), data = bearDataFall)

ModListFallBear[[6]] <- glm.nb.Fstems <-
  glm.nb(nSeqs ~ log10(Num_Stems) + offset(log(Deploy.Duration)), data = bearDataFall)

ModListFallBear[[7]] <- glm.nb.Fedd <-
  glm.nb(nSeqs ~ Summer.Fall.EDD + offset(log(Deploy.Duration)), data = bearDataFall)

ModListFallBear[[8]] <- glm.nb.Flog_stems <-
  glm.nb(nSeqs ~ Log.in.View + log10(Num_Stems) + offset(log(Deploy.Duration)), data = bearDataFall)

ModListFallBear[[9]] <-glm.nb.Flog_Edd <-
  glm.nb(nSeqs ~ Log.in.View + Summer.Fall.EDD + offset(log(Deploy.Duration)), data = bearDataFall)

ModListFallBear[[10]] <- glm.nb.Flog_Oak <-
  glm.nb(nSeqs ~ Log.in.View + OakDBH + offset(log(Deploy.Duration)), data = bearDataFall)

ModListFallBear[[11]] <- glm.nb.Flog_Hgt <-
  glm.nb(nSeqs ~ Log.in.View + Height_cm + offset(log(Deploy.Duration)), data = bearDataFall)

ModListFallBear[[12]] <- glm.nb.FEdd_Hgt <-
  glm.nb(nSeqs ~ Summer.Fall.EDD + Height_cm + offset(log(Deploy.Duration)), data = bearDataFall)

ModListFallBear[[13]] <-  glm.nb.FEdd_stems <-
  glm.nb(nSeqs ~ Summer.Fall.EDD + log10(Num_Stems) + offset(log(Deploy.Duration)), data = bearDataFall)

ModListFallBear[[14]] <- glm.nb.FEdd_oak <-
  glm.nb(nSeqs ~ Summer.Fall.EDD + OakDBH + offset(log(Deploy.Duration)), data = bearDataFall)

ModListFallBear[[15]] <- glm.nb.Fstems_hgt <-
  glm.nb(nSeqs ~ Height_cm + log10(Num_Stems) + offset(log(Deploy.Duration)), data = bearDataFall)

ModListFallBear[[16]] <- glm.nb.Fstems_oak <-
  glm.nb(nSeqs ~ OakDBH + log10(Num_Stems) + offset(log(Deploy.Duration)), data = bearDataFall)

ModListFallBear[[17]] <- glm.nb.Fhgt_oak <-
  glm.nb(nSeqs ~ Height_cm + OakDBH + offset(log(Deploy.Duration)), data = bearDataFall)

ModListFallBear[[18]] <- glm.nb.FlogStemsHgt <-
  glm.nb(
    nSeqs ~ Log.in.View + Height_cm + log10(Num_Stems) +
      offset(log(Deploy.Duration)),
    data = bearDataFall
  )

ModListFallBear[[19]] <- glm.nb.FlogStemsOak <-
  glm.nb(
    nSeqs ~ Log.in.View + OakDBH + log10(Num_Stems) +
      offset(log(Deploy.Duration)),
    data = bearDataFall
  )

ModListFallBear[[20]] <- glm.nb.FlogStemsEdd <-
  glm.nb(
    nSeqs ~ Log.in.View + Summer.Fall.EDD + log10(Num_Stems) +
      offset(log(Deploy.Duration)),
    data = bearDataFall
  )

ModListFallBear[[21]] <- glm.nb.FHgtStemsEdd <-
  glm.nb(
    nSeqs ~ Height_cm + Summer.Fall.EDD + log10(Num_Stems) +
      offset(log(Deploy.Duration)),
    data = bearDataFall
  )

ModListFallBear[[22]] <- glm.nb.FHgtStemsOak <-
  glm.nb(
    nSeqs ~ Height_cm + OakDBH + log10(Num_Stems) +
      offset(log(Deploy.Duration)),
    data = bearDataFall
  )

ModListFallBear[[23]] <- glm.nb.FHgtEddOak <-
  glm.nb(
    nSeqs ~ Height_cm + OakDBH + Summer.Fall.EDD +
      offset(log(Deploy.Duration)),
    data = bearDataFall
  )
ModListFallBear[[24]] <- glm.nb.FHgtEddLog <-
  glm.nb(
    nSeqs ~ Height_cm + Log.in.View + Summer.Fall.EDD +
      offset(log(Deploy.Duration)),
    data = bearDataFall
  )
ModListFallBear[[25]] <- glm.nb.FHgtOakLog <-
  glm.nb(
    nSeqs ~ Height_cm + Log.in.View + OakDBH +
      offset(log(Deploy.Duration)),
    data = bearDataFall
  )
modnamesF <- c("Full","Intercept","Log","Height","Oak","Stems","EDD","Log + Stems","Log + EDD","Log + Oak","Log + Height","EDD + Height","EDD + Stems","EDD + Oak","Stems + Height","Stems + Oak","Height + Oak","Log + Stems + Height","Log + Stems + Oak","Log + Stems + EDD","Height + Stems + EDD","Height + Stems + Oak","Height + EDD + Oak","Height + EDD + Log","Height + Oak + Log")
modtabFallBear <- aictab(cand.set = ModListFallBear, modnames = modnamesF)
modtabFallBear

summary(glm.nb.FEdd_oak) # Seems like best model, indicates bear more likely to be photographed in front of cameras with more/bigger oaks in front of it, in the Fall

#Summed Model Weights Bears in Fall from full table. Note interaction not included at the moment
FallSW_edd <- sum(modtabFallBear$AICcWt[grep("EDD", modtabFallBear$Modnames)])
FallSW_OakDBH <- sum(modtabFallBear$AICcWt[grep("Oak", modtabFallBear$Modnames)])
FallSW_Stems <- sum(modtabFallBear$AICcWt[grep("Stems", modtabFallBear$Modnames)])
FallSW_Log <- sum(modtabFallBear$AICcWt[grep("Log", modtabFallBear$Modnames)])
FallSW_Height <- sum(modtabFallBear$AICcWt[grep("Height", modtabFallBear$Modnames)])

#Explained deviance of best model
1 - glm.nb.FEdd_oak$deviance / glm.nb.FEdd_oak$null.deviance #0.2167876
#Explained deviance of full model
1 - glm.nb.fullF$deviance / glm.nb.fullF$null.deviance #0.2523777

#Response Curve of Oak DBH
summary(bearDataFall$OakDBH)
seq.oak <- seq(min(bearDataFall$OakDBH), max(bearDataFall$OakDBH), length.out = 100)
nd.seq.oak <- data.frame(Deploy.Duration = 61, OakDBH = seq.oak, Summer.Fall.EDD = median(bearDataFall$Summer.Fall.EDD))
pred.oak <- predict(glm.nb.FEdd_oak, newdata = nd.seq.oak, se=TRUE, type = "link")
plot(nSeqs ~ OakDBH, data=bearDataFall, pch=16, col=rgb(.75,.25,0,.5), las=1, xlab = "Total DBH Oak", ylab = "# Bear Sequences in Fall")
lines(nd.seq.oak$OakDBH, exp(pred.oak$fit), col="dark olive green")
lines(nd.seq.oak$OakDBH, exp(pred.oak$fit + 1.96*pred.oak$se.fit), lty=2, col="dark olive green")
lines(nd.seq.oak$OakDBH, exp(pred.oak$fit - 1.96*pred.oak$se.fit), lty=2, col="dark olive green")

#Response Curve of EDD
summary(bearDataFall$OakDBH)
seq.edd <- seq(min(bearDataFall$Summer.Fall.EDD), max(bearDataFall$Summer.Fall.EDD), length.out = 100)
nd.seq.edd <- data.frame(Deploy.Duration = 61, OakDBH = median(bearDataFall$OakDBH), Summer.Fall.EDD = seq.edd)
pred.edd <- predict(glm.nb.FEdd_oak, newdata = nd.seq.edd, se=TRUE, type = "link")
plot(nSeqs ~ Summer.Fall.EDD, data=bearDataFall, pch=16, col=rgb(.75,.25,0,.5), las=1, xlab = "Effective Detection Distance (m)", ylab = "# Bear Sequences in Fall")
lines(nd.seq.edd$Summer.Fall.EDD, exp(pred.edd$fit), col="dark olive green")
lines(nd.seq.edd$Summer.Fall.EDD, exp(pred.edd$fit + 1.96*pred.edd$se.fit), lty=2, col="dark olive green")
lines(nd.seq.edd$Summer.Fall.EDD, exp(pred.edd$fit - 1.96*pred.edd$se.fit), lty=2, col="dark olive green")


# Bear in Spring NB Full Model Selection ----------------------------------

ModListSprBear <- list(NA)
ModListSprBear[[1]] <- glm.nb.fullSP <-
  glm.nb(
    nSeqs ~ Log.in.View + Winter.Spring.EDD + Height_cm + log10(Num_Stems) + OakDBH +
      offset(log(Deploy.Duration)),
    data = bearDataSpr
  )

ModListSprBear[[2]] <- glm.nb.SPI <- glm.nb(nSeqs ~ 1 + offset(log(Deploy.Duration)), data = bearDataSpr)

ModListSprBear[[3]] <- glm.nb.SPlog <-
  glm.nb(nSeqs ~ Log.in.View + offset(log(Deploy.Duration)), data = bearDataSpr)

ModListSprBear[[4]] <- glm.nb.SPhgt <-
  glm.nb(nSeqs ~ Height_cm + offset(log(Deploy.Duration)), data = bearDataSpr)

ModListSprBear[[5]] <- glm.nb.SPoak <-
  glm.nb(nSeqs ~ OakDBH + offset(log(Deploy.Duration)), data = bearDataSpr)

ModListSprBear[[6]] <- glm.nb.SPstems <-
  glm.nb(nSeqs ~ log10(Num_Stems) + offset(log(Deploy.Duration)), data = bearDataSpr)

ModListSprBear[[7]] <- glm.nb.SPedd <-
  glm.nb(nSeqs ~ Winter.Spring.EDD + offset(log(Deploy.Duration)), data = bearDataSpr)

ModListSprBear[[8]] <- glm.nb.SPlog_stems <-
  glm.nb(nSeqs ~ Log.in.View + log10(Num_Stems) + offset(log(Deploy.Duration)), data = bearDataSpr)

ModListSprBear[[9]] <-glm.nb.SPlog_Edd <-
  glm.nb(nSeqs ~ Log.in.View + Winter.Spring.EDD + offset(log(Deploy.Duration)), data = bearDataSpr)

ModListSprBear[[10]] <- glm.nb.SPlog_Oak <-
  glm.nb(nSeqs ~ Log.in.View + OakDBH + offset(log(Deploy.Duration)), data = bearDataSpr)

ModListSprBear[[11]] <- glm.nb.SPlog_Hgt <-
  glm.nb(nSeqs ~ Log.in.View + Height_cm + offset(log(Deploy.Duration)), data = bearDataSpr)

ModListSprBear[[12]] <- glm.nb.SPEdd_Hgt <-
  glm.nb(nSeqs ~ Winter.Spring.EDD + Height_cm + offset(log(Deploy.Duration)), data = bearDataSpr)

ModListSprBear[[13]] <-  glm.nb.SPEdd_stems <-
  glm.nb(nSeqs ~ Winter.Spring.EDD + log10(Num_Stems) + offset(log(Deploy.Duration)), data = bearDataSpr)

ModListSprBear[[14]] <- glm.nb.SPEdd_oak <-
  glm.nb(nSeqs ~ Winter.Spring.EDD + OakDBH + offset(log(Deploy.Duration)), data = bearDataSpr)

ModListSprBear[[15]] <- glm.nb.SPstems_hgt <-
  glm.nb(nSeqs ~ Height_cm + log10(Num_Stems) + offset(log(Deploy.Duration)), data = bearDataSpr)

ModListSprBear[[16]] <- glm.nb.SPstems_oak <-
  glm.nb(nSeqs ~ OakDBH + log10(Num_Stems) + offset(log(Deploy.Duration)), data = bearDataSpr)

ModListSprBear[[17]] <- glm.nb.SPhgt_oak <-
  glm.nb(nSeqs ~ Height_cm + OakDBH + offset(log(Deploy.Duration)), data = bearDataSpr)

ModListSprBear[[18]] <- glm.nb.SPlogStemsHgt <-
  glm.nb(
    nSeqs ~ Log.in.View + Height_cm + log10(Num_Stems) +
      offset(log(Deploy.Duration)),
    data = bearDataSpr
  )

ModListSprBear[[19]] <- glm.nb.SPlogStemsOak <-
  glm.nb(
    nSeqs ~ Log.in.View + OakDBH + log10(Num_Stems) +
      offset(log(Deploy.Duration)),
    data = bearDataSpr
  )

ModListSprBear[[20]] <- glm.nb.SPlogStemsEdd <-
  glm.nb(
    nSeqs ~ Log.in.View + Winter.Spring.EDD + log10(Num_Stems) +
      offset(log(Deploy.Duration)),
    data = bearDataSpr
  )

ModListSprBear[[21]] <- glm.nb.SPHgtStemsEdd <-
  glm.nb(
    nSeqs ~ Height_cm + Winter.Spring.EDD + log10(Num_Stems) +
      offset(log(Deploy.Duration)),
    data = bearDataSpr
  )

ModListSprBear[[22]] <- glm.nb.SPHgtStemsOak <-
  glm.nb(
    nSeqs ~ Height_cm + OakDBH + log10(Num_Stems) +
      offset(log(Deploy.Duration)),
    data = bearDataSpr
  )

ModListSprBear[[23]] <- glm.nb.SPHgtEddOak <-
  glm.nb(
    nSeqs ~ Height_cm + OakDBH + Winter.Spring.EDD +
      offset(log(Deploy.Duration)),
    data = bearDataSpr
  )
ModListSprBear[[24]] <- glm.nb.SPHgtEddLog <-
  glm.nb(
    nSeqs ~ Height_cm + Log.in.View + Winter.Spring.EDD +
      offset(log(Deploy.Duration)),
    data = bearDataSpr
  )
ModListSprBear[[25]] <- glm.nb.SPHgtOakLog <-
  glm.nb(
    nSeqs ~ Height_cm + Log.in.View + OakDBH +
      offset(log(Deploy.Duration)),
    data = bearDataSpr
  )
modnamesSP <- c("Full","Intercept","Log","Height","Oak","Stems","EDD","Log + Stems","Log + EDD","Log + Oak","Log + Height","EDD + Height","EDD + Stems","EDD + Oak","Stems + Height","Stems + Oak","Height + Oak","Log + Stems + Height","Log + Stems + Oak","Log + Stems + EDD","Height + Stems + EDD","Height + Stems + Oak","Height + EDD + Oak","Height + EDD + Log","Height + Oak + Log")
modtabSprBear <- aictab(cand.set = ModListSprBear, modnames = modnamesSP)
modtabSprBear

#Summed Model Weights for Bears in Spring from full table. 
SprSW_edd <- sum(modtabSprBear$AICcWt[grep("EDD", modtabSprBear$Modnames)])
SprSW_OakDBH <- sum(modtabSprBear$AICcWt[grep("Oak", modtabSprBear$Modnames)])
SprSW_Stems <- sum(modtabSprBear$AICcWt[grep("Stems", modtabSprBear$Modnames)])
SprSW_Log <- sum(modtabSprBear$AICcWt[grep("Log", modtabSprBear$Modnames)])
SprSW_Height <- sum(modtabSprBear$AICcWt[grep("Height", modtabSprBear$Modnames)])

#Nothing useful in Spring. Intercept is best model.
#Explained deviance of best model
1 - glm.nb.SPI$deviance / glm.nb.SPI$null.deviance #0.0
#Explained deviance of full model
1 - glm.nb.fullSP$deviance / glm.nb.fullSP$null.deviance #0.1462535