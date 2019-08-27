###Regression Analysis of SCBI Camera Grid Data
rm(list = ls())
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

summary(deerDataSum) #no missing values
hist(deerDataSum$Num_Stems)
summary(deerDataSum$Num_Stems) #no zeros
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

#Looking at how many cameras picked up deer in each season, and how often
barplot(deerDataSum$nSeqs, xlab = "Summer", ylab = "# of Seqs per camera")
barplot(deerDataFall$nSeqs, xlab = "Fall", ylab = "# of Seqs per camera")
barplot(deerDataWin$nSeqs, xlab = "Winter", ylab = "# of Seqs per camera")
barplot(deerDataSpr$nSeqs, xlab = "Spring", ylab = "# of Seqs per camera")

# Full Poisson GLM and Overdispersion -------------------------------------

#Exploratory Poisson model for Summer data to see if a negative binomial model is needed due to potential overdispersion. 5 covariates for 27 data points is too much, but just checking for overdispersion
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
# Summer Deer NB Full Model Comparison ----------------------------------------------

#Here I am running 27 models. All the combinations of our 5 parameters, with a maximum of three variables in each model, plus the full model. But I'm not considering the full model in the AIC table. Here I manually created the AIC table, because the aictab function doesn't do work on glm.nb at this point. Afterward they updated the package so it works with nb! So I don't need all this crap in the end.


glm.nb.fullS <-
  glm.nb(
    nSeqs ~ Log.in.View + Summer.Fall.EDD + Height_cm + log10(Num_Stems) + OakDBH +
      offset(log(Deploy.Duration)),
    data = deerDataSum
  )
summary(glm.nb.fullS)
AICc(glm.nb.fullS)
logLik(glm.nb.fullS)

#Full 26 model AIC table comparison

modtabSumDeer <- data.frame(matrix(NA, nrow = 26, ncol = 8))
names(modtabSumDeer) <- c("Model_Name", "K", "AIC", "AICc", "DeltaAICc", "ModelL", "ModelW", "LL")


ModListSumDeer <- list(NA)
i = 1

modtabSumDeer$Model_Name[i] <- "Intercept"
glm.nb.SI <- glm.nb(nSeqs ~ 1 + offset(log(Deploy.Duration)), data = deerDataSum)

ModListSumDeer[[i]] <- summary(glm.nb.SI)
modtabSumDeer$AICc[i] <- AICc(glm.nb.SI)
modtabSumDeer$LL[i] <- logLik(glm.nb.SI)
i = i+1

modtabSumDeer$Model_Name[i] <- "Log"
glm.nb.Slog <-
  glm.nb(nSeqs ~ Log.in.View + offset(log(Deploy.Duration)), data = deerDataSum)
ModListSumDeer[[i]] <- summary(glm.nb.Slog)
modtabSumDeer$AICc[i] <- AICc(glm.nb.Slog)
modtabSumDeer$LL[i] <- logLik(glm.nb.Slog)
i = i+1

modtabSumDeer$Model_Name[i] <- "Height"

glm.nb.Shgt <-
  glm.nb(nSeqs ~ Height_cm + offset(log(Deploy.Duration)), data = deerDataSum)
ModListSumDeer[[i]] <- summary(glm.nb.Shgt)
modtabSumDeer$AICc[i] <- AICc(glm.nb.Shgt)
modtabSumDeer$LL[i] <- logLik(glm.nb.Shgt)
i = i+1

modtabSumDeer$Model_Name[i] <- "Oak"
glm.nb.Soak <-
  glm.nb(nSeqs ~ OakDBH + offset(log(Deploy.Duration)), data = deerDataSum)
ModListSumDeer[[i]] <- summary(glm.nb.Soak)
modtabSumDeer$AICc[i] <- AICc(glm.nb.Soak)
modtabSumDeer$LL[i] <- logLik(glm.nb.Soak)
i = i+1

modtabSumDeer$Model_Name[i] <- "Stems"
glm.nb.Sstems <-
  glm.nb(nSeqs ~ log10(Num_Stems) + offset(log(Deploy.Duration)), data = deerDataSum)
ModListSumDeer[[i]] <- summary(glm.nb.Sstems)
modtabSumDeer$AICc[i] <- AICc(glm.nb.Sstems)
modtabSumDeer$LL[i] <- logLik(glm.nb.Sstems)
i = i+1

modtabSumDeer$Model_Name[i] <- "EDD"
glm.nb.Sedd <-
  glm.nb(nSeqs ~ Summer.Fall.EDD + offset(log(Deploy.Duration)), data = deerDataSum)
ModListSumDeer[[i]] <- summary(glm.nb.Sedd)
modtabSumDeer$AICc[i] <- AICc(glm.nb.Sedd)
modtabSumDeer$LL[i] <- logLik(glm.nb.Sedd)
i = i+1

modtabSumDeer$Model_Name[i] <- "Log + Stems"
glm.nb.Slog_stems <-
  glm.nb(nSeqs ~ Log.in.View + log10(Num_Stems) + offset(log(Deploy.Duration)), data = deerDataSum)
ModListSumDeer[[i]] <- summary(glm.nb.Slog_stems)
modtabSumDeer$AICc[i] <- AICc(glm.nb.Slog_stems)
modtabSumDeer$LL[i] <- logLik(glm.nb.Slog_stems)
i = i+1

modtabSumDeer$Model_Name[i] <- "Log + EDD"
glm.nb.Slog_Edd <-
  glm.nb(nSeqs ~ Log.in.View + Summer.Fall.EDD + offset(log(Deploy.Duration)), data = deerDataSum)
ModListSumDeer[[i]] <- summary(glm.nb.Slog_Edd)
modtabSumDeer$AICc[i] <- AICc(glm.nb.Slog_Edd)
modtabSumDeer$LL[i] <- logLik(glm.nb.Slog_Edd)
i = i+1

modtabSumDeer$Model_Name[i] <- "Log + Oak"
glm.nb.Slog_Oak <-
  glm.nb(nSeqs ~ Log.in.View + OakDBH + offset(log(Deploy.Duration)), data = deerDataSum)
ModListSumDeer[[i]] <- summary(glm.nb.Slog_Oak)
modtabSumDeer$AICc[i] <- AICc(glm.nb.Slog_Oak)
modtabSumDeer$LL[i] <- logLik(glm.nb.Slog_Oak)
i = i+1

modtabSumDeer$Model_Name[i] <- "Log + Height"
glm.nb.Slog_Hgt <-
  glm.nb(nSeqs ~ Log.in.View + Height_cm + offset(log(Deploy.Duration)), data = deerDataSum)
ModListSumDeer[[i]] <- summary(glm.nb.Slog_Hgt)
modtabSumDeer$AICc[i] <- AICc(glm.nb.Slog_Hgt)
modtabSumDeer$LL[i] <- logLik(glm.nb.Slog_Hgt)
i = i+1

modtabSumDeer$Model_Name[i] <- "EDD + Height"
glm.nb.SEdd_Hgt <-
  glm.nb(nSeqs ~ Summer.Fall.EDD + Height_cm + offset(log(Deploy.Duration)), data = deerDataSum)
ModListSumDeer[[i]] <- summary(glm.nb.SEdd_Hgt)
modtabSumDeer$AICc[i] <- AICc(glm.nb.SEdd_Hgt)
modtabSumDeer$LL[i] <- logLik(glm.nb.SEdd_Hgt)
i = i+1

modtabSumDeer$Model_Name[i] <- "EDD + Stems"
glm.nb.SEdd_stems <-
  glm.nb(nSeqs ~ Summer.Fall.EDD + log10(Num_Stems) + offset(log(Deploy.Duration)), data = deerDataSum)
ModListSumDeer[[i]] <- summary(glm.nb.SEdd_Hgt)
modtabSumDeer$AICc[i] <- AICc(glm.nb.SEdd_Hgt)
modtabSumDeer$LL[i] <- logLik(glm.nb.SEdd_Hgt)
i = i+1

modtabSumDeer$Model_Name[i] <- "EDD + Oak"
glm.nb.SEdd_oak <-
  glm.nb(nSeqs ~ Summer.Fall.EDD + OakDBH + offset(log(Deploy.Duration)), data = deerDataSum)
ModListSumDeer[[i]] <- summary(glm.nb.SEdd_oak)
modtabSumDeer$AICc[i] <- AICc(glm.nb.SEdd_oak)
modtabSumDeer$LL[i] <- logLik(glm.nb.SEdd_oak)
i = i+1

modtabSumDeer$Model_Name[i] <- "Stems + Height"
glm.nb.Sstems_hgt <-
  glm.nb(nSeqs ~ Height_cm + log10(Num_Stems) + offset(log(Deploy.Duration)), data = deerDataSum)
ModListSumDeer[[i]] <- summary(glm.nb.Sstems_hgt)
modtabSumDeer$AICc[i] <- AICc(glm.nb.Sstems_hgt)
modtabSumDeer$LL[i] <- logLik(glm.nb.Sstems_hgt)
i = i+1

modtabSumDeer$Model_Name[i] <- "Stems + Oak"
glm.nb.Sstems_oak <-
  glm.nb(nSeqs ~ OakDBH + log10(Num_Stems) + offset(log(Deploy.Duration)), data = deerDataSum)
ModListSumDeer[[i]] <- summary(glm.nb.Sstems_oak)
modtabSumDeer$AICc[i] <- AICc(glm.nb.Sstems_oak)
modtabSumDeer$LL[i] <- logLik(glm.nb.Sstems_oak)
i = i+1

modtabSumDeer$Model_Name[i] <- "Height + Oak"
glm.nb.Shgt_oak <-
  glm.nb(nSeqs ~ Height_cm + OakDBH + offset(log(Deploy.Duration)), data = deerDataSum)
ModListSumDeer[[i]] <- summary(glm.nb.Shgt_oak)
modtabSumDeer$AICc[i] <- AICc(glm.nb.Shgt_oak)
modtabSumDeer$LL[i] <- logLik(glm.nb.Shgt_oak)
i = i+1


modtabSumDeer$Model_Name[i] <- "Log + Stems + Height"

glm.nb.SlogStemsHgt <-
  glm.nb(
    nSeqs ~ Log.in.View + Height_cm + log10(Num_Stems) +
      offset(log(Deploy.Duration)),
    data = deerDataSum
  )
ModListSumDeer[[i]] <- summary(glm.nb.SlogStemsHgt)
modtabSumDeer$AICc[i] <- AICc(glm.nb.SlogStemsHgt)
modtabSumDeer$LL[i] <- logLik(glm.nb.SlogStemsHgt)
i = i+1


modtabSumDeer$Model_Name[i] <- "Log + Stems + Oak"
glm.nb.SlogStemsOak <-
  glm.nb(
    nSeqs ~ Log.in.View + OakDBH + log10(Num_Stems) +
      offset(log(Deploy.Duration)),
    data = deerDataSum
  )
ModListSumDeer[[i]] <- summary(glm.nb.SlogStemsOak)
modtabSumDeer$AICc[i] <- AICc(glm.nb.SlogStemsOak)
modtabSumDeer$LL[i] <- logLik(glm.nb.SlogStemsOak)
i = i+1

modtabSumDeer$Model_Name[i] <- "Log + Stems + EDD"
glm.nb.SlogStemsEdd <-
  glm.nb(
    nSeqs ~ Log.in.View + Summer.Fall.EDD + log10(Num_Stems) +
      offset(log(Deploy.Duration)),
    data = deerDataSum
  )
ModListSumDeer[[i]] <- summary(glm.nb.SlogStemsEdd)
modtabSumDeer$AICc[i] <- AICc(glm.nb.SlogStemsEdd)
modtabSumDeer$LL[i] <- logLik(glm.nb.SlogStemsEdd)
i = i+1

modtabSumDeer$Model_Name[i] <- "Height + Stems + EDD"
glm.nb.SHgtStemsEdd <-
  glm.nb(
    nSeqs ~ Height_cm + Summer.Fall.EDD + log10(Num_Stems) +
      offset(log(Deploy.Duration)),
    data = deerDataSum
  )
ModListSumDeer[[i]] <- summary(glm.nb.SHgtStemsEdd)
modtabSumDeer$AICc[i] <- AICc(glm.nb.SHgtStemsEdd)
modtabSumDeer$LL[i] <- logLik(glm.nb.SHgtStemsEdd)
i = i+1

modtabSumDeer$Model_Name[i] <- "Height + Stems + Oak"
glm.nb.SHgtStemsOak <-
  glm.nb(
    nSeqs ~ Height_cm + OakDBH + log10(Num_Stems) +
      offset(log(Deploy.Duration)),
    data = deerDataSum
  )
ModListSumDeer[[i]] <- summary(glm.nb.SHgtStemsOak)
modtabSumDeer$AICc[i] <- AICc(glm.nb.SHgtStemsOak)
modtabSumDeer$LL[i] <- logLik(glm.nb.SHgtStemsOak)
i = i+1

modtabSumDeer$Model_Name[i] <- "Height + EDD + Oak"
glm.nb.SHgtEddOak <-
  glm.nb(
    nSeqs ~ Height_cm + OakDBH + Summer.Fall.EDD +
      offset(log(Deploy.Duration)),
    data = deerDataSum
  )
ModListSumDeer[[i]] <- summary(glm.nb.SHgtEddOak)
modtabSumDeer$AICc[i] <- AICc(glm.nb.SHgtEddOak)
modtabSumDeer$LL[i] <- logLik(glm.nb.SHgtEddOak)
i = i+1

modtabSumDeer$Model_Name[i] <- "Height + EDD + Log"
glm.nb.SHgtEddLog <-
  glm.nb(
    nSeqs ~ Height_cm + Log.in.View + Summer.Fall.EDD +
      offset(log(Deploy.Duration)),
    data = deerDataSum
  )
ModListSumDeer[[i]] <- summary(glm.nb.SHgtEddLog)
modtabSumDeer$AICc[i] <- AICc(glm.nb.SHgtEddLog)
modtabSumDeer$LL[i] <- logLik(glm.nb.SHgtEddLog)
i = i+1

modtabSumDeer$Model_Name[i] <- "Height + Oak + Log"
glm.nb.SHgtOakLog <-
  glm.nb(
    nSeqs ~ Height_cm + Log.in.View + OakDBH +
      offset(log(Deploy.Duration)),
    data = deerDataSum
  )
ModListSumDeer[[i]] <- summary(glm.nb.SHgtOakLog)
modtabSumDeer$AICc[i] <- AICc(glm.nb.SHgtOakLog)
modtabSumDeer$LL[i] <- logLik(glm.nb.SHgtOakLog)
i = i+1

modtabSumDeer$Model_Name[i] <- "Stems + Oak + EDD"
glm.nb.SOakStemsEdd <-
  glm.nb(
    nSeqs ~ log10(Num_Stems) + Summer.Fall.EDD + OakDBH +
      offset(log(Deploy.Duration)),
    data = deerDataSum
  )
ModListSumDeer[[i]] <- summary(glm.nb.SOakStemsEdd)
modtabSumDeer$AICc[i] <- AICc(glm.nb.SOakStemsEdd)
modtabSumDeer$LL[i] <- logLik(glm.nb.SOakStemsEdd)
i = i+1

modtabSumDeer$Model_Name[i] <- "Log + Oak + EDD"
glm.nb.SLogOakEdd <-
  glm.nb(
    nSeqs ~ Log.in.View + Summer.Fall.EDD + OakDBH +
      offset(log(Deploy.Duration)),
    data = deerDataSum
  )
ModListSumDeer[[i]] <- summary(glm.nb.SLogOakEdd)
modtabSumDeer$AICc[i] <- AICc(glm.nb.SLogOakEdd)
modtabSumDeer$LL[i] <- logLik(glm.nb.SLogOakEdd)



for (i in 1:26) {
  modtabSumDeer$AIC[i] <- ModListSumDeer[[i]]$aic
  modtabSumDeer$K[i] <- length(ModListSumDeer[[i]]$coefficients[,1])
  
}

#Calculate rest of table, based on AICc values
modtabSumDeer <- modtabSumDeer[order(modtabSumDeer$AICc),]
modtabSumDeer$DeltaAICc <- modtabSumDeer$AICc - min(modtabSumDeer$AICc)
modtabSumDeer$ModelL <- exp(-0.5*(modtabSumDeer$DeltaAICc))  
modtabSumDeer$ModelW <- modtabSumDeer$ModelL/sum(modtabSumDeer$ModelL)
modtabSumDeer
#write.csv(modtabSumDeer, "DeerSummerModTable.csv", row.names = F)

#The best model has only log as a covariate. Presence of a log reduces # of sequences for deer
summary(glm.nb.Slog)

#Explained Deviance of best model
1 - glm.nb.Slog$deviance / glm.nb.Slog$null.deviance #0.3784
#Explained Deviance of best 3 variable model
1 - glm.nb.SlogStemsEdd$deviance / glm.nb.SlogStemsEdd$null.deviance #0.430
#Explained Deviance of full model
1 - glm.nb.fullS$deviance / glm.nb.fullS$null.deviance #0.454 


#Calculate summed model weights for the 5 different covariates
SumSW_edd <- sum(modtabSumDeer$ModelW[grep("EDD", modtabSumDeer$Model_Name)])
SumSW_OakDBH <- sum(modtabSumDeer$ModelW[grep("Oak", modtabSumDeer$Model_Name)])
SumSW_Stems <- sum(modtabSumDeer$ModelW[grep("Stems", modtabSumDeer$Model_Name)])
SumSW_Log <- sum(modtabSumDeer$ModelW[grep("Log", modtabSumDeer$Model_Name)])
SumSW_Height <- sum(modtabSumDeer$ModelW[grep("Height", modtabSumDeer$Model_Name)])


# Fall Deer NB Full Model Comparison --------------------------------------

glm.nb.fullF <-
  glm.nb(
    nSeqs ~ Log.in.View + Summer.Fall.EDD + Height_cm + log10(Num_Stems) + OakDBH +
      offset(log(Deploy.Duration)),
    data = deerDataFall
  )

ModListFallDeer <- list(NA)

ModListFallDeer[[1]] <- glm.nb.FI <- glm.nb(nSeqs ~ 1 + offset(log(Deploy.Duration)), data = deerDataFall)

ModListFallDeer[[2]] <- glm.nb.Flog <-
  glm.nb(nSeqs ~ Log.in.View + offset(log(Deploy.Duration)), data = deerDataFall)

ModListFallDeer[[3]] <- glm.nb.Fhgt <-
  glm.nb(nSeqs ~ Height_cm + offset(log(Deploy.Duration)), data = deerDataFall)

ModListFallDeer[[4]] <- glm.nb.Foak <-
  glm.nb(nSeqs ~ OakDBH + offset(log(Deploy.Duration)), data = deerDataFall)

ModListFallDeer[[5]] <- glm.nb.Fstems <-
  glm.nb(nSeqs ~ log10(Num_Stems) + offset(log(Deploy.Duration)), data = deerDataFall)

ModListFallDeer[[6]] <- glm.nb.Fedd <-
  glm.nb(nSeqs ~ Summer.Fall.EDD + offset(log(Deploy.Duration)), data = deerDataFall)

ModListFallDeer[[7]] <- glm.nb.Flog_stems <-
  glm.nb(nSeqs ~ Log.in.View + log10(Num_Stems) + offset(log(Deploy.Duration)), data = deerDataFall)

ModListFallDeer[[8]] <-glm.nb.Flog_Edd <-
  glm.nb(nSeqs ~ Log.in.View + Summer.Fall.EDD + offset(log(Deploy.Duration)), data = deerDataFall)

ModListFallDeer[[9]] <- glm.nb.Flog_Oak <-
  glm.nb(nSeqs ~ Log.in.View + OakDBH + offset(log(Deploy.Duration)), data = deerDataFall)

ModListFallDeer[[10]] <- glm.nb.Flog_Hgt <-
  glm.nb(nSeqs ~ Log.in.View + Height_cm + offset(log(Deploy.Duration)), data = deerDataFall)

ModListFallDeer[[11]] <- glm.nb.FEdd_Hgt <-
  glm.nb(nSeqs ~ Summer.Fall.EDD + Height_cm + offset(log(Deploy.Duration)), data = deerDataFall)

ModListFallDeer[[12]] <-  glm.nb.FEdd_stems <-
  glm.nb(nSeqs ~ Summer.Fall.EDD + log10(Num_Stems) + offset(log(Deploy.Duration)), data = deerDataFall)

ModListFallDeer[[13]] <- glm.nb.FEdd_oak <-
  glm.nb(nSeqs ~ Summer.Fall.EDD + OakDBH + offset(log(Deploy.Duration)), data = deerDataFall)

ModListFallDeer[[14]] <- glm.nb.Fstems_hgt <-
  glm.nb(nSeqs ~ Height_cm + log10(Num_Stems) + offset(log(Deploy.Duration)), data = deerDataFall)

ModListFallDeer[[15]] <- glm.nb.Fstems_oak <-
  glm.nb(nSeqs ~ OakDBH + log10(Num_Stems) + offset(log(Deploy.Duration)), data = deerDataFall)

ModListFallDeer[[16]] <- glm.nb.Fhgt_oak <-
  glm.nb(nSeqs ~ Height_cm + OakDBH + offset(log(Deploy.Duration)), data = deerDataFall)

ModListFallDeer[[17]] <- glm.nb.FlogStemsHgt <-
  glm.nb(
    nSeqs ~ Log.in.View + Height_cm + log10(Num_Stems) +
      offset(log(Deploy.Duration)),
    data = deerDataFall
  )

ModListFallDeer[[18]] <- glm.nb.FlogStemsOak <-
  glm.nb(
    nSeqs ~ Log.in.View + OakDBH + log10(Num_Stems) +
      offset(log(Deploy.Duration)),
    data = deerDataFall
  )

ModListFallDeer[[19]] <- glm.nb.FlogStemsEdd <-
  glm.nb(
    nSeqs ~ Log.in.View + Summer.Fall.EDD + log10(Num_Stems) +
      offset(log(Deploy.Duration)),
    data = deerDataFall
  )

ModListFallDeer[[20]] <- glm.nb.FHgtStemsEdd <-
  glm.nb(
    nSeqs ~ Height_cm + Summer.Fall.EDD + log10(Num_Stems) +
      offset(log(Deploy.Duration)),
    data = deerDataFall
  )

ModListFallDeer[[21]] <- glm.nb.FHgtStemsOak <-
  glm.nb(
    nSeqs ~ Height_cm + OakDBH + log10(Num_Stems) +
      offset(log(Deploy.Duration)),
    data = deerDataFall
  )

ModListFallDeer[[22]] <- glm.nb.FHgtEddOak <-
  glm.nb(
    nSeqs ~ Height_cm + OakDBH + Summer.Fall.EDD +
      offset(log(Deploy.Duration)),
    data = deerDataFall
  )
ModListFallDeer[[23]] <- glm.nb.FHgtEddLog <-
  glm.nb(
    nSeqs ~ Height_cm + Log.in.View + Summer.Fall.EDD +
      offset(log(Deploy.Duration)),
    data = deerDataFall
  )
ModListFallDeer[[24]] <- glm.nb.FHgtOakLog <-
  glm.nb(
    nSeqs ~ Height_cm + Log.in.View + OakDBH +
      offset(log(Deploy.Duration)),
    data = deerDataFall
  )

ModListFallDeer[[25]] <- glm.nb.SOakStemsEdd <-
  glm.nb(
    nSeqs ~ log10(Num_Stems) + Summer.Fall.EDD + OakDBH +
      offset(log(Deploy.Duration)),
    data = deerDataFall
  )

ModListFallDeer[[26]] <- glm.nb.SLogOakEdd <-
  glm.nb(
    nSeqs ~ Log.in.View + Summer.Fall.EDD + OakDBH +
      offset(log(Deploy.Duration)),
    data = deerDataFall
  )

modnamesF <- c("Intercept","Log","Height","Oak","Stems","EDD","Log + Stems","Log + EDD","Log + Oak","Log + Height","EDD + Height","EDD + Stems","EDD + Oak","Stems + Height","Stems + Oak","Height + Oak","Log + Stems + Height","Log + Stems + Oak","Log + Stems + EDD","Height + Stems + EDD","Height + Stems + Oak","Height + EDD + Oak","Height + EDD + Log","Height + Oak + Log", "Stems + EDD + Oak", "Log + EDD + Oak")
modtabFallDeer <- aictab(cand.set = ModListFallDeer, modnames = modnamesF)
modtabFallDeer #Intercept is the best model!

#Summed Model Weights for Deer in Fall from Full Table
FallSW_edd <- sum(modtabFallDeer$AICcWt[grep("EDD", modtabFallDeer$Modnames)])
FallSW_OakDBH <- sum(modtabFallDeer$AICcWt[grep("Oak", modtabFallDeer$Modnames)])
FallSW_Stems <- sum(modtabFallDeer$AICcWt[grep("Stems", modtabFallDeer$Modnames)])
FallSW_Log <- sum(modtabFallDeer$AICcWt[grep("Log", modtabFallDeer$Modnames)])
FallSW_Height <- sum(modtabFallDeer$AICcWt[grep("Height", modtabFallDeer$Modnames)])

#Explained Deviance of best model (Intercept)
1 - glm.nb.FI$deviance / glm.nb.FI$null.deviance #0.0
#Explained Deviance of Full Model
1 - glm.nb.fullF$deviance / glm.nb.fullF$null.deviance #0.1172986


# Winter Deer NB Full Model Selection -------------------------------------

glm.nb.fullW <-
  glm.nb(
    nSeqs ~ Log.in.View + Winter.Spring.EDD + Height_cm + log10(Num_Stems) + OakDBH +
      offset(log(Deploy.Duration)),
    data = deerDataWin
  )

ModListWinDeer <- list(NA)


ModListWinDeer[[1]] <- glm.nb.WI <- glm.nb(nSeqs ~ 1 + offset(log(Deploy.Duration)), data = deerDataWin)

ModListWinDeer[[2]] <- glm.nb.Wlog <-
  glm.nb(nSeqs ~ Log.in.View + offset(log(Deploy.Duration)), data = deerDataWin)

ModListWinDeer[[3]] <- glm.nb.Whgt <-
  glm.nb(nSeqs ~ Height_cm + offset(log(Deploy.Duration)), data = deerDataWin)

ModListWinDeer[[4]] <- glm.nb.Woak <-
  glm.nb(nSeqs ~ OakDBH + offset(log(Deploy.Duration)), data = deerDataWin)

ModListWinDeer[[5]] <- glm.nb.Wstems <-
  glm.nb(nSeqs ~ log10(Num_Stems) + offset(log(Deploy.Duration)), data = deerDataWin)

ModListWinDeer[[6]] <- glm.nb.Wedd <-
  glm.nb(nSeqs ~ Winter.Spring.EDD + offset(log(Deploy.Duration)), data = deerDataWin)

ModListWinDeer[[7]] <- glm.nb.Wlog_stems <-
  glm.nb(nSeqs ~ Log.in.View + log10(Num_Stems) + offset(log(Deploy.Duration)), data = deerDataWin)

ModListWinDeer[[8]] <-glm.nb.Wlog_Edd <-
  glm.nb(nSeqs ~ Log.in.View + Winter.Spring.EDD + offset(log(Deploy.Duration)), data = deerDataWin)

ModListWinDeer[[9]] <- glm.nb.Wlog_Oak <-
  glm.nb(nSeqs ~ Log.in.View + OakDBH + offset(log(Deploy.Duration)), data = deerDataWin)

ModListWinDeer[[10]] <- glm.nb.Wlog_Hgt <-
  glm.nb(nSeqs ~ Log.in.View + Height_cm + offset(log(Deploy.Duration)), data = deerDataWin)

ModListWinDeer[[11]] <- glm.nb.WEdd_Hgt <-
  glm.nb(nSeqs ~ Winter.Spring.EDD + Height_cm + offset(log(Deploy.Duration)), data = deerDataWin)

ModListWinDeer[[12]] <-  glm.nb.WEdd_stems <-
  glm.nb(nSeqs ~ Winter.Spring.EDD + log10(Num_Stems) + offset(log(Deploy.Duration)), data = deerDataWin)

ModListWinDeer[[13]] <- glm.nb.WEdd_oak <-
  glm.nb(nSeqs ~ Winter.Spring.EDD + OakDBH + offset(log(Deploy.Duration)), data = deerDataWin)

ModListWinDeer[[14]] <- glm.nb.Wstems_hgt <-
  glm.nb(nSeqs ~ Height_cm + log10(Num_Stems) + offset(log(Deploy.Duration)), data = deerDataWin)

ModListWinDeer[[15]] <- glm.nb.Wstems_oak <-
  glm.nb(nSeqs ~ OakDBH + log10(Num_Stems) + offset(log(Deploy.Duration)), data = deerDataWin)

ModListWinDeer[[16]] <- glm.nb.Whgt_oak <-
  glm.nb(nSeqs ~ Height_cm + OakDBH + offset(log(Deploy.Duration)), data = deerDataWin)

ModListWinDeer[[17]] <- glm.nb.WlogStemsHgt <-
  glm.nb(
    nSeqs ~ Log.in.View + Height_cm + log10(Num_Stems) +
      offset(log(Deploy.Duration)),
    data = deerDataWin
  )

ModListWinDeer[[18]] <- glm.nb.WlogStemsOak <-
  glm.nb(
    nSeqs ~ Log.in.View + OakDBH + log10(Num_Stems) +
      offset(log(Deploy.Duration)),
    data = deerDataWin
  )

ModListWinDeer[[19]] <- glm.nb.WlogStemsEdd <-
  glm.nb(
    nSeqs ~ Log.in.View + Winter.Spring.EDD + log10(Num_Stems) +
      offset(log(Deploy.Duration)),
    data = deerDataWin
  )

ModListWinDeer[[20]] <- glm.nb.WHgtStemsEdd <-
  glm.nb(
    nSeqs ~ Height_cm + Winter.Spring.EDD + log10(Num_Stems) +
      offset(log(Deploy.Duration)),
    data = deerDataWin
  )

ModListWinDeer[[21]] <- glm.nb.WHgtStemsOak <-
  glm.nb(
    nSeqs ~ Height_cm + OakDBH + log10(Num_Stems) +
      offset(log(Deploy.Duration)),
    data = deerDataWin
  )

ModListWinDeer[[22]] <- glm.nb.WHgtEddOak <-
  glm.nb(
    nSeqs ~ Height_cm + OakDBH + Winter.Spring.EDD +
      offset(log(Deploy.Duration)),
    data = deerDataWin
  )
ModListWinDeer[[23]] <- glm.nb.WHgtEddLog <-
  glm.nb(
    nSeqs ~ Height_cm + Log.in.View + Winter.Spring.EDD +
      offset(log(Deploy.Duration)),
    data = deerDataWin
  )
ModListWinDeer[[24]] <- glm.nb.WHgtOakLog <-
  glm.nb(
    nSeqs ~ Height_cm + Log.in.View + OakDBH +
      offset(log(Deploy.Duration)),
    data = deerDataWin
  )

ModListWinDeer[[25]] <- glm.nb.WOakStemsEdd <-
  glm.nb(
    nSeqs ~ log10(Num_Stems) + Winter.Spring.EDD + OakDBH +
      offset(log(Deploy.Duration)),
    data = deerDataWin
  )

ModListWinDeer[[26]] <- glm.nb.WLogOakEdd <-
  glm.nb(
    nSeqs ~ Log.in.View + Winter.Spring.EDD + OakDBH +
      offset(log(Deploy.Duration)),
    data = deerDataWin
  )
modnamesW <- c("Intercept","Log","Height","Oak","Stems","EDD","Log + Stems","Log + EDD","Log + Oak","Log + Height","EDD + Height","EDD + Stems","EDD + Oak","Stems + Height","Stems + Oak","Height + Oak","Log + Stems + Height","Log + Stems + Oak","Log + Stems + EDD","Height + Stems + EDD","Height + Stems + Oak","Height + EDD + Oak","Height + EDD + Log","Height + Oak + Log", "Stems + EDD + Oak", "Log + EDD + Oak")
modtabWinDeer <- aictab(cand.set = ModListWinDeer, modnames = modnamesW)
modtabWinDeer #intercept model is the best

#Summed Model Weights for Deer in Winter from full table. 
WinSW_edd <- sum(modtabWinDeer$AICcWt[grep("EDD", modtabWinDeer$Modnames)])
WinSW_OakDBH <- sum(modtabWinDeer$AICcWt[grep("Oak", modtabWinDeer$Modnames)])
WinSW_Stems <- sum(modtabWinDeer$AICcWt[grep("Stems", modtabWinDeer$Modnames)])
WinSW_Log <- sum(modtabWinDeer$AICcWt[grep("Log", modtabWinDeer$Modnames)])
WinSW_Height <- sum(modtabWinDeer$AICcWt[grep("Height", modtabWinDeer$Modnames)])

#Explained Deviance of best model (Intercept)
1 - glm.nb.WI$deviance / glm.nb.WI$null.deviance #0.0
#Explained Deviance of Full Model
1 - glm.nb.fullW$deviance / glm.nb.fullW$null.deviance #0.1040


# Spring Deer NB Full Model Selection -------------------------------------

glm.nb.fullSP <-
  glm.nb(
    nSeqs ~ Log.in.View + Winter.Spring.EDD + Height_cm + log10(Num_Stems) + OakDBH +
      offset(log(Deploy.Duration)),
    data = deerDataSpr
  )


ModListSprDeer <- list(NA)


ModListSprDeer[[1]] <- glm.nb.SPI <- glm.nb(nSeqs ~ 1 + offset(log(Deploy.Duration)), data = deerDataSpr)

ModListSprDeer[[2]] <- glm.nb.SPlog <-
  glm.nb(nSeqs ~ Log.in.View + offset(log(Deploy.Duration)), data = deerDataSpr)

ModListSprDeer[[3]] <- glm.nb.SPhgt <-
  glm.nb(nSeqs ~ Height_cm + offset(log(Deploy.Duration)), data = deerDataSpr)

ModListSprDeer[[4]] <- glm.nb.SPoak <-
  glm.nb(nSeqs ~ OakDBH + offset(log(Deploy.Duration)), data = deerDataSpr)

ModListSprDeer[[5]] <- glm.nb.SPstems <-
  glm.nb(nSeqs ~ log10(Num_Stems) + offset(log(Deploy.Duration)), data = deerDataSpr)

ModListSprDeer[[6]] <- glm.nb.SPedd <-
  glm.nb(nSeqs ~ Winter.Spring.EDD + offset(log(Deploy.Duration)), data = deerDataSpr)

ModListSprDeer[[7]] <- glm.nb.SPlog_stems <-
  glm.nb(nSeqs ~ Log.in.View + log10(Num_Stems) + offset(log(Deploy.Duration)), data = deerDataSpr)

ModListSprDeer[[8]] <-glm.nb.SPlog_Edd <-
  glm.nb(nSeqs ~ Log.in.View + Winter.Spring.EDD + offset(log(Deploy.Duration)), data = deerDataSpr)

ModListSprDeer[[9]] <- glm.nb.SPlog_Oak <-
  glm.nb(nSeqs ~ Log.in.View + OakDBH + offset(log(Deploy.Duration)), data = deerDataSpr)

ModListSprDeer[[10]] <- glm.nb.SPlog_Hgt <-
  glm.nb(nSeqs ~ Log.in.View + Height_cm + offset(log(Deploy.Duration)), data = deerDataSpr)

ModListSprDeer[[11]] <- glm.nb.SPEdd_Hgt <-
  glm.nb(nSeqs ~ Winter.Spring.EDD + Height_cm + offset(log(Deploy.Duration)), data = deerDataSpr)

ModListSprDeer[[12]] <-  glm.nb.SPEdd_stems <-
  glm.nb(nSeqs ~ Winter.Spring.EDD + log10(Num_Stems) + offset(log(Deploy.Duration)), data = deerDataSpr)

ModListSprDeer[[13]] <- glm.nb.SPEdd_oak <-
  glm.nb(nSeqs ~ Winter.Spring.EDD + OakDBH + offset(log(Deploy.Duration)), data = deerDataSpr)

ModListSprDeer[[14]] <- glm.nb.SPstems_hgt <-
  glm.nb(nSeqs ~ Height_cm + log10(Num_Stems) + offset(log(Deploy.Duration)), data = deerDataSpr)

ModListSprDeer[[15]] <- glm.nb.SPstems_oak <-
  glm.nb(nSeqs ~ OakDBH + log10(Num_Stems) + offset(log(Deploy.Duration)), data = deerDataSpr)

ModListSprDeer[[16]] <- glm.nb.SPhgt_oak <-
  glm.nb(nSeqs ~ Height_cm + OakDBH + offset(log(Deploy.Duration)), data = deerDataSpr)

ModListSprDeer[[17]] <- glm.nb.SPlogStemsHgt <-
  glm.nb(
    nSeqs ~ Log.in.View + Height_cm + log10(Num_Stems) +
      offset(log(Deploy.Duration)),
    data = deerDataSpr
  )

ModListSprDeer[[18]] <- glm.nb.SPlogStemsOak <-
  glm.nb(
    nSeqs ~ Log.in.View + OakDBH + log10(Num_Stems) +
      offset(log(Deploy.Duration)),
    data = deerDataSpr
  )

ModListSprDeer[[19]] <- glm.nb.SPlogStemsEdd <-
  glm.nb(
    nSeqs ~ Log.in.View + Winter.Spring.EDD + log10(Num_Stems) +
      offset(log(Deploy.Duration)),
    data = deerDataSpr
  )

ModListSprDeer[[20]] <- glm.nb.SPHgtStemsEdd <-
  glm.nb(
    nSeqs ~ Height_cm + Winter.Spring.EDD + log10(Num_Stems) +
      offset(log(Deploy.Duration)),
    data = deerDataSpr
  )

ModListSprDeer[[21]] <- glm.nb.SPHgtStemsOak <-
  glm.nb(
    nSeqs ~ Height_cm + OakDBH + log10(Num_Stems) +
      offset(log(Deploy.Duration)),
    data = deerDataSpr
  )

ModListSprDeer[[22]] <- glm.nb.SPHgtEddOak <-
  glm.nb(
    nSeqs ~ Height_cm + OakDBH + Winter.Spring.EDD +
      offset(log(Deploy.Duration)),
    data = deerDataSpr
  )
ModListSprDeer[[23]] <- glm.nb.SPHgtEddLog <-
  glm.nb(
    nSeqs ~ Height_cm + Log.in.View + Winter.Spring.EDD +
      offset(log(Deploy.Duration)),
    data = deerDataSpr
  )
ModListSprDeer[[24]] <- glm.nb.SPHgtOakLog <-
  glm.nb(
    nSeqs ~ Height_cm + Log.in.View + OakDBH +
      offset(log(Deploy.Duration)),
    data = deerDataSpr
  )

ModListSprDeer[[25]] <- glm.nb.SPOakStemsEdd <-
  glm.nb(
    nSeqs ~ log10(Num_Stems) + Winter.Spring.EDD + OakDBH +
      offset(log(Deploy.Duration)),
    data = deerDataSpr
  )

ModListSprDeer[[26]] <- glm.nb.SPLogOakEdd <-
  glm.nb(
    nSeqs ~ Log.in.View + Winter.Spring.EDD + OakDBH +
      offset(log(Deploy.Duration)),
    data = deerDataSpr
  )
modnamesSP <- c("Intercept","Log","Height","Oak","Stems","EDD","Log + Stems","Log + EDD","Log + Oak","Log + Height","EDD + Height","EDD + Stems","EDD + Oak","Stems + Height","Stems + Oak","Height + Oak","Log + Stems + Height","Log + Stems + Oak","Log + Stems + EDD","Height + Stems + EDD","Height + Stems + Oak","Height + EDD + Oak","Height + EDD + Log","Height + Oak + Log", "Stems + EDD + Oak", "Log + EDD + Oak")
modtabSprDeer <- aictab(cand.set = ModListSprDeer, modnames = modnamesSP)
modtabSprDeer #Intercept Only model is best 

#Summed Model Weights for Deer in Spring from full table. 
SprSW_edd <- sum(modtabSprDeer$AICcWt[grep("EDD", modtabSprDeer$Modnames)])
SprSW_OakDBH <- sum(modtabSprDeer$AICcWt[grep("Oak", modtabSprDeer$Modnames)])
SprSW_Stems <- sum(modtabSprDeer$AICcWt[grep("Stems", modtabSprDeer$Modnames)])
SprSW_Log <- sum(modtabSprDeer$AICcWt[grep("Log", modtabSprDeer$Modnames)])
SprSW_Height <- sum(modtabSprDeer$AICcWt[grep("Height", modtabSprDeer$Modnames)])

#Explained Deviance of best model (Intercept)
1 - glm.nb.SPI$deviance / glm.nb.SPI$null.deviance #0.0
#Explained Deviance of Full Model
1 - glm.nb.fullSP$deviance / glm.nb.fullSP$null.deviance #0.1267786

#For Deer, the only thing that comes out important is Logs in Summer, which reduce captures of deer when present. That's it!