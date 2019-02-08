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


# Negative Binomial Regression - Black Bear --------------------------------------
#Summer

glm.nb.fullS <-
  glm.nb(
    nSeqs ~ Log.in.View + Summer.Fall.EDD + Height_cm + log10(Num_Stems) + OakDBH +
      offset(log(Deploy.Duration)),
    data = bearDataSum
  )
summary(glm.nb.fullS)

glm.nb.SI <- glm.nb(nSeqs ~ 1 + offset(log(Deploy.Duration)), data = bearDataSum)

glm.nb.Slog <-
  glm.nb(nSeqs ~ Log.in.View + offset(log(Deploy.Duration)), data = bearDataSum)

glm.nb.Shgt <-
  glm.nb(nSeqs ~ Height_cm + offset(log(Deploy.Duration)), data = bearDataSum)

glm.nb.Soak <-
  glm.nb(nSeqs ~ OakDBH + offset(log(Deploy.Duration)), data = bearDataSum)

glm.nb.Sstems <-
  glm.nb(nSeqs ~ Num_Stems + offset(log(Deploy.Duration)), data = bearDataSum)

glm.nb.Sedd <-
  glm.nb(nSeqs ~ Summer.Fall.EDD + offset(log(Deploy.Duration)), data = bearDataSum)

glm.nb.Slog_edd <-
  glm.nb(nSeqs ~ Summer.Fall.EDD + Log.in.View + offset(log(Deploy.Duration)), data = bearDataSum)

mod.namesBS <- c("Full", "log", "height", "oaks", "stems", "EDD", "log + EDD", "Intercept")
glmNBBS <- cbind(mod.namesBS, c(AICc(glm.nb.fullS), AICc(glm.nb.Slog), AICc(glm.nb.Shgt), AICc(glm.nb.Soak), AICc(glm.nb.Sstems), AICc(glm.nb.Sedd), AICc(glm.nb.Slog_edd), AICc(glm.nb.SI)))
glmNBBS

summary(glm.nb.Slog) #This is likely the best model. Indicates bear are captured less frequenctly when there is a log in view, in summer.

#Fall
glm.nb.fullF <- glm.nb(nSeqs ~ Log.in.View + Summer.Fall.EDD + Height_cm + log10(Num_Stems) + OakDBH +offset(log(Deploy.Duration)), data = bearDataFall)
summary(glm.nb.fullF)

glm.nb.FI <- glm.nb(nSeqs ~ 1 + offset(log(Deploy.Duration)), data = bearDataFall)

glm.nb.Flog <-
  glm.nb(nSeqs ~ Log.in.View + offset(log(Deploy.Duration)), data = bearDataFall)

glm.nb.Fhgt <-
  glm.nb(nSeqs ~ Height_cm + offset(log(Deploy.Duration)), data = bearDataFall)

glm.nb.FhgtLog <-
  glm.nb(nSeqs ~ Height_cm + Log.in.View + offset(log(Deploy.Duration)), data = bearDataFall)

glm.nb.Foak <-
  glm.nb(nSeqs ~ OakDBH + offset(log(Deploy.Duration)), data = bearDataFall)

glm.nb.Fstems <-
  glm.nb(nSeqs ~ Num_Stems + offset(log(Deploy.Duration)), data = bearDataFall)

glm.nb.Fedd <-
  glm.nb(nSeqs ~ Summer.Fall.EDD + offset(log(Deploy.Duration)), data = bearDataFall)

mod.namesBF <- c("Full", "log", "height", "log + height", "oaks", "stems", "EDD", "Intercept")
glmNBBF <- cbind(mod.namesBF, c(AICc(glm.nb.fullF), AICc(glm.nb.Flog), AICc(glm.nb.Fhgt), AICc(glm.nb.FhgtLog), AICc(glm.nb.Foak), AICc(glm.nb.Fstems), AICc(glm.nb.Fedd), AICc(glm.nb.FI)))
glmNBBF

summary(glm.nb.Foak) # Seems like best model, indicates bear more likely to be photographed in front of cameras with more/bigger oaks in front of it, in the Fall


#Spring
glm.nb.fullSp <- glm.nb(nSeqs ~ Log.in.View + Winter.Spring.EDD + Height_cm + log10(Num_Stems) + OakDBH +offset(log(Deploy.Duration)), data = bearDataSpr)
summary(glm.nb.fullSp)

glm.nb.SpI <- glm.nb(nSeqs ~ 1 + offset(log(Deploy.Duration)), data = bearDataSpr)

glm.nb.Splog <-
  glm.nb(nSeqs ~ Log.in.View + offset(log(Deploy.Duration)), data = bearDataSpr)

glm.nb.Sphgt <-
  glm.nb(nSeqs ~ Height_cm + offset(log(Deploy.Duration)), data = bearDataSpr)

glm.nb.SphgtLog <-
  glm.nb(nSeqs ~ Height_cm + Log.in.View + offset(log(Deploy.Duration)), data = bearDataSpr)

glm.nb.Spoak <-
  glm.nb(nSeqs ~ OakDBH + offset(log(Deploy.Duration)), data = bearDataSpr)

glm.nb.Spstems <-
  glm.nb(nSeqs ~ Num_Stems + offset(log(Deploy.Duration)), data = bearDataSpr)

glm.nb.Spedd <-
  glm.nb(nSeqs ~ Winter.Spring.EDD + offset(log(Deploy.Duration)), data = bearDataSpr)

mod.namesBSp <- c("Full", "log", "height", "log + height", "oaks", "stems", "EDD", "Intercept")
glmNBBSp <- cbind(mod.namesBSp, c(AICc(glm.nb.fullSp), AICc(glm.nb.Splog), AICc(glm.nb.Sphgt), AICc(glm.nb.SphgtLog), AICc(glm.nb.Spoak), AICc(glm.nb.Spstems), AICc(glm.nb.Spedd), AICc(glm.nb.SpI)))
glmNBBSp #No useful parameters for Spring and Black Bear. I think maybe we didn't have many sequences overall

#____________________________________
#Response Curves for Covariates
#____________________________________


#Bear in Fall, influence of dbh of oaks in front of camera on detection rates
summary(glm.nb.Foak)


#Response Curves for OakDBH. No other parameters were important.
summary(bearDataFall$OakDBH)
seq.oak <- seq(min(bearDataFall$OakDBH), max(bearDataFall$OakDBH), length.out = 100)
nd.seq.oak <- data.frame(Deploy.Duration = 61, OakDBH = seq.oak)
pred.oak <- predict(glm.nb.Foak, newdata = nd.seq.oak, se=TRUE, type = "link")
plot(nSeqs ~ OakDBH, data=bearDataFall, pch=16, col=rgb(.75,.25,0,.5), las=1, xlab = "Total DBH Oak", ylab = "# Bear Sequences in Fall")
lines(nd.seq.oak$OakDBH, exp(pred.oak$fit), col="dark olive green")
lines(nd.seq.oak$OakDBH, exp(pred.oak$fit + 1.96*pred.oak$se.fit), lty=2, col="dark olive green")
lines(nd.seq.oak$OakDBH, exp(pred.oak$fit - 1.96*pred.oak$se.fit), lty=2, col="dark olive green")
