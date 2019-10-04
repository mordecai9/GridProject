#Assessement of optimal detection probalitity models for each of Fox Squirrel using occupancy models in Rpresence. Data comes from high resolution 27 camera grid is Posey Hollow. Detection histories were prepared using the script "PrepforDetectionHistory", where R objects of detection hustories were created for each species and each season. These are loaded in here. The script also requires data files prepared for a separate GLM analysis, since these include the site covariates. 


# Clear Environment and Load Packages -------------------------------------

rm(list = ls())
library(tidyverse)
library(RPresence)

# Load Detection Histories Fox Squirrel ------------------------------------------------

load("results/DetHistFallSciurus niger")
DHFallSN <- temp

load("results/DetHistSumSciurus niger")
DHSumSN <- temp

load("results/DetHistWinSciurus niger")
DHWinSN <- temp

load("results/DetHistSpSciurus niger")
DHSpSN <- temp

rm(temp)

# Fox Squirrel Occupancy Models in Spring ---------------------------------

sitenamesFull <- row.names(as.data.frame(DHSpSN)) 
nsitesFull=nrow(DHSpSN) 
nSURVEYsSpSN=ncol(DHSpSN)  #  set number of sites,surveys from det. history data

#create covariates file

load("data/sqDataSpr.RData")
covSpSN <- sqDataSpr %>%
  select(Deployment_Name, Height_cm, Num_Stems, OakDBH, Log.in.View, Squirrel_EDD_WSp)

#Create input file for RPresence
pSpSN <- createPao(DHSpSN,unitcov = covSpSN,title="Fox Squirrel Spring",unitnames=sitenamesFull)

#Test to see if quadratic effect on height improves the height only model

modsHtest <- list(); 
  modsHtest[[1]] <- occMod(model=list(psi ~ 1, p ~ Height_cm),data = pSpSN,type = "so")
  modsHtest[[2]] <- occMod(model=list(psi ~ 1, p ~ poly(Height_cm, 2, raw = T)),data = pSpSN,type = "so")

resultsHtest <- createAicTable(modsHtest)
resultsHtest$table #there is evidence of a quadratic effect on camera height, so it should be in all models with camera height

#Test to see if there is an interaction between EDD and whether there is a log in view

modsEDDtest <- list(); 
modsEDDtest[[1]] <- occMod(model=list(psi ~ 1, p ~ Log.in.View + Squirrel_EDD_WSp),data = pSpSN,type = "so")
modsEDDtest[[2]] <- occMod(model=list(psi ~ 1, p ~ Log.in.View * Squirrel_EDD_WSp),data = pSpSN,type = "so")

resultsEDDtest <- createAicTable(modsEDDtest)
resultsEDDtest$table #there is evidence of a interaction effect here, so this needs to be in all models below with both of these parameters
  
  
mods=list(); i=1
  mods[[i]]=occMod(model=list(psi~1, p~1),    data=pSpSN,type="so");i=i+1
  mods[[i]]=occMod(model=list(psi~1, p~poly(Height_cm, 2, raw = T)),data=pSpSN,type="so");i=i+1
  mods[[i]]=occMod(model=list(psi~1, p~log10(Num_Stems)),data=pSpSN,type="so");i=i+1
  mods[[i]]=occMod(model=list(psi~1, p~OakDBH),data=pSpSN,type="so");i=i+1
  mods[[i]]=occMod(model=list(psi~1, p~Log.in.View),data=pSpSN,type="so");i=i+1
  mods[[i]]=occMod(model=list(psi~1, p~Squirrel_EDD_WSp),data=pSpSN,type="so");i=i+1
  mods[[i]]=occMod(model=list(psi~1, p~poly(Height_cm, 2, raw = T) + log10(Num_Stems)),data=pSpSN,type="so");i=i+1
  mods[[i]]=occMod(model=list(psi~1, p~poly(Height_cm, 2, raw = T) + OakDBH),data=pSpSN,type="so");i=i+1
  mods[[i]]=occMod(model=list(psi~1, p~poly(Height_cm, 2, raw = T) + Log.in.View),data=pSpSN,type="so");i=i+1
  mods[[i]]=occMod(model=list(psi~1, p~poly(Height_cm, 2, raw = T) + Squirrel_EDD_WSp),data=pSpSN,type="so");i=i+1
  mods[[i]]=occMod(model=list(psi~1, p~log10(Num_Stems) + OakDBH),data=pSpSN,type="so");i=i+1
  mods[[i]]=occMod(model=list(psi~1, p~log10(Num_Stems) + Log.in.View),data=pSpSN,type="so");i=i+1
  mods[[i]]=occMod(model=list(psi~1, p~log10(Num_Stems) + Squirrel_EDD_WSp),data=pSpSN,type="so");i=i+1
  mods[[i]]=occMod(model=list(psi~1, p~OakDBH + Log.in.View),data=pSpSN,type="so");i=i+1
  mods[[i]]=occMod(model=list(psi~1, p~OakDBH + Squirrel_EDD_WSp),data=pSpSN,type="so");i=i+1
  mods[[i]]=occMod(model=list(psi~1, p~Log.in.View * Squirrel_EDD_WSp),data=pSpSN,type="so");i=i+1
  mods[[i]]=occMod(model=list(psi~1, p~poly(Height_cm, 2, raw = T) + log10(Num_Stems) + OakDBH),data=pSpSN,type="so");i=i+1
  mods[[i]]=occMod(model=list(psi~1, p~poly(Height_cm, 2, raw = T) + log10(Num_Stems) + Log.in.View),data=pSpSN,type="so");i=i+1
  mods[[i]]=occMod(model=list(psi~1, p~poly(Height_cm, 2, raw = T) + log10(Num_Stems) + Squirrel_EDD_WSp),data=pSpSN,type="so");i=i+1
  mods[[i]]=occMod(model=list(psi~1, p~poly(Height_cm, 2, raw = T) + OakDBH + Log.in.View),data=pSpSN,type="so");i=i+1
  mods[[i]]=occMod(model=list(psi~1, p~poly(Height_cm, 2, raw = T) + OakDBH + Squirrel_EDD_WSp),data=pSpSN,type="so");i=i+1
  mods[[i]]=occMod(model=list(psi~1, p~poly(Height_cm, 2, raw = T) + Log.in.View * Squirrel_EDD_WSp),data=pSpSN,type="so");i=i+1
  mods[[i]]=occMod(model=list(psi~1, p~log10(Num_Stems) + OakDBH + Log.in.View),data=pSpSN,type="so");i=i+1
  mods[[i]]=occMod(model=list(psi~1, p~log10(Num_Stems) + OakDBH + Squirrel_EDD_WSp),data=pSpSN,type="so");i=i+1
  mods[[i]]=occMod(model=list(psi~1, p~log10(Num_Stems) + Log.in.View * Squirrel_EDD_WSp),data=pSpSN,type="so");i=i+1
  mods[[i]]=occMod(model=list(psi~1, p~OakDBH + Log.in.View * Squirrel_EDD_WSp),data=pSpSN,type="so");i=i+1
 
  
#     create AIC table of model results and print
resultsSpSN <- createAicTable(mods)
resultsSpSN$table

bestSpSN <- occMod(model=list(psi~1, p~poly(Height_cm, 2, raw = T) + Squirrel_EDD_WSp),data=pSpSN,type="so") 
unique(bestSpSN$real$p) 
unique(bestSpSN$real$psi)
bestSpSN$beta$p

#Sum up model weights for each of the 5 covariates
mnames=resultsSpSN$table$Model;
for (s in c("Height","Oak", "Stems", "Log", "EDD")) {
  i <- grep(s,mnames); #  get model names with each covariate name
  wgt=sum(resultsSpSN$table$wgt[i]); #  add up wgts of those models
  cat('CumWgt(',s,')=',wgt,'\n')     #  print sum of weights
}

# Response Curve - Height Fox Sq Spring -----------------------------------

seq.SpHgt <- seq(min(sqDataSpr$Height_cm), max(sqDataSpr$Height_cm), length.out = 100)

nd.seq.SpHgt <- data.frame(Height_cm = seq.SpHgt, Squirrel_EDD_WSp = median(sqDataSpr$Squirrel_EDD_WSp) )

pred.SpHgt <- predict(bestSpSN, newdata = nd.seq.SpHgtL, param = "p", conf= 0.95)

predPlotFSpHgt <- data.frame(nd.seq.SpHgt, pred.SpHgt) 

#Not sure why we get some NAs for standard errors...should double check summaries of all models to make sure they ran normally.
load(file = "data/responseTheme.Rdata")

FSqSprHgtPlot <- ggplot(predPlotFSpHgt, aes(x=Height_cm, y=est)) +
  # Confidence region
  geom_ribbon(aes(ymin=lower_0.95, ymax=upper_0.95), alpha=0.25, show.legend = F) +
  # Prediction Lines
  geom_line() +
  xlab("Camera height (cm)") +   
  ylab("Estimated detection probability") +
  geom_text(aes(x = 34, y = .40, label = "A"), size = 8)+
  myTheme

# Fox Squirrel Occupancy Models in Winter ---------------------------------


sitenamesShort <- row.names(as.data.frame(DHWinSN)) 
nsitesShort=nrow(DHWinSN) 
nSURVEYsWinSN=ncol(DHWinSN)  #  set number of sites,surveys from det. history data

#create covariates file

load("data/sqDataWin.RData") #I'm using general squirrel data here because the site covariates are the same for both species
covWinSN <- sqDataWin %>%
  select(Deployment_Name, Height_cm, Num_Stems, OakDBH, Log.in.View, Squirrel_EDD_WSp)

#Create input file for RPresence
pWinSN <- createPao(DHWinSN,unitcov = covWinSN,title="Fox Squirrel Winter",unitnames=sitenamesShort)

#Test to see if quadratic effect on height improves the height only model

modsHtest <- list(); 
modsHtest[[1]] <- occMod(model=list(psi ~ 1, p ~ Height_cm),data = pWinSN,type = "so")
modsHtest[[2]] <- occMod(model=list(psi ~ 1, p ~ poly(Height_cm, 2, raw = T)),data = pWinSN,type = "so")

resultsHtest <- createAicTable(modsHtest)
resultsHtest$table #there is no evidence of a quadratic effect on camera height

#Test to see if there is an interaction between EDD and whether there is a log in view

modsEDDtest <- list(); 
modsEDDtest[[1]] <- occMod(model=list(psi ~ 1, p ~ Log.in.View + Squirrel_EDD_WSp),data = pWinSN,type = "so")
modsEDDtest[[2]] <- occMod(model=list(psi ~ 1, p ~ Log.in.View * Squirrel_EDD_WSp),data = pWinSN,type = "so")

resultsEDDtest <- createAicTable(modsEDDtest)
resultsEDDtest$table #there is evidence of a interaction effect here, so this needs to be in all models below with both of these parameters


mods=list(); i=1
mods[[i]]=occMod(model=list(psi~1, p~1),    data=pWinSN,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~Height_cm),data=pWinSN,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~log10(Num_Stems)),data=pWinSN,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~OakDBH),data=pWinSN,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~Log.in.View),data=pWinSN,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~Squirrel_EDD_WSp),data=pWinSN,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~Height_cm + log10(Num_Stems)),data=pWinSN,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~Height_cm + OakDBH),data=pWinSN,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~Height_cm + Log.in.View),data=pWinSN,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~Height_cm + Squirrel_EDD_WSp),data=pWinSN,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~log10(Num_Stems) + OakDBH),data=pWinSN,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~log10(Num_Stems) + Log.in.View),data=pWinSN,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~log10(Num_Stems) + Squirrel_EDD_WSp),data=pWinSN,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~OakDBH + Log.in.View),data=pWinSN,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~OakDBH + Squirrel_EDD_WSp),data=pWinSN,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~Log.in.View * Squirrel_EDD_WSp),data=pWinSN,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~Height_cm + log10(Num_Stems) + OakDBH),data=pWinSN,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~Height_cm + log10(Num_Stems) + Log.in.View),data=pWinSN,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~Height_cm + log10(Num_Stems) + Squirrel_EDD_WSp),data=pWinSN,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~Height_cm + OakDBH + Log.in.View),data=pWinSN,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~Height_cm + OakDBH + Squirrel_EDD_WSp),data=pWinSN,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~Height_cm + Log.in.View * Squirrel_EDD_WSp),data=pWinSN,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~log10(Num_Stems) + OakDBH + Log.in.View),data=pWinSN,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~log10(Num_Stems) + OakDBH + Squirrel_EDD_WSp),data=pWinSN,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~log10(Num_Stems) + Log.in.View * Squirrel_EDD_WSp),data=pWinSN,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~OakDBH + Log.in.View * Squirrel_EDD_WSp),data=pWinSN,type="so");i=i+1


#     create AIC table of model results and print
resultsWinSN <- createAicTable(mods)
resultsWinSN$table

bestWinSN <- occMod(model=list(psi~1, p ~ Log.in.View * Squirrel_EDD_WSp),data=pWinSN,type="so")

#Sum up model weights for each of the 5 covariates
mnames=resultsWinSN$table$Model;
for (s in c("Height","Oak", "Stems", "Log", "EDD")) {
  i <- grep(s,mnames); #  get model names with each covariate name
  wgt=sum(resultsWinSN$table$wgt[i]); #  add up wgts of those models
  cat('CumWgt(',s,')=',wgt,'\n')     #  print sum of weights
}

# Fox Squirrel Occupancy Models in Summer ---------------------------------


nSURVEYsSumSN=ncol(DHSumSN)  #  set number of sites,surveys from det. history data

#create covariates file

load("data/sqDataSum.RData") #I'm using general squirrel data here because the site covariates are the same for both species
covSumSN <- sqDataSum %>%
  select(Deployment_Name, Height_cm, Num_Stems, OakDBH, Log.in.View, Squirrel_EDD_4S)

#Create input file for RPresence
pSumSN <- createPao(DHSumSN,unitcov = covSumSN,title="Fox Squirrel Summer",unitnames=sitenamesFull)

#Test to see if quadratic effect on height improves the height only model

modsHtest <- list(); 
modsHtest[[1]] <- occMod(model=list(psi ~ 1, p ~ Height_cm),data = pSumSN,type = "so")
modsHtest[[2]] <- occMod(model=list(psi ~ 1, p ~ poly(Height_cm, 2, raw = T)),data = pSumSN,type = "so")
modsHtest[[2]]$beta$p


resultsHtest <- createAicTable(modsHtest)
resultsHtest$table #there is evidence of a quadratic effect on camera height here so need quadratic in all models for height

#Test to see if there is an interaction between EDD and whether there is a log in view

modsEDDtest <- list(); 
modsEDDtest[[1]] <- occMod(model=list(psi ~ 1, p ~ Log.in.View + Squirrel_EDD_4S),data = pSumSN,type = "so")
modsEDDtest[[2]] <- occMod(model=list(psi ~ 1, p ~ Log.in.View * Squirrel_EDD_4S),data = pSumSN,type = "so")

resultsEDDtest <- createAicTable(modsEDDtest)
resultsEDDtest$table #there is no evidence of a interaction effect here

mods=list(); i=1
mods[[i]]=occMod(model=list(psi~1, p~1),    data=pSumSN,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~poly(Height_cm, 2, raw = T)),data=pSumSN,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~log10(Num_Stems)),data=pSumSN,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~OakDBH),data=pSumSN,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~Log.in.View),data=pSumSN,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~Squirrel_EDD_4S),data=pSumSN,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~poly(Height_cm, 2, raw = T) + log10(Num_Stems)),data=pSumSN,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~poly(Height_cm, 2, raw = T) + OakDBH),data=pSumSN,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~poly(Height_cm, 2, raw = T) + Log.in.View),data=pSumSN,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~poly(Height_cm, 2, raw = T) + Squirrel_EDD_4S),data=pSumSN,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~log10(Num_Stems) + OakDBH),data=pSumSN,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~log10(Num_Stems) + Log.in.View),data=pSumSN,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~log10(Num_Stems) + Squirrel_EDD_4S),data=pSumSN,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~OakDBH + Log.in.View),data=pSumSN,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~OakDBH + Squirrel_EDD_4S),data=pSumSN,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~Log.in.View + Squirrel_EDD_4S),data=pSumSN,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~poly(Height_cm, 2, raw = T) + log10(Num_Stems) + OakDBH),data=pSumSN,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~poly(Height_cm, 2, raw = T) + log10(Num_Stems) + Log.in.View),data=pSumSN,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~poly(Height_cm, 2, raw = T) + log10(Num_Stems) + Squirrel_EDD_4S),data=pSumSN,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~poly(Height_cm, 2, raw = T) + OakDBH + Log.in.View),data=pSumSN,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~poly(Height_cm, 2, raw = T) + OakDBH + Squirrel_EDD_4S),data=pSumSN,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~poly(Height_cm, 2, raw = T) + Log.in.View + Squirrel_EDD_4S),data=pSumSN,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~log10(Num_Stems) + OakDBH + Log.in.View),data=pSumSN,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~log10(Num_Stems) + OakDBH + Squirrel_EDD_4S),data=pSumSN,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~log10(Num_Stems) + Log.in.View + Squirrel_EDD_4S),data=pSumSN,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~OakDBH + Log.in.View + Squirrel_EDD_4S),data=pSumSN,type="so");i=i+1


#     create AIC table of model results and print
resultsSumSN <- createAicTable(mods)
resultsSumSN$table

bestSumSN <- occMod(model=list(psi~1, p~poly(Height_cm, 2, raw = T) + Log.in.View + Squirrel_EDD_4S),data=pSumSN,type="so")
summary(bestSumSN)
unique(bestSumSN$real$p) #the NAs indicate maybe this model didn't work...
unique(bestSumSN$real$psi)
bestSumSN$beta$p

#Testing model with only camera height quadratic term, since it looks like all height models had trouble
FsSumHgt <- occMod(model=list(psi~1, p~poly(Height_cm, 2, raw = T)),data=pSumSN,type="so")
FsSumHgt$beta$p
FsSumHgt$aic

#Trying without raw
FsSumHgt2 <- occMod(model=list(psi~1, p~poly(Height_cm, 2)),data=pSumSN,type="so")
FsSumHgt2$beta$psi
FsSumHgt2$aic

#Trying manually. Same results as raw with poly.
FsSumHgt3 <- occMod(model=list(psi~1, p~Height_cm + I(Height_cm^2)),data=pSumSN,type="so")
FsSumHgt3$beta$psi
FsSumHgt3$beta$p
FsSumHgt3$aic

#Sum up model weights for each of the 5 covariates
mnames=resultsSumSN$table$Model;
for (s in c("Height","Oak", "Stems", "Log", "EDD")) {
  i <- grep(s,mnames); #  get model names with each covariate name
  wgt=sum(resultsSumSN$table$wgt[i]); #  add up wgts of those models
  cat('CumWgt(',s,')=',wgt,'\n')     #  print sum of weights
}


# Response Curve - Height Fox Sq Summer -----------------------------------

seq.SHgt <- seq(min(sqDataSum$Height_cm), max(sqDataSum$Height_cm), length.out = 100)

nd.seq.SHgtL <- data.frame(Height_cm = seq.SHgt, Log.in.View = factor("YES", levels=c("NO", "YES")), Squirrel_EDD_4S = median(sqDataSum$Squirrel_EDD_4S) )

pred.SHgtL <- predict(bestSumSN, newdata = nd.seq.SHgtL, param = "p", conf= 0.95)

predPlotFSumHgt <- data.frame(nd.seq.SHgtL, pred.SHgtL) 

#Not sure why we get some NAs for standard errors...should double check summaries of all models to make sure they ran normally. Looks like nearly all models with height quadratic term had convergence problems.
load(file = "data/responseTheme.Rdata")

FSqSumHgtPlot <- ggplot(predPlotFSumHgt, aes(x=Height_cm, y=est)) +
  # Confidence region
  geom_ribbon(aes(ymin=lower_0.95, ymax=upper_0.95), alpha=0.25, show.legend = F) +
  # Prediction Lines
  geom_line() +
  xlab("Camera Height (cm)") +   
  ylab("Estimated detection probability") +
  geom_text(aes(x = 34, y = .22, label = "B"), size = 8)+
  myTheme

# Fox Squirrel Occupancy Models in Fall ---------------------------------


nSURVEYsFallSN=ncol(DHFallSN)  #  set number of sites,surveys from det. history data

#create covariates file

load("data/sqDataFall.RData") #I'm using general squirrel data here because the site covariates are the same for both species
covFallSN <- sqDataFall %>%
  select(Deployment_Name, Height_cm, Num_Stems, OakDBH, Log.in.View, Squirrel_EDD_4S)

#Create input file for RPresence
pFallSN <- createPao(DHFallSN,unitcov = covFallSN,title="Fox Squirrel Fall",unitnames=sitenamesShort)

#Test to see if quadratic effect on height improves the height only model

modsHtest <- list(); 
modsHtest[[1]] <- occMod(model=list(psi ~ 1, p ~ Height_cm),data = pFallSN,type = "so")
modsHtest[[2]] <- occMod(model=list(psi ~ 1, p ~ poly(Height_cm, 2, raw = T)),data = pFallSN,type = "so")

resultsHtest <- createAicTable(modsHtest)
resultsHtest$table #there is no evidence of a quadratic effect on camera height here

#Test to see if there is an interaction between EDD and whether there is a log in view

modsEDDtest <- list(); 
modsEDDtest[[1]] <- occMod(model=list(psi ~ 1, p ~ Log.in.View + Squirrel_EDD_4S),data = pFallSN,type = "so")
modsEDDtest[[2]] <- occMod(model=list(psi ~ 1, p ~ Log.in.View * Squirrel_EDD_4S),data = pFallSN,type = "so")

resultsEDDtest <- createAicTable(modsEDDtest)
resultsEDDtest$table #there is no evidence of a interaction effect here


mods=list(); i=1
mods[[i]]=occMod(model=list(psi~1, p~1),    data=pFallSN,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~Height_cm),data=pFallSN,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~log10(Num_Stems)),data=pFallSN,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~OakDBH),data=pFallSN,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~Log.in.View),data=pFallSN,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~Squirrel_EDD_4S),data=pFallSN,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~Height_cm + log10(Num_Stems)),data=pFallSN,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~Height_cm + OakDBH),data=pFallSN,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~Height_cm + Log.in.View),data=pFallSN,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~Height_cm + Squirrel_EDD_4S),data=pFallSN,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~log10(Num_Stems) + OakDBH),data=pFallSN,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~log10(Num_Stems) + Log.in.View),data=pFallSN,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~log10(Num_Stems) + Squirrel_EDD_4S),data=pFallSN,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~OakDBH + Log.in.View),data=pFallSN,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~OakDBH + Squirrel_EDD_4S),data=pFallSN,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~Log.in.View + Squirrel_EDD_4S),data=pFallSN,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~Height_cm + log10(Num_Stems) + OakDBH),data=pFallSN,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~Height_cm + log10(Num_Stems) + Log.in.View),data=pFallSN,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~Height_cm + log10(Num_Stems) + Squirrel_EDD_4S),data=pFallSN,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~Height_cm + OakDBH + Log.in.View),data=pFallSN,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~Height_cm + OakDBH + Squirrel_EDD_4S),data=pFallSN,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~Height_cm + Log.in.View + Squirrel_EDD_4S),data=pFallSN,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~log10(Num_Stems) + OakDBH + Log.in.View),data=pFallSN,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~log10(Num_Stems) + OakDBH + Squirrel_EDD_4S),data=pFallSN,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~log10(Num_Stems) + Log.in.View + Squirrel_EDD_4S),data=pFallSN,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~OakDBH + Log.in.View + Squirrel_EDD_4S),data=pFallSN,type="so");i=i+1


#     create AIC table of model results and print
resultsFallSN <- createAicTable(mods)
resultsFallSN$table

bestFallSN <- occMod(model=list(psi~1, p~log10(Num_Stems) + Log.in.View + Squirrel_EDD_4S),data=pFallSN,type="so") 

#Sum up model weights for each of the 5 covariates
mnames=resultsFallSN$table$Model;
for (s in c("Height","Oak", "Stems", "Log", "EDD")) {
  i <- grep(s,mnames); #  get model names with each covariate name
  wgt=sum(resultsFallSN$table$wgt[i]); #  add up wgts of those models
  cat('CumWgt(',s,')=',wgt,'\n')     #  print sum of weights
}
