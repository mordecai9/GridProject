#Assessement of optimal detection probalitity models for each of Gray Squirrel using occupancy models in Rpresence. Data comes from high resolution 27 camera grid is Posey Hollow. Detection histories were prepared using the script "PrepforDetectionHistory", where R objects of detection hustories were created for each species and each season. These are loaded in here. The script also requires data files prepared for a separate GLM analysis, since these include the site covariates. 


# Clear Environment and Load Packages -------------------------------------

rm(list = ls())
library(tidyverse)
library(RPresence)

# Load Detection Histories Gray Squirrel ------------------------------------------------

load("results/DetHistFallSciurus carolinensis")
DHFallSC <- temp

load("results/DetHistSumSciurus carolinensis")
DHSumSC <- temp

load("results/DetHistWinSciurus carolinensis")
DHWinSC <- temp

load("results/DetHistSpSciurus carolinensis")
DHSpSC <- temp

rm(temp)

# Gray Squirrel Occupancy Models in Spring ---------------------------------

sitenamesFull <- row.names(as.data.frame(DHSpSC)) 
nsitesFull=nrow(DHSpSC) 
nSURVEYsSpSC=ncol(DHSpSC)  #  set number of sites,surveys from det. history data

#create covariates file

load("data/sqDataSpr.RData")
covSpSC <- sqDataSpr %>%
  select(Deployment_Name, Height_cm, Num_Stems, OakDBH, Log.in.View, Squirrel_EDD_WSp)

#Create input file for RPresence
pSpSC <- createPao(DHSpSC,unitcov = covSpSC,title="Gray Squirrel Spring",unitnames=sitenamesFull)

#Test to see if quadratic effect on height improves the height only model

modsHtest <- list(); 
modsHtest[[1]] <- occMod(model=list(psi ~ 1, p ~ Height_cm),data = pSpSC,type = "so")
modsHtest[[2]] <- occMod(model=list(psi ~ 1, p ~ poly(Height_cm, 2, raw = T)),data = pSpSC,type = "so")

resultsHtest <- createAicTable(modsHtest)
resultsHtest$table #there is evidence of a quadratic effect on camera height, so it should be in all models with camera height

#Test to see if there is an interaction between EDD and whether there is a log in view

modsEDDtest <- list(); 
modsEDDtest[[1]] <- occMod(model=list(psi ~ 1, p ~ Log.in.View + Squirrel_EDD_WSp),data = pSpSC,type = "so")
modsEDDtest[[2]] <- occMod(model=list(psi ~ 1, p ~ Log.in.View * Squirrel_EDD_WSp),data = pSpSC,type = "so")

resultsEDDtest <- createAicTable(modsEDDtest)
resultsEDDtest$table #there is evidence of a interaction effect here, so this needs to be in all models below with both of these parameters


mods=list(); i=1
mods[[i]]=occMod(model=list(psi~1, p~1),    data=pSpSC,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~poly(Height_cm, 2, raw = T)),data=pSpSC,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~log10(Num_Stems)),data=pSpSC,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~OakDBH),data=pSpSC,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~Log.in.View),data=pSpSC,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~Squirrel_EDD_WSp),data=pSpSC,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~poly(Height_cm, 2, raw = T) + log10(Num_Stems)),data=pSpSC,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~poly(Height_cm, 2, raw = T) + OakDBH),data=pSpSC,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~poly(Height_cm, 2, raw = T) + Log.in.View),data=pSpSC,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~poly(Height_cm, 2, raw = T) + Squirrel_EDD_WSp),data=pSpSC,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~log10(Num_Stems) + OakDBH),data=pSpSC,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~log10(Num_Stems) + Log.in.View),data=pSpSC,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~log10(Num_Stems) + Squirrel_EDD_WSp),data=pSpSC,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~OakDBH + Log.in.View),data=pSpSC,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~OakDBH + Squirrel_EDD_WSp),data=pSpSC,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~Log.in.View * Squirrel_EDD_WSp),data=pSpSC,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~poly(Height_cm, 2, raw = T) + log10(Num_Stems) + OakDBH),data=pSpSC,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~poly(Height_cm, 2, raw = T) + log10(Num_Stems) + Log.in.View),data=pSpSC,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~poly(Height_cm, 2, raw = T) + log10(Num_Stems) + Squirrel_EDD_WSp),data=pSpSC,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~poly(Height_cm, 2, raw = T) + OakDBH + Log.in.View),data=pSpSC,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~poly(Height_cm, 2, raw = T) + OakDBH + Squirrel_EDD_WSp),data=pSpSC,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~poly(Height_cm, 2, raw = T) + Log.in.View * Squirrel_EDD_WSp),data=pSpSC,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~log10(Num_Stems) + OakDBH + Log.in.View),data=pSpSC,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~log10(Num_Stems) + OakDBH + Squirrel_EDD_WSp),data=pSpSC,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~log10(Num_Stems) + Log.in.View * Squirrel_EDD_WSp),data=pSpSC,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~OakDBH + Log.in.View * Squirrel_EDD_WSp),data=pSpSC,type="so");i=i+1


#     create AIC table of model results and print
resultsSpSC <- createAicTable(mods)
resultsSpSC$table

bestSpSC <- occMod(model=list(psi~1, p~OakDBH + Log.in.View * Squirrel_EDD_WSp),data=pSpSC,type="so") 
bestSpSC$beta$p

#Sum up model weights for each of the 5 covariates
mnames=resultsSpSC$table$Model;
for (s in c("Height","Oak", "Stems", "Log", "EDD")) {
  i <- grep(s,mnames); #  get model names with each covariate name
  wgt=sum(resultsSpSC$table$wgt[i]); #  add up wgts of those models
  cat('CumWgt(',s,')=',wgt,'\n')     #  print sum of weights
}


# Gray Squirrel Occupancy Models in Winter ---------------------------------


sitenamesShort <- row.names(as.data.frame(DHWinSC)) 
nsitesShort=nrow(DHWinSC) 
nSURVEYsWinSC=ncol(DHWinSC)  #  set number of sites,surveys from det. history data

#create covariates file

load("data/sqDataWin.RData") #I'm using general squirrel data here because the site covariates are the same for both species
covWinSC <- sqDataWin %>%
  select(Deployment_Name, Height_cm, Num_Stems, OakDBH, Log.in.View, Squirrel_EDD_WSp)

#Create input file for RPresence
pWinSC <- createPao(DHWinSC,unitcov = covWinSC,title="Gray Squirrel Winter",unitnames=sitenamesShort)

#Test to see if quadratic effect on height improves the height only model

modsHtest <- list(); 
modsHtest[[1]] <- occMod(model=list(psi ~ 1, p ~ Height_cm),data = pWinSC,type = "so")
modsHtest[[2]] <- occMod(model=list(psi ~ 1, p ~ poly(Height_cm, 2, raw = T)),data = pWinSC,type = "so")

resultsHtest <- createAicTable(modsHtest)
resultsHtest$table #there is evidence of a quadratic effect on camera height, so this needs to be in models below

#Test to see if there is an interaction between EDD and whether there is a log in view

modsEDDtest <- list(); 
modsEDDtest[[1]] <- occMod(model=list(psi ~ 1, p ~ Log.in.View + Squirrel_EDD_WSp),data = pWinSC,type = "so")
modsEDDtest[[2]] <- occMod(model=list(psi ~ 1, p ~ Log.in.View * Squirrel_EDD_WSp),data = pWinSC,type = "so")

resultsEDDtest <- createAicTable(modsEDDtest)
resultsEDDtest$table #there is evidence of a interaction effect here, so this needs to be in all models below with both of these parameters


mods=list(); i=1
mods[[i]]=occMod(model=list(psi~1, p~1),    data=pWinSC,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~poly(Height_cm, 2, raw = T)),data=pWinSC,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~log10(Num_Stems)),data=pWinSC,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~OakDBH),data=pWinSC,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~Log.in.View),data=pWinSC,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~Squirrel_EDD_WSp),data=pWinSC,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~poly(Height_cm, 2, raw = T) + log10(Num_Stems)),data=pWinSC,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~poly(Height_cm, 2, raw = T) + OakDBH),data=pWinSC,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~poly(Height_cm, 2, raw = T) + Log.in.View),data=pWinSC,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~poly(Height_cm, 2, raw = T) + Squirrel_EDD_WSp),data=pWinSC,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~log10(Num_Stems) + OakDBH),data=pWinSC,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~log10(Num_Stems) + Log.in.View),data=pWinSC,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~log10(Num_Stems) + Squirrel_EDD_WSp),data=pWinSC,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~OakDBH + Log.in.View),data=pWinSC,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~OakDBH + Squirrel_EDD_WSp),data=pWinSC,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~Log.in.View * Squirrel_EDD_WSp),data=pWinSC,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~poly(Height_cm, 2, raw = T) + log10(Num_Stems) + OakDBH),data=pWinSC,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~poly(Height_cm, 2, raw = T) + log10(Num_Stems) + Log.in.View),data=pWinSC,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~poly(Height_cm, 2, raw = T) + log10(Num_Stems) + Squirrel_EDD_WSp),data=pWinSC,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~poly(Height_cm, 2, raw = T) + OakDBH + Log.in.View),data=pWinSC,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~poly(Height_cm, 2, raw = T) + OakDBH + Squirrel_EDD_WSp),data=pWinSC,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~poly(Height_cm, 2, raw = T) + Log.in.View * Squirrel_EDD_WSp),data=pWinSC,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~log10(Num_Stems) + OakDBH + Log.in.View),data=pWinSC,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~log10(Num_Stems) + OakDBH + Squirrel_EDD_WSp),data=pWinSC,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~log10(Num_Stems) + Log.in.View * Squirrel_EDD_WSp),data=pWinSC,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~OakDBH + Log.in.View * Squirrel_EDD_WSp),data=pWinSC,type="so");i=i+1


#     create AIC table of model results and print
resultsWinSC <- createAicTable(mods)
resultsWinSC$table

bestWinSC <- occMod(model=list(psi~1, p ~ poly(Height_cm, 2, raw = T) + Log.in.View * Squirrel_EDD_WSp),data=pWinSC,type="so")
bestWinSC$beta$p

#Sum up model weights for each of the 5 covariates
mnames=resultsWinSC$table$Model;
for (s in c("Height","Oak", "Stems", "Log", "EDD")) {
  i <- grep(s,mnames); #  get model names with each covariate name
  wgt=sum(resultsWinSC$table$wgt[i]); #  add up wgts of those models
  cat('CumWgt(',s,')=',wgt,'\n')     #  print sum of weights
}

# Gray Squirrel Occupancy Models in Summer ---------------------------------


nSURVEYsSumSC=ncol(DHSumSC)  #  set number of sites,surveys from det. history data

#create covariates file

load("data/sqDataSum.RData") #I'm using general squirrel data here because the site covariates are the same for both species
covSumSC <- sqDataSum %>%
  select(Deployment_Name, Height_cm, Num_Stems, OakDBH, Log.in.View, Squirrel_EDD_4S)

#Create input file for RPresence
pSumSC <- createPao(DHSumSC,unitcov = covSumSC,title="Gray Squirrel Summer",unitnames=sitenamesFull)

#Test to see if quadratic effect on height improves the height only model

modsHtest <- list(); 
modsHtest[[1]] <- occMod(model=list(psi ~ 1, p ~ Height_cm),data = pSumSC,type = "so")
modsHtest[[2]] <- occMod(model=list(psi ~ 1, p ~ poly(Height_cm, 2, raw = T)),data = pSumSC,type = "so")

resultsHtest <- createAicTable(modsHtest)
resultsHtest$table #there is evidence of a quadratic effect on camera height here so need quadratic in all models for height

#Test to see if there is an interaction between EDD and whether there is a log in view

modsEDDtest <- list(); 
modsEDDtest[[1]] <- occMod(model=list(psi ~ 1, p ~ Log.in.View + Squirrel_EDD_4S),data = pSumSC,type = "so")
modsEDDtest[[2]] <- occMod(model=list(psi ~ 1, p ~ Log.in.View * Squirrel_EDD_4S),data = pSumSC,type = "so")

resultsEDDtest <- createAicTable(modsEDDtest)
resultsEDDtest$table #there is no evidence of a interaction effect here

mods=list(); i=1
mods[[i]]=occMod(model=list(psi~1, p~1),    data=pSumSC,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~poly(Height_cm, 2, raw = T)),data=pSumSC,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~log10(Num_Stems)),data=pSumSC,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~OakDBH),data=pSumSC,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~Log.in.View),data=pSumSC,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~Squirrel_EDD_4S),data=pSumSC,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~poly(Height_cm, 2, raw = T) + log10(Num_Stems)),data=pSumSC,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~poly(Height_cm, 2, raw = T) + OakDBH),data=pSumSC,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~poly(Height_cm, 2, raw = T) + Log.in.View),data=pSumSC,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~poly(Height_cm, 2, raw = T) + Squirrel_EDD_4S),data=pSumSC,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~log10(Num_Stems) + OakDBH),data=pSumSC,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~log10(Num_Stems) + Log.in.View),data=pSumSC,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~log10(Num_Stems) + Squirrel_EDD_4S),data=pSumSC,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~OakDBH + Log.in.View),data=pSumSC,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~OakDBH + Squirrel_EDD_4S),data=pSumSC,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~Log.in.View + Squirrel_EDD_4S),data=pSumSC,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~poly(Height_cm, 2, raw = T) + log10(Num_Stems) + OakDBH),data=pSumSC,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~poly(Height_cm, 2, raw = T) + log10(Num_Stems) + Log.in.View),data=pSumSC,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~poly(Height_cm, 2, raw = T) + log10(Num_Stems) + Squirrel_EDD_4S),data=pSumSC,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~poly(Height_cm, 2, raw = T) + OakDBH + Log.in.View),data=pSumSC,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~poly(Height_cm, 2, raw = T) + OakDBH + Squirrel_EDD_4S),data=pSumSC,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~poly(Height_cm, 2, raw = T) + Log.in.View + Squirrel_EDD_4S),data=pSumSC,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~log10(Num_Stems) + OakDBH + Log.in.View),data=pSumSC,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~log10(Num_Stems) + OakDBH + Squirrel_EDD_4S),data=pSumSC,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~log10(Num_Stems) + Log.in.View + Squirrel_EDD_4S),data=pSumSC,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~OakDBH + Log.in.View + Squirrel_EDD_4S),data=pSumSC,type="so");i=i+1


#     create AIC table of model results and print
resultsSumSC <- createAicTable(mods)
resultsSumSC$table

bestSumSC <- occMod(model = list(psi~1, p ~ OakDBH + Log.in.View + Squirrel_EDD_4S),data = pSumSC,type = "so")
bestSumSC$beta$p

#Sum up model weights for each of the 5 covariates
mnames=resultsSumSC$table$Model;
for (s in c("Height","Oak", "Stems", "Log", "EDD")) {
  i <- grep(s,mnames); #  get model names with each covariate name
  wgt=sum(resultsSumSC$table$wgt[i]); #  add up wgts of those models
  cat('CumWgt(',s,')=',wgt,'\n')     #  print sum of weights
}

# Gray Squirrel Occupancy Models in Fall ---------------------------------


nSURVEYsFallSC=ncol(DHFallSC)  #  set number of sites,surveys from det. history data

#create covariates file

load("data/sqDataFall.RData") #I'm using general squirrel data here because the site covariates are the same for both species
covFallSC <- sqDataFall %>%
  select(Deployment_Name, Height_cm, Num_Stems, OakDBH, Log.in.View, Squirrel_EDD_4S)

#Create input file for RPresence
pFallSC <- createPao(DHFallSC,unitcov = covFallSC,title="Gray Squirrel Fall",unitnames=sitenamesShort)

#Test to see if quadratic effect on height improves the height only model

modsHtest <- list(); 
modsHtest[[1]] <- occMod(model=list(psi ~ 1, p ~ Height_cm),data = pFallSC,type = "so")
modsHtest[[2]] <- occMod(model=list(psi ~ 1, p ~ poly(Height_cm, 2, raw = T)),data = pFallSC,type = "so")

resultsHtest <- createAicTable(modsHtest)
resultsHtest$table #there is no evidence of a quadratic effect on camera height here

#Test to see if there is an interaction between EDD and whether there is a log in view

modsEDDtest <- list(); 
modsEDDtest[[1]] <- occMod(model=list(psi ~ 1, p ~ Log.in.View + Squirrel_EDD_4S),data = pFallSC,type = "so")
modsEDDtest[[2]] <- occMod(model=list(psi ~ 1, p ~ Log.in.View * Squirrel_EDD_4S),data = pFallSC,type = "so")

resultsEDDtest <- createAicTable(modsEDDtest)
resultsEDDtest$table #there is evidence of a interaction effect here so it needs to be reflected in models below


mods=list(); i=1
mods[[i]]=occMod(model=list(psi~1, p~1),    data=pFallSC,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~Height_cm),data=pFallSC,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~log10(Num_Stems)),data=pFallSC,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~OakDBH),data=pFallSC,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~Log.in.View),data=pFallSC,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~Squirrel_EDD_4S),data=pFallSC,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~Height_cm + log10(Num_Stems)),data=pFallSC,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~Height_cm + OakDBH),data=pFallSC,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~Height_cm + Log.in.View),data=pFallSC,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~Height_cm + Squirrel_EDD_4S),data=pFallSC,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~log10(Num_Stems) + OakDBH),data=pFallSC,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~log10(Num_Stems) + Log.in.View),data=pFallSC,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~log10(Num_Stems) + Squirrel_EDD_4S),data=pFallSC,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~OakDBH + Log.in.View),data=pFallSC,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~OakDBH + Squirrel_EDD_4S),data=pFallSC,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~Log.in.View * Squirrel_EDD_4S),data=pFallSC,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~Height_cm + log10(Num_Stems) + OakDBH),data=pFallSC,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~Height_cm + log10(Num_Stems) + Log.in.View),data=pFallSC,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~Height_cm + log10(Num_Stems) + Squirrel_EDD_4S),data=pFallSC,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~Height_cm + OakDBH + Log.in.View),data=pFallSC,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~Height_cm + OakDBH + Squirrel_EDD_4S),data=pFallSC,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~Height_cm + Log.in.View * Squirrel_EDD_4S),data=pFallSC,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~log10(Num_Stems) + OakDBH + Log.in.View),data=pFallSC,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~log10(Num_Stems) + OakDBH + Squirrel_EDD_4S),data=pFallSC,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~log10(Num_Stems) + Log.in.View * Squirrel_EDD_4S),data=pFallSC,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~OakDBH + Log.in.View * Squirrel_EDD_4S),data=pFallSC,type="so");i=i+1


#     create AIC table of model results and print
resultsFallSC <- createAicTable(mods)
resultsFallSC$table

bestFallSC <- occMod(model=list(psi~1, p~OakDBH + Log.in.View * Squirrel_EDD_4S),data=pFallSC,type="so") 
bestFallSC$beta$p

#Sum up model weights for each of the 5 covariates
mnames=resultsFallSC$table$Model;
for (s in c("Height","Oak", "Stems", "Log", "EDD")) {
  i <- grep(s,mnames); #  get model names with each covariate name
  wgt=sum(resultsFallSC$table$wgt[i]); #  add up wgts of those models
  cat('CumWgt(',s,')=',wgt,'\n')     #  print sum of weights
}
