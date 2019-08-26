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

#Sum up model weights for each of the 5 covariates
mnames=resultsSpSN$table$Model;
for (s in c("Height","Oak", "Stems", "Log", "EDD")) {
  i <- grep(s,mnames); #  get model names with each covariate name
  wgt=sum(resultsSpSN$table$wgt[i]); #  add up wgts of those models
  cat('CumWgt(',s,')=',wgt,'\n')     #  print sum of weights
}


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

#Sum up model weights for each of the 5 covariates
mnames=resultsSumSN$table$Model;
for (s in c("Height","Oak", "Stems", "Log", "EDD")) {
  i <- grep(s,mnames); #  get model names with each covariate name
  wgt=sum(resultsSumSN$table$wgt[i]); #  add up wgts of those models
  cat('CumWgt(',s,')=',wgt,'\n')     #  print sum of weights
}

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
