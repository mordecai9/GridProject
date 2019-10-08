#Assessment of optimal detection probalitity models for each of Deerusing occupancy models in Rpresence. Data comes from high resolution 27 camera grid is Posey Hollow. Detection histories were prepared using the OCript "PrepforDetectionHistory", where R objects of detection hustories were created for each species and each season. These are loaded in here. The OCript also requires data files prepared for a separate GLM analysis, since these include the site covariates. 


# Clear Environment and Load Packages -------------------------------------

rm(list = ls())
library(tidyverse)
library(RPresence)

# Load Detection Histories Deer------------------------------------------------

load("results/DetHistFallOdocoileus virginianus")
DHFallOV <- temp

load("results/DetHistSumOdocoileus virginianus")
DHSumOV <- temp

load("results/DetHistWinOdocoileus virginianus")
DHWinOV <- temp

load("results/DetHistSpOdocoileus virginianus")
DHSpOV <- temp

rm(temp)

# Deer Occupancy Models in Spring ---------------------------------

sitenamesFull <- row.names(as.data.frame(DHSpOV)) 
nsitesFull=nrow(DHSpOV) 
nSURVEYsSpOV=ncol(DHSpOV)  #  set number of sites,surveys from det. history data

#create covariates file

load("data/deerDataSpr.RData")
covSpOV <- deerDataSpr %>%
  select(Deployment_Name, Height_cm, Num_Stems, OakDBH, Log.in.View, Winter.Spring.EDD)

#Created scaled version of the covariates file
covSpOV2 <- covSpOV %>%
  mutate(Height_cm = scale(Height_cm),
         Num_Stems = scale(Num_Stems),
         OakDBH = scale(OakDBH),
         Winter.Spring.EDD = scale(Winter.Spring.EDD))

#Create input file for RPresence
pSpOV <- createPao(DHSpOV,unitcov = covSpOV,title="DeerSpring",unitnames=sitenamesFull)
pSpOV2 <- createPao(DHSpOV,unitcov = covSpOV2,title="DeerSpring-scaled",unitnames=sitenamesFull)


#Test to see if there is an interaction between EDD and whether there is a log in view - did not have this as a rule in protocol for paper.

modsEDDtest <- list(); 
modsEDDtest[[1]] <- occMod(model=list(psi ~ 1, p ~ Log.in.View + Winter.Spring.EDD),data = pSpOV2,type = "so")
modsEDDtest[[2]] <- occMod(model=list(psi ~ 1, p ~ Log.in.View * Winter.Spring.EDD),data = pSpOV2,type = "so")

resultsEDDtest <- createAicTable(modsEDDtest)
resultsEDDtest$table #there is no evidence of a interaction effect here


mods=list(); i=1
mods[[i]]=occMod(model=list(psi~1, p~1),    data=pSpOV2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~Height_cm),data=pSpOV2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~Num_Stems),data=pSpOV2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~OakDBH),data=pSpOV2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~Log.in.View),data=pSpOV2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~Winter.Spring.EDD),data=pSpOV2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~Height_cm + Num_Stems),data=pSpOV2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~Height_cm + OakDBH),data=pSpOV2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~Height_cm + Log.in.View),data=pSpOV2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~Height_cm + Winter.Spring.EDD),data=pSpOV2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~Num_Stems + OakDBH),data=pSpOV2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~Num_Stems + Log.in.View),data=pSpOV2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~Num_Stems + Winter.Spring.EDD),data=pSpOV2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~OakDBH + Log.in.View),data=pSpOV2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~OakDBH + Winter.Spring.EDD),data=pSpOV2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~Log.in.View + Winter.Spring.EDD),data=pSpOV2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~Height_cm + Num_Stems + OakDBH),data=pSpOV2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~Height_cm + Num_Stems + Log.in.View),data=pSpOV2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~Height_cm + Num_Stems + Winter.Spring.EDD),data=pSpOV2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~Height_cm + OakDBH + Log.in.View),data=pSpOV2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~Height_cm + OakDBH + Winter.Spring.EDD),data=pSpOV2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~Height_cm + Log.in.View + Winter.Spring.EDD),data=pSpOV2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~Num_Stems + OakDBH + Log.in.View),data=pSpOV2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~Num_Stems + OakDBH + Winter.Spring.EDD),data=pSpOV2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~Num_Stems + Log.in.View + Winter.Spring.EDD),data=pSpOV2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~OakDBH + Log.in.View + Winter.Spring.EDD),data=pSpOV2,type="so");i=i+1


#     create AIC table of model results and print
resultsSpOV <- createAicTable(mods)
resultsSpOV$table

bestSpOV <- occMod(model=list(psi~1, p~Num_Stems + Log.in.View),data=pSpOV2,type="so") 
bestSpOV$beta$p

#Sum up model weights for each of the 5 covariates
mnames=resultsSpOV$table$Model;
for (s in c("Height","Oak", "Stems", "Log", "EDD")) {
  i <- grep(s,mnames); #  get model names with each covariate name
  wgt=sum(resultsSpOV$table$wgt[i]); #  add up wgts of those models
  cat('CumWgt(',s,')=',wgt,'\n')     #  print sum of weights
}


# Deer Occupancy Models in Winter ---------------------------------


sitenamesShort <- row.names(as.data.frame(DHWinOV)) 
nsitesShort=nrow(DHWinOV) 
nSURVEYsWinOV=ncol(DHWinOV)  #  set number of sites,surveys from det. history data

#create covariates file

load("data/deerDataWin.RData") #I'm using general squirrel data here because the site covariates are the same for both species
covWinOV <- deerDataWin %>%
  select(Deployment_Name, Height_cm, Num_Stems, OakDBH, Log.in.View, Winter.Spring.EDD)

covWinOV2 <- covWinOV %>%
  mutate(Height_cm = scale(Height_cm),
         Num_Stems = scale(Num_Stems),
         OakDBH = scale(OakDBH),
         Winter.Spring.EDD = scale(Winter.Spring.EDD))

#Create input file for RPresence
pWinOV <- createPao(DHWinOV,unitcov = covWinOV,title="DeerWinter",unitnames=sitenamesShort)
pWinOV2 <- createPao(DHWinOV,unitcov = covWinOV2,title="DeerWinter-scaled",unitnames=sitenamesShort)


#Test to see if there is an interaction between EDD and whether there is a log in view

#modsEDDtest <- list(); 
#modsEDDtest[[1]] <- occMod(model=list(psi ~ 1, p ~ Log.in.View + Winter.Spring.EDD),data = pWinOV,type = "so")
#modsEDDtest[[2]] <- occMod(model=list(psi ~ 1, p ~ Log.in.View * Winter.Spring.EDD),data = pWinOV,type = "so")

#resultsEDDtest <- createAicTable(modsEDDtest)
#resultsEDDtest$table #there is evidence of an interaction effect here, so this needs to be in all models below with both of these parameters


mods=list(); i=1
mods[[i]]=occMod(model=list(psi~1, p~1),    data=pWinOV2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~Height_cm),data=pWinOV2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~Num_Stems),data=pWinOV2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~OakDBH),data=pWinOV2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~Log.in.View),data=pWinOV2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~Winter.Spring.EDD),data=pWinOV2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~Height_cm + Num_Stems),data=pWinOV2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~Height_cm + OakDBH),data=pWinOV2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~Height_cm + Log.in.View),data=pWinOV2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~Height_cm + Winter.Spring.EDD),data=pWinOV2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~Num_Stems + OakDBH),data=pWinOV2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~Num_Stems + Log.in.View),data=pWinOV2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~Num_Stems + Winter.Spring.EDD),data=pWinOV2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~OakDBH + Log.in.View),data=pWinOV2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~OakDBH + Winter.Spring.EDD),data=pWinOV2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~Log.in.View + Winter.Spring.EDD),data=pWinOV2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~Height_cm + Num_Stems + OakDBH),data=pWinOV2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~Height_cm + Num_Stems + Log.in.View),data=pWinOV2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~Height_cm + Num_Stems + Winter.Spring.EDD),data=pWinOV2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~Height_cm + OakDBH + Log.in.View),data=pWinOV2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~Height_cm + OakDBH + Winter.Spring.EDD),data=pWinOV2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~Height_cm + Log.in.View + Winter.Spring.EDD),data=pWinOV2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~Num_Stems + OakDBH + Log.in.View),data=pWinOV2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~Num_Stems + OakDBH + Winter.Spring.EDD),data=pWinOV2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~Num_Stems + Log.in.View + Winter.Spring.EDD),data=pWinOV2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~OakDBH + Log.in.View + Winter.Spring.EDD),data=pWinOV2,type="so");i=i+1


#     create AIC table of model results and print
resultsWinOV <- createAicTable(mods)
resultsWinOV$table

bestWinOV <- occMod(model=list(psi~1, p ~ Num_Stems + Winter.Spring.EDD),data=pWinOV2,type="so")
bestWinOV$beta$p

#Sum up model weights for each of the 5 covariates
mnames=resultsWinOV$table$Model;
for (s in c("Height","Oak", "Stems", "Log", "EDD")) {
  i <- grep(s,mnames); #  get model names with each covariate name
  wgt=sum(resultsWinOV$table$wgt[i]); #  add up wgts of those models
  cat('CumWgt(',s,')=',wgt,'\n')     #  print sum of weights
}

# Deer Occupancy Models in Summer ---------------------------------


nSURVEYsSumOV=ncol(DHSumOV)  #  set number of sites,surveys from det. history data

#create covariates file

load("data/deerDataSum.RData") #I'm using general squirrel data here because the site covariates are the same for both species
covSumOV <- deerDataSum %>%
  select(Deployment_Name, Height_cm, Num_Stems, OakDBH, Log.in.View, Summer.Fall.EDD)

#Created scaled version of the covariates file
covSumOV2 <- covSumOV %>%
  mutate(Height_cm = scale(Height_cm),
         Num_Stems = scale(Num_Stems),
         OakDBH = scale(OakDBH),
         Summer.Fall.EDD = scale(Summer.Fall.EDD))

#Create input file for RPresence
pSumOV <- createPao(DHSumOV,unitcov = covSumOV,title="DeerSummer",unitnames=sitenamesFull)
pSumOV2 <- createPao(DHSumOV,unitcov = covSumOV2,title="DeerSummer-scaled",unitnames=sitenamesFull)



mods=list(); i=1
mods[[i]]=occMod(model=list(psi~1, p~1),    data=pSumOV2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~Height_cm),data=pSumOV2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~Num_Stems),data=pSumOV2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~OakDBH),data=pSumOV2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~Log.in.View),data=pSumOV2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~Summer.Fall.EDD),data=pSumOV2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~Height_cm + Num_Stems),data=pSumOV2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~Height_cm + OakDBH),data=pSumOV2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~Height_cm + Log.in.View),data=pSumOV2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~Height_cm + Summer.Fall.EDD),data=pSumOV2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~Num_Stems + OakDBH),data=pSumOV2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~Num_Stems + Log.in.View),data=pSumOV2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~Num_Stems + Summer.Fall.EDD),data=pSumOV2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~OakDBH + Log.in.View),data=pSumOV2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~OakDBH + Summer.Fall.EDD),data=pSumOV2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~Log.in.View + Summer.Fall.EDD),data=pSumOV2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~Height_cm + Num_Stems + OakDBH),data=pSumOV2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~Height_cm + Num_Stems + Log.in.View),data=pSumOV2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~Height_cm + Num_Stems + Summer.Fall.EDD),data=pSumOV2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~Height_cm + OakDBH + Log.in.View),data=pSumOV2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~Height_cm + OakDBH + Summer.Fall.EDD),data=pSumOV2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~Height_cm + Log.in.View + Summer.Fall.EDD),data=pSumOV2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~Num_Stems + OakDBH + Log.in.View),data=pSumOV2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~Num_Stems + OakDBH + Summer.Fall.EDD),data=pSumOV2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~Num_Stems + Log.in.View + Summer.Fall.EDD),data=pSumOV2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~OakDBH + Log.in.View + Summer.Fall.EDD),data=pSumOV2,type="so");i=i+1


#     create AIC table of model results and print
resultsSumOV <- createAicTable(mods)
resultsSumOV$table

bestSumOV <- occMod(model = list(psi~1, p ~ Log.in.View + Summer.Fall.EDD),data = pSumOV2,type = "so")
bestSumOV$beta$p

#Sum up model weights for each of the 5 covariates
mnames=resultsSumOV$table$Model;
for (s in c("Height","Oak", "Stems", "Log", "EDD")) {
  i <- grep(s,mnames); #  get model names with each covariate name
  wgt=sum(resultsSumOV$table$wgt[i]); #  add up wgts of those models
  cat('CumWgt(',s,')=',wgt,'\n')     #  print sum of weights
}

# Deer Occupancy Models in Fall ---------------------------------


nSURVEYsFallOV=ncol(DHFallOV)  #  set number of sites,surveys from det. history data

#create covariates file

load("data/deerDataFall.RData") #I'm using general squirrel data here because the site covariates are the same for both species
covFallOV <- deerDataFall %>%
  select(Deployment_Name, Height_cm, Num_Stems, OakDBH, Log.in.View, Summer.Fall.EDD)

#Created scaled version of the covariates file
covFallOV2 <- covFallOV %>%
  mutate(Height_cm = scale(Height_cm),
         Num_Stems = scale(Num_Stems),
         OakDBH = scale(OakDBH),
         Summer.Fall.EDD = scale(Summer.Fall.EDD))

#Create input file for RPresence
pFallOV <- createPao(DHFallOV,unitcov = covFallOV,title="DeerFall",unitnames=sitenamesShort)
pFallOV2 <- createPao(DHFallOV,unitcov = covFallOV2,title="DeerFall-scaled",unitnames=sitenamesShort)



mods=list(); i=1
mods[[i]]=occMod(model=list(psi~1, p~1),    data=pFallOV2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~Height_cm),data=pFallOV2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~Num_Stems),data=pFallOV2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~OakDBH),data=pFallOV2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~Log.in.View),data=pFallOV2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~Summer.Fall.EDD),data=pFallOV2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~Height_cm + Num_Stems),data=pFallOV2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~Height_cm + OakDBH),data=pFallOV2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~Height_cm + Log.in.View),data=pFallOV2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~Height_cm + Summer.Fall.EDD),data=pFallOV2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~Num_Stems + OakDBH),data=pFallOV2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~Num_Stems + Log.in.View),data=pFallOV2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~Num_Stems + Summer.Fall.EDD),data=pFallOV2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~OakDBH + Log.in.View),data=pFallOV2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~OakDBH + Summer.Fall.EDD),data=pFallOV2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~Log.in.View + Summer.Fall.EDD),data=pFallOV2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~Height_cm + Num_Stems + OakDBH),data=pFallOV2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~Height_cm + Num_Stems + Log.in.View),data=pFallOV2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~Height_cm + Num_Stems + Summer.Fall.EDD),data=pFallOV2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~Height_cm + OakDBH + Log.in.View),data=pFallOV2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~Height_cm + OakDBH + Summer.Fall.EDD),data=pFallOV2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~Height_cm + Log.in.View + Summer.Fall.EDD),data=pFallOV2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~Num_Stems + OakDBH + Log.in.View),data=pFallOV2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~Num_Stems + OakDBH + Summer.Fall.EDD),data=pFallOV2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~Num_Stems + Log.in.View + Summer.Fall.EDD),data=pFallOV2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~OakDBH + Log.in.View + Summer.Fall.EDD),data=pFallOV2,type="so");i=i+1


#     create AIC table of model results and print
resultsFallOV <- createAicTable(mods)
resultsFallOV$table

bestFallOV <- occMod(model = list(psi ~ 1, p ~ 1),data = pFallOV2,type = "so") 
bestFallOV$beta$p

#Sum up model weights for each of the 5 covariates
mnames=resultsFallOV$table$Model;
for (s in c("Height","Oak", "Stems", "Log", "EDD")) {
  i <- grep(s,mnames); #  get model names with each covariate name
  wgt=sum(resultsFallOV$table$wgt[i]); #  add up wgts of those models
  cat('CumWgt(',s,')=',wgt,'\n')     #  print sum of weights
}
