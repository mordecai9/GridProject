#Assessment of optimal detection probalitity models for Raccoons using occupancy models in Rpresence. Data comes from high resolution 27 camera grid is Posey Hollow. Detection histories were prepared using the script "PrepforDetectionHistory", where R objects of detection histories were created for each species and each season. These are loaded in here. The script also requires data files prepared for a separate GLM analysis, since these include the site covariates. 


# Clear Environment and Load Packages -------------------------------------

rm(list = ls())
library(tidyverse)
library(RPresence)

# Load Detection Histories Deer------------------------------------------------

load("results/DetHistFallProcyon lotor")
DHFallPL <- temp

load("results/DetHistSumProcyon lotor")
DHSumPL <- temp

load("results/DetHistWinProcyon lotor")
DHWinPL <- temp

load("results/DetHistSpProcyon lotor")
DHSpPL <- temp

rm(temp)

# Raccoon Occupancy Models in Spring ---------------------------------

sitenamesFull <- row.names(as.data.frame(DHSpPL)) 
nsitesFull=nrow(DHSpPL) 
nSURVEYsSpPL=ncol(DHSpPL)  #  set number of sites,surveys from det. history data

#create covariates file

load("data/racDataSpr.RData")
covSpPL <- racDataSpr %>%
  select(Deployment_Name, Height_cm, Num_Stems, OakDBH, Log.in.View, Raccoon.EDD)

#Created scaled version of the covariates file. 
covSpPL2 <- covSpPL %>%
  mutate(Height_cm = scale(Height_cm),
         Num_Stems = scale(Num_Stems),
         OakDBH = scale(OakDBH),
         Raccoon.EDD = scale(Raccoon.EDD))

#Create input file for RPresence
pSpPL <- createPao(DHSpPL,unitcov = covSpPL,title="RaccoonSpring",unitnames=sitenamesFull)
pSpPL2 <- createPao(DHSpPL,unitcov = covSpPL2,title="RaccoonSpring - scaled",unitnames=sitenamesFull)


mods=list(); i=1
mods[[i]]=occMod(model=list(psi~1, p~1),    data=pSpPL2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~Height_cm),data=pSpPL2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~Num_Stems),data=pSpPL2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~OakDBH),data=pSpPL2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~Log.in.View),data=pSpPL2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~Raccoon.EDD),data=pSpPL2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~Height_cm + Num_Stems),data=pSpPL2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~Height_cm + OakDBH),data=pSpPL2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~Height_cm + Log.in.View),data=pSpPL2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~Height_cm + Raccoon.EDD),data=pSpPL2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~Num_Stems + OakDBH),data=pSpPL2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~Num_Stems + Log.in.View),data=pSpPL2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~Num_Stems + Raccoon.EDD),data=pSpPL2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~OakDBH + Log.in.View),data=pSpPL2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~OakDBH + Raccoon.EDD),data=pSpPL2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~Log.in.View + Raccoon.EDD),data=pSpPL2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~Height_cm + Num_Stems + OakDBH),data=pSpPL2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~Height_cm + Num_Stems + Log.in.View),data=pSpPL2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~Height_cm + Num_Stems + Raccoon.EDD),data=pSpPL2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~Height_cm + OakDBH + Log.in.View),data=pSpPL2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~Height_cm + OakDBH + Raccoon.EDD),data=pSpPL2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~Height_cm + Log.in.View + Raccoon.EDD),data=pSpPL2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~Num_Stems + OakDBH + Log.in.View),data=pSpPL2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~Num_Stems + OakDBH + Raccoon.EDD),data=pSpPL2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~Num_Stems + Log.in.View + Raccoon.EDD),data=pSpPL2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~OakDBH + Log.in.View + Raccoon.EDD),data=pSpPL2,type="so");i=i+1


#     create AIC table of model results and print
resultsSpPL <- createAicTable(mods)
resultsSpPL$table

bestSpPL <- occMod(model=list(psi~1, p~OakDBH),data=pSpPL2,type="so") 
bestSpPL$beta$p

#Average detection probability (Oak)
MeanpRac <- predict(bestSpPL, newdata = data.frame(OakDBH = mean(covSpPL2$OakDBH)), param = "p", conf= 0.95)
row.names(MeanpRac) <- "Raccoon Spring"

#Sum up model weights for each of the 5 covariates
mnames=resultsSpPL$table$Model;
for (s in c("Height","Oak", "Stems", "Log", "EDD")) {
  i <- grep(s,mnames); #  get model names with each covariate name
  wgt=sum(resultsSpPL$table$wgt[i]); #  add up wgts of those models
  cat('CumWgt(',s,')=',wgt,'\n')     #  print sum of weights
}


# Raccoon Occupancy Models in Winter ---------------------------------


sitenamesShort <- row.names(as.data.frame(DHWinPL)) 
nsitesShort=nrow(DHWinPL) 
nSURVEYsWinPL=ncol(DHWinPL)  #  set number of sites,surveys from det. history data

#create covariates file

load("data/racDataWin.RData") #I'm using general squirrel data here because the site covariates are the same for both species
covWinPL <- racDataWin %>%
  select(Deployment_Name, Height_cm, Num_Stems, OakDBH, Log.in.View, Raccoon.EDD)
covWinPL2 <- covWinPL %>%
  mutate(Height_cm = scale(Height_cm),
         Num_Stems = scale(Num_Stems),
         OakDBH = scale(OakDBH),
         Raccoon.EDD = scale(Raccoon.EDD))

#Create input file for RPresence
pWinPL <- createPao(DHWinPL,unitcov = covWinPL,title="RaccoonWinter",unitnames=sitenamesShort)
pWinPL2 <- createPao(DHWinPL,unitcov = covWinPL2,title="RaccoonWinter - Scaled",unitnames=sitenamesShort)


mods=list(); i=1
mods[[i]]=occMod(model=list(psi~1, p~1),    data=pWinPL2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~Height_cm),data=pWinPL2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~Num_Stems),data=pWinPL2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~OakDBH),data=pWinPL2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~Log.in.View),data=pWinPL2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~Raccoon.EDD),data=pWinPL2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~Height_cm + Num_Stems),data=pWinPL2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~Height_cm + OakDBH),data=pWinPL2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~Height_cm + Log.in.View),data=pWinPL2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~Height_cm + Raccoon.EDD),data=pWinPL2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~Num_Stems + OakDBH),data=pWinPL2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~Num_Stems + Log.in.View),data=pWinPL2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~Num_Stems + Raccoon.EDD),data=pWinPL2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~OakDBH + Log.in.View),data=pWinPL2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~OakDBH + Raccoon.EDD),data=pWinPL2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~Log.in.View + Raccoon.EDD),data=pWinPL2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~Height_cm + Num_Stems + OakDBH),data=pWinPL2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~Height_cm + Num_Stems + Log.in.View),data=pWinPL2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~Height_cm + Num_Stems + Raccoon.EDD),data=pWinPL2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~Height_cm + OakDBH + Log.in.View),data=pWinPL2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~Height_cm + OakDBH + Raccoon.EDD),data=pWinPL2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~Height_cm + Log.in.View + Raccoon.EDD),data=pWinPL2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~Num_Stems + OakDBH + Log.in.View),data=pWinPL2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~Num_Stems + OakDBH + Raccoon.EDD),data=pWinPL2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~Num_Stems + Log.in.View + Raccoon.EDD),data=pWinPL2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~OakDBH + Log.in.View + Raccoon.EDD),data=pWinPL2,type="so");i=i+1


#     create AIC table of model results and print
resultsWinPL <- createAicTable(mods)
resultsWinPL$table

bestWinPL <- occMod(model=list(psi~1, p ~ Height_cm + Num_Stems + Log.in.View),data=pWinPL2,type="so")
bestWinPL$beta$p

#Average detection probability (Stems, Log, Height_cm)
MeanpRac[2,] <- predict(bestWinPL, newdata = data.frame(Num_Stems = mean(covWinPL2$Num_Stems), Log.in.View = factor("NO", levels=c("NO", "YES")), Height_cm = mean(covWinPL2$Height_cm)), param = "p", conf= 0.95)
row.names(MeanpRac)[2] <- "Raccoon Winter"


#Sum up model weights for each of the 5 covariates
mnames=resultsWinPL$table$Model;
for (s in c("Height","Oak", "Stems", "Log", "EDD")) {
  i <- grep(s,mnames); #  get model names with each covariate name
  wgt=sum(resultsWinPL$table$wgt[i]); #  add up wgts of those models
  cat('CumWgt(',s,')=',wgt,'\n')     #  print sum of weights
}

# Raccoon Occupancy Models in Summer ---------------------------------


nSURVEYsSumPL=ncol(DHSumPL)  #  set number of sites,surveys from det. history data

#create covariates file

load("data/racDataSum.RData") #I'm using general squirrel data here because the site covariates are the same for both species
covSumPL <- racDataSum %>%
  select(Deployment_Name, Height_cm, Num_Stems, OakDBH, Log.in.View, Raccoon.EDD)

#Created scaled version of the covariates file. 
covSumPL2 <- covSumPL %>%
  mutate(Height_cm = scale(Height_cm),
         Num_Stems = scale(Num_Stems),
         OakDBH = scale(OakDBH),
         Raccoon.EDD = scale(Raccoon.EDD))

#Create input file for RPresence
pSumPL <- createPao(DHSumPL,unitcov = covSumPL,title="RaccoonSummer",unitnames=sitenamesFull)
pSumPL2 <- createPao(DHSumPL,unitcov = covSumPL2,title="RaccoonSummer-scaled",unitnames=sitenamesFull)


mods=list(); i=1
mods[[i]]=occMod(model=list(psi~1, p~1),    data=pSumPL2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~Height_cm),data=pSumPL2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~Num_Stems),data=pSumPL2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~OakDBH),data=pSumPL2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~Log.in.View),data=pSumPL2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~Raccoon.EDD),data=pSumPL2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~Height_cm + Num_Stems),data=pSumPL2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~Height_cm + OakDBH),data=pSumPL2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~Height_cm + Log.in.View),data=pSumPL2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~Height_cm + Raccoon.EDD),data=pSumPL2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~Num_Stems + OakDBH),data=pSumPL2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~Num_Stems + Log.in.View),data=pSumPL2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~Num_Stems + Raccoon.EDD),data=pSumPL2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~OakDBH + Log.in.View),data=pSumPL2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~OakDBH + Raccoon.EDD),data=pSumPL2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~Log.in.View + Raccoon.EDD),data=pSumPL2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~Height_cm + Num_Stems + OakDBH),data=pSumPL2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~Height_cm + Num_Stems + Log.in.View),data=pSumPL2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~Height_cm + Num_Stems + Raccoon.EDD),data=pSumPL2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~Height_cm + OakDBH + Log.in.View),data=pSumPL2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~Height_cm + OakDBH + Raccoon.EDD),data=pSumPL2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~Height_cm + Log.in.View + Raccoon.EDD),data=pSumPL2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~Num_Stems + OakDBH + Log.in.View),data=pSumPL2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~Num_Stems + OakDBH + Raccoon.EDD),data=pSumPL2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~Num_Stems + Log.in.View + Raccoon.EDD),data=pSumPL2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~OakDBH + Log.in.View + Raccoon.EDD),data=pSumPL2,type="so");i=i+1


#     create AIC table of model results and print
resultsSumPL <- createAicTable(mods)
resultsSumPL$table

bestSumPL <- occMod(model = list(psi~1, p ~ Num_Stems),data = pSumPL2,type = "so")
bestSumPL$beta$p

#Average detection probability (Height, Log, Squirrel EDD 4S)
MeanpRac[3,] <- predict(bestSumPL, newdata = data.frame(Num_Stems = mean(covSumPL2$Num_Stems)), param = "p", conf= 0.95)
row.names(MeanpRac)[3] <- "Raccoon Summer"

#Sum up model weights for each of the 5 covariates
mnames=resultsSumPL$table$Model;
for (s in c("Height","Oak", "Stems", "Log", "EDD")) {
  i <- grep(s,mnames); #  get model names with each covariate name
  wgt=sum(resultsSumPL$table$wgt[i]); #  add up wgts of those models
  cat('CumWgt(',s,')=',wgt,'\n')     #  print sum of weights
}

# Raccoon Occupancy Models in Fall ---------------------------------


nSURVEYsFallPL=ncol(DHFallPL)  #  set number of sites,surveys from det. history data

#create covariates file

load("data/racDataFall.RData") #I'm using general squirrel data here because the site covariates are the same for both species
covFallPL <- racDataFall %>%
  select(Deployment_Name, Height_cm, Num_Stems, OakDBH, Log.in.View, Raccoon.EDD)

#Created scaled version of the covariates file. 
covFallPL2 <- covFallPL %>%
  mutate(Height_cm = scale(Height_cm),
         Num_Stems = scale(Num_Stems),
         OakDBH = scale(OakDBH),
         Raccoon.EDD = scale(Raccoon.EDD))


#Create input file for RPresence
pFallPL <- createPao(DHFallPL,unitcov = covFallPL,title="RaccoonFall",unitnames=sitenamesShort)
pFallPL2 <- createPao(DHFallPL,unitcov = covFallPL2,title="RaccoonFall-scaled",unitnames=sitenamesShort)



mods=list(); i=1
mods[[i]]=occMod(model=list(psi~1, p~1),    data=pFallPL2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~Height_cm),data=pFallPL2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~Num_Stems),data=pFallPL2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~OakDBH),data=pFallPL2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~Log.in.View),data=pFallPL2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~Raccoon.EDD),data=pFallPL2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~Height_cm + Num_Stems),data=pFallPL2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~Height_cm + OakDBH),data=pFallPL2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~Height_cm + Log.in.View),data=pFallPL2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~Height_cm + Raccoon.EDD),data=pFallPL2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~Num_Stems + OakDBH),data=pFallPL2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~Num_Stems + Log.in.View),data=pFallPL2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~Num_Stems + Raccoon.EDD),data=pFallPL2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~OakDBH + Log.in.View),data=pFallPL2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~OakDBH + Raccoon.EDD),data=pFallPL2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~Log.in.View + Raccoon.EDD),data=pFallPL2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~Height_cm + Num_Stems + OakDBH),data=pFallPL2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~Height_cm + Num_Stems + Log.in.View),data=pFallPL2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~Height_cm + Num_Stems + Raccoon.EDD),data=pFallPL2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~Height_cm + OakDBH + Log.in.View),data=pFallPL2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~Height_cm + OakDBH + Raccoon.EDD),data=pFallPL2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~Height_cm + Log.in.View + Raccoon.EDD),data=pFallPL2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~Num_Stems + OakDBH + Log.in.View),data=pFallPL2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~Num_Stems + OakDBH + Raccoon.EDD),data=pFallPL2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~Num_Stems + Log.in.View + Raccoon.EDD),data=pFallPL2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~OakDBH + Log.in.View + Raccoon.EDD),data=pFallPL2,type="so");i=i+1


#     create AIC table of model results and print
resultsFallPL <- createAicTable(mods)
resultsFallPL$table

bestFallPL <- occMod(model = list(psi ~ 1, p ~ Num_Stems + OakDBH + Log.in.View),data = pFallPL2,type = "so") 
bestFallPL$beta$p

#Average detection probability (Stems, Log, Oak)
MeanpRac[4,] <- predict(bestFallPL, newdata = data.frame(Log.in.View = factor("NO", levels=c("NO", "YES")), Num_Stems = mean(covFallPL2$Num_Stems), OakDBH = mean(covFallPL2$OakDBH)), param = "p", conf= 0.95)
row.names(MeanpRac)[4] <- "Raccoon Fall"
save(MeanpRac, file = "results/MeanPSeasons_Rac")

#Sum up model weights for each of the 5 covariates
mnames=resultsFallPL$table$Model;
for (s in c("Height","Oak", "Stems", "Log", "EDD")) {
  i <- grep(s,mnames); #  get model names with each covariate name
  wgt=sum(resultsFallPL$table$wgt[i]); #  add up wgts of those models
  cat('CumWgt(',s,')=',wgt,'\n')
  }#  print sum of weights