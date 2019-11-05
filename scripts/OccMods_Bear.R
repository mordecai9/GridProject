#Assessment of optimal detection probalitity models for each of Deerusing occupancy models in Rpresence. Data comes from high resolution 27 camera grid is Posey Hollow. Detection histories were prepared using the OCript "PrepforDetectionHistory", where R objects of detection hustories were created for each species and each season. These are loaded in here. The script also requires data files prepared for a separate GLM analysis, since these include the site covariates. 


# Clear Environment and Load Packages -------------------------------------

rm(list = ls())
library(tidyverse)
library(RPresence)

# Load Detection Histories Bear------------------------------------------------

load("results/DetHistFallUrsus americanus")
DHFallUA <- temp

load("results/DetHistSumUrsus americanus")
DHSumUA <- temp

rm(temp)


# Bear Occupancy Models in Summer ---------------------------------

sitenamesFull <- row.names(as.data.frame(DHSumUA)) 
nSURVEYsSumUA=ncol(DHSumUA)  #  set number of sites,surveys from det. history data

#create covariates file

load("data/deerDataSum.RData") #I'm using general deer data here because the site covariates are the same for deer and bear species
covSumUA <- deerDataSum %>%
  select(Deployment_Name, Height_cm, Num_Stems, OakDBH, Log.in.View, Summer.Fall.EDD)

#Created scaled version of the covariates file
covSumUA2 <- covSumUA %>%
  mutate(Height_cm = scale(Height_cm),
         Num_Stems = scale(Num_Stems),
         OakDBH = scale(OakDBH),
         Summer.Fall.EDD = scale(Summer.Fall.EDD))

#Create input file for RPresence
pSumUA <- createPao(DHSumUA,unitcov = covSumUA,title="Bear in Summer",unitnames=sitenamesFull)
pSumUA2 <- createPao(DHSumUA,unitcov = covSumUA2,title="BearSummer-scaled",unitnames=sitenamesFull)


mods=list(); i=1
mods[[i]]=occMod(model=list(psi~1, p~1),    data=pSumUA2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~Height_cm),data=pSumUA2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~Num_Stems),data=pSumUA2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~OakDBH),data=pSumUA2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~Log.in.View),data=pSumUA2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~Summer.Fall.EDD),data=pSumUA2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~Height_cm + Num_Stems),data=pSumUA2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~Height_cm + OakDBH),data=pSumUA2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~Height_cm + Log.in.View),data=pSumUA2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~Height_cm + Summer.Fall.EDD),data=pSumUA2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~Num_Stems + OakDBH),data=pSumUA2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~Num_Stems + Log.in.View),data=pSumUA2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~Num_Stems + Summer.Fall.EDD),data=pSumUA2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~OakDBH + Log.in.View),data=pSumUA2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~OakDBH + Summer.Fall.EDD),data=pSumUA2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~Log.in.View + Summer.Fall.EDD),data=pSumUA2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~Height_cm + Num_Stems + OakDBH),data=pSumUA2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~Height_cm + Num_Stems + Log.in.View),data=pSumUA2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~Height_cm + Num_Stems + Summer.Fall.EDD),data=pSumUA2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~Height_cm + OakDBH + Log.in.View),data=pSumUA2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~Height_cm + OakDBH + Summer.Fall.EDD),data=pSumUA2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~Height_cm + Log.in.View + Summer.Fall.EDD),data=pSumUA2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~Num_Stems + OakDBH + Log.in.View),data=pSumUA2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~Num_Stems + OakDBH + Summer.Fall.EDD),data=pSumUA2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~Num_Stems + Log.in.View + Summer.Fall.EDD),data=pSumUA2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~OakDBH + Log.in.View + Summer.Fall.EDD),data=pSumUA2,type="so");i=i+1


#     create AIC table of model results and print
resultsSumUA <- createAicTable(mods)
resultsSumUA$table

bestSumUA <- occMod(model = list(psi~1, p ~ Log.in.View + Summer.Fall.EDD),data = pSumUA2,type = "so")
bestSumUA$beta$p

#Average detection probability (Log in View, Summer.Fall EDD)

MeanpBear <- predict(bestSumUA, newdata = data.frame(Log.in.View = factor("NO", levels=c("NO", "YES")), Summer.Fall.EDD = mean(covSumUA2$Summer.Fall.EDD)), param = "p", conf= 0.95)
row.names(MeanpBear)[1] <- "Bear Summer"

#Sum up model weights for each of the 5 covariates
mnames=resultsSumUA$table$Model;
for (s in c("Height","Oak", "Stems", "Log", "EDD")) {
  i <- grep(s,mnames); #  get model names with each cUAariate name
  wgt=sum(resultsSumUA$table$wgt[i]); #  add up wgts of those models
  cat('CumWgt(',s,')=',wgt,'\n')     #  print sum of weights
}

# Bear Occupancy Models in Fall ---------------------------------

sitenamesShort <- row.names(as.data.frame(DHFallUA)) 
nSURVEYsFallUA=ncol(DHFallUA)  #  set number of sites,surveys from det. history data

#create covariates file

load("data/deerDataFall.RData") #I'm using general deer data here because the site covariates are the same for both species
covFallUA <- deerDataFall %>%
  select(Deployment_Name, Height_cm, Num_Stems, OakDBH, Log.in.View, Summer.Fall.EDD)

#Created scaled version of the covariates file
covFallUA2 <- covFallUA %>%
  mutate(Height_cm = scale(Height_cm),
         Num_Stems = scale(Num_Stems),
         OakDBH = scale(OakDBH),
         Summer.Fall.EDD = scale(Summer.Fall.EDD))
#Create input file for RPresence
pFallUA <- createPao(DHFallUA,unitcov = covFallUA,title="BearFall",unitnames=sitenamesShort)
pFallUA2 <- createPao(DHFallUA,unitcov = covFallUA2,title="BearFall-scaled",unitnames=sitenamesShort)



mods=list(); i=1
mods[[i]]=occMod(model=list(psi~1, p~1),    data=pFallUA2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~Height_cm),data=pFallUA2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~Num_Stems),data=pFallUA2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~OakDBH),data=pFallUA2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~Log.in.View),data=pFallUA2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~Summer.Fall.EDD),data=pFallUA2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~Height_cm + Num_Stems),data=pFallUA2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~Height_cm + OakDBH),data=pFallUA2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~Height_cm + Log.in.View),data=pFallUA2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~Height_cm + Summer.Fall.EDD),data=pFallUA2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~Num_Stems + OakDBH),data=pFallUA2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~Num_Stems + Log.in.View),data=pFallUA2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~Num_Stems + Summer.Fall.EDD),data=pFallUA2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~OakDBH + Log.in.View),data=pFallUA2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~OakDBH + Summer.Fall.EDD),data=pFallUA2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~Log.in.View + Summer.Fall.EDD),data=pFallUA2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~Height_cm + Num_Stems + OakDBH),data=pFallUA2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~Height_cm + Num_Stems + Log.in.View),data=pFallUA2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~Height_cm + Num_Stems + Summer.Fall.EDD),data=pFallUA2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~Height_cm + OakDBH + Log.in.View),data=pFallUA2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~Height_cm + OakDBH + Summer.Fall.EDD),data=pFallUA2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~Height_cm + Log.in.View + Summer.Fall.EDD),data=pFallUA2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~Num_Stems + OakDBH + Log.in.View),data=pFallUA2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~Num_Stems + OakDBH + Summer.Fall.EDD),data=pFallUA2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~Num_Stems + Log.in.View + Summer.Fall.EDD),data=pFallUA2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~OakDBH + Log.in.View + Summer.Fall.EDD),data=pFallUA2,type="so");i=i+1


#     create AIC table of model results and print
resultsFallUA <- createAicTable(mods)
resultsFallUA$table

bestFallUA <- occMod(model = list(psi ~ 1, p ~ OakDBH + Summer.Fall.EDD),data = pFallUA2,type = "so") 
bestFallUA$beta$p

#Average detection probability (OakDBH, Summer.Fall EDD)

MeanpBear[2,] <- predict(bestFallUA, newdata = data.frame(OakDBH = mean(covFallUA2$OakDBH), Summer.Fall.EDD = mean(covFallUA2$Summer.Fall.EDD)), param = "p", conf= 0.95)
row.names(MeanpBear)[2] <- "Bear Fall"
save(MeanpBear, file = "results/MeanPSeasonsBear")

#Sum up model weights for each of the 5 covariates
mnames=resultsFallUA$table$Model;
for (s in c("Height","Oak", "Stems", "Log", "EDD")) {
  i <- grep(s,mnames); #  get model names with each covariate name
  wgt=sum(resultsFallUA$table$wgt[i]); #  add up wgts of those models
  cat('CumWgt(',s,')=',wgt,'\n') 
  }#  print sum of weights
