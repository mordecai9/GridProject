#Assessment of optimal detection probalitity models for each of Gray Squirrel using occupancy models in Rpresence. Data comes from high resolution 27 camera grid is Posey Hollow. Detection histories were prepared using the script "PrepforDetectionHistory", where R objects of detection hustories were created for each species and each season. These are loaded in here. The script also requires data files prepared for a separate GLM analysis, since these include the site covariates. 


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
#Created scaled version of the covariates file
covSpSC2 <- covSpSC %>%
  mutate(Height_cm = scale(Height_cm),
         Num_Stems = scale(Num_Stems),
         OakDBH = scale(OakDBH),
         Squirrel_EDD_WSp = scale(Squirrel_EDD_WSp))

#Create input file for RPresence
pSpSC <- createPao(DHSpSC,unitcov = covSpSC,title="Gray Squirrel Spring",unitnames=sitenamesFull)
pSpSC2 <- createPao(DHSpSC,unitcov = covSpSC2,title="Gray Squirrel Spring - scaled",unitnames=sitenamesFull)


#Test to see if quadratic effect on height improves the height only model

modsHtest <- list(); 
modsHtest[[1]] <- occMod(model=list(psi ~ 1, p ~ Height_cm),data = pSpSC2,type = "so")
modsHtest[[2]] <- occMod(model=list(psi ~ 1, p ~ poly(Height_cm, 2, raw = T)),data = pSpSC2,type = "so")

resultsHtest <- createAicTable(modsHtest)
resultsHtest$table #there is evidence of a quadratic effect on camera height, so it should be in all models with camera height

#Test to see if there is an interaction between EDD and whether there is a log in view

modsEDDtest <- list(); 
modsEDDtest[[1]] <- occMod(model=list(psi ~ 1, p ~ Log.in.View + Squirrel_EDD_WSp),data = pSpSC2,type = "so")
modsEDDtest[[2]] <- occMod(model=list(psi ~ 1, p ~ Log.in.View * Squirrel_EDD_WSp),data = pSpSC2,type = "so")

resultsEDDtest <- createAicTable(modsEDDtest)
resultsEDDtest$table #there is evidence of a interaction effect here, so this needs to be in all models below with both of these parameters


mods=list(); i=1
mods[[i]]=occMod(model=list(psi~1, p~1),    data=pSpSC2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~poly(Height_cm, 2, raw = T)),data=pSpSC2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~Num_Stems),data=pSpSC2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~OakDBH),data=pSpSC2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~Log.in.View),data=pSpSC2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~Squirrel_EDD_WSp),data=pSpSC2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~poly(Height_cm, 2, raw = T) + Num_Stems),data=pSpSC2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~poly(Height_cm, 2, raw = T) + OakDBH),data=pSpSC2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~poly(Height_cm, 2, raw = T) + Log.in.View),data=pSpSC2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~poly(Height_cm, 2, raw = T) + Squirrel_EDD_WSp),data=pSpSC2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~Num_Stems + OakDBH),data=pSpSC2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~Num_Stems + Log.in.View),data=pSpSC2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~Num_Stems + Squirrel_EDD_WSp),data=pSpSC2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~OakDBH + Log.in.View),data=pSpSC2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~OakDBH + Squirrel_EDD_WSp),data=pSpSC2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~Log.in.View * Squirrel_EDD_WSp),data=pSpSC2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~poly(Height_cm, 2, raw = T) + Num_Stems + OakDBH),data=pSpSC2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~poly(Height_cm, 2, raw = T) + Num_Stems + Log.in.View),data=pSpSC2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~poly(Height_cm, 2, raw = T) + Num_Stems + Squirrel_EDD_WSp),data=pSpSC2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~poly(Height_cm, 2, raw = T) + OakDBH + Log.in.View),data=pSpSC2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~poly(Height_cm, 2, raw = T) + OakDBH + Squirrel_EDD_WSp),data=pSpSC2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~poly(Height_cm, 2, raw = T) + Log.in.View * Squirrel_EDD_WSp),data=pSpSC2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~Num_Stems + OakDBH + Log.in.View),data=pSpSC2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~Num_Stems + OakDBH + Squirrel_EDD_WSp),data=pSpSC2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~Num_Stems + Log.in.View * Squirrel_EDD_WSp),data=pSpSC2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~OakDBH + Log.in.View * Squirrel_EDD_WSp),data=pSpSC2,type="so");i=i+1


#     create AIC table of model results and print
resultsSpSC <- createAicTable(mods)
resultsSpSC$table

bestSpSC <- occMod(model=list(psi~1, p~OakDBH + Log.in.View * Squirrel_EDD_WSp),data=pSpSC2,type="so") 
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
#scaled covariate set
covWinSC2 <- covWinSC %>%
  mutate(Height_cm = scale(Height_cm),
         Num_Stems = scale(Num_Stems),
         OakDBH = scale(OakDBH),
         Squirrel_EDD_WSp = scale(Squirrel_EDD_WSp))
#Create input file for RPresence
pWinSC <- createPao(DHWinSC,unitcov = covWinSC,title="Gray Squirrel Winter",unitnames=sitenamesShort)
pWinSC2 <- createPao(DHWinSC,unitcov = covWinSC2,title="Gray Squirrel Winter - scaled",unitnames=sitenamesShort)

#Test to see if quadratic effect on height improves the height only model

modsHtest <- list(); 
modsHtest[[1]] <- occMod(model=list(psi ~ 1, p ~ Height_cm),data = pWinSC2,type = "so")
modsHtest[[2]] <- occMod(model=list(psi ~ 1, p ~ poly(Height_cm, 2, raw = T)),data = pWinSC2,type = "so")

resultsHtest <- createAicTable(modsHtest)
resultsHtest$table #there is evidence of a quadratic effect on camera height, so this needs to be in models below

#Test to see if there is an interaction between EDD and whether there is a log in view

modsEDDtest <- list(); 
modsEDDtest[[1]] <- occMod(model=list(psi ~ 1, p ~ Log.in.View + Squirrel_EDD_WSp),data = pWinSC2,type = "so")
modsEDDtest[[2]] <- occMod(model=list(psi ~ 1, p ~ Log.in.View * Squirrel_EDD_WSp),data = pWinSC2,type = "so")

resultsEDDtest <- createAicTable(modsEDDtest)
resultsEDDtest$table #there is evidence of a interaction effect here, so this needs to be in all models below with both of these parameters


mods=list(); i=1
mods[[i]]=occMod(model=list(psi~1, p~1),    data=pWinSC2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~poly(Height_cm, 2, raw = T)),data=pWinSC2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~Num_Stems),data=pWinSC2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~OakDBH),data=pWinSC2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~Log.in.View),data=pWinSC2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~Squirrel_EDD_WSp),data=pWinSC2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~poly(Height_cm, 2, raw = T) + Num_Stems),data=pWinSC2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~poly(Height_cm, 2, raw = T) + OakDBH),data=pWinSC2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~poly(Height_cm, 2, raw = T) + Log.in.View),data=pWinSC2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~poly(Height_cm, 2, raw = T) + Squirrel_EDD_WSp),data=pWinSC2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~Num_Stems + OakDBH),data=pWinSC2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~Num_Stems + Log.in.View),data=pWinSC2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~Num_Stems + Squirrel_EDD_WSp),data=pWinSC2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~OakDBH + Log.in.View),data=pWinSC2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~OakDBH + Squirrel_EDD_WSp),data=pWinSC2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~Log.in.View * Squirrel_EDD_WSp),data=pWinSC2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~poly(Height_cm, 2, raw = T) + Num_Stems + OakDBH),data=pWinSC2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~poly(Height_cm, 2, raw = T) + Num_Stems + Log.in.View),data=pWinSC2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~poly(Height_cm, 2, raw = T) + Num_Stems + Squirrel_EDD_WSp),data=pWinSC2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~poly(Height_cm, 2, raw = T) + OakDBH + Log.in.View),data=pWinSC2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~poly(Height_cm, 2, raw = T) + OakDBH + Squirrel_EDD_WSp),data=pWinSC2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~poly(Height_cm, 2, raw = T) + Log.in.View * Squirrel_EDD_WSp),data=pWinSC2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~Num_Stems + OakDBH + Log.in.View),data=pWinSC2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~Num_Stems + OakDBH + Squirrel_EDD_WSp),data=pWinSC2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~Num_Stems + Log.in.View * Squirrel_EDD_WSp),data=pWinSC2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~OakDBH + Log.in.View * Squirrel_EDD_WSp),data=pWinSC2,type="so");i=i+1


#     create AIC table of model results and print
resultsWinSC <- createAicTable(mods)
resultsWinSC$table

bestWinSC <- occMod(model=list(psi~1, p ~ poly(Height_cm, 2, raw = T) + Log.in.View * Squirrel_EDD_WSp),data=pWinSC2,type="so")
bestWinSC$beta$p

#Sum up model weights for each of the 5 covariates
mnames=resultsWinSC$table$Model;
for (s in c("Height","Oak", "Stems", "Log", "EDD")) {
  i <- grep(s,mnames); #  get model names with each covariate name
  wgt=sum(resultsWinSC$table$wgt[i]); #  add up wgts of those models
  cat('CumWgt(',s,')=',wgt,'\n')     #  print sum of weights
}

# Response Curve - Height Gray Sq Winter -----------------------------------
#Could think about doing two curves, one with and one without logs in view
seq.WHgt <- seq(min(covWinSC2$Height_cm), max(covWinSC2$Height_cm), length.out = 100)

nd.seq.WHgt <- data.frame(Height_cm = seq.WHgt, Squirrel_EDD_WSp = median(covWinSC2$Squirrel_EDD_WSp),Log.in.View = factor("YES", levels=c("NO", "YES")) )

pred.WHgt <- predict(bestWinSC, newdata = nd.seq.WHgt, param = "p", conf= 0.95)

predPlotWHgt <- data.frame(nd.seq.WHgt, pred.WHgt) 

load(file = "data/responseTheme.Rdata")

GSqWinHgtPlot <- ggplot(predPlotWHgt, aes(x=(Height_cm*sd(covWinSC$Height_cm))+mean(covWinSC$Height_cm), y=est)) +
  # Confidence region
  geom_ribbon(aes(ymin=lower_0.95, ymax=upper_0.95), alpha=0.25, show.legend = F) +
  # Prediction Lines
  geom_line() +
  xlab("Camera height (cm)") +   
  ylab("Estimated detection probability") +
  geom_text(aes(x = 34, y = .40, label = "A"), size = 8)+
  myTheme

# Gray Squirrel Occupancy Models in Summer ---------------------------------


nSURVEYsSumSC=ncol(DHSumSC)  #  set number of sites,surveys from det. history data

#create covariates file

load("data/sqDataSum.RData") #I'm using general squirrel data here because the site covariates are the same for both species
covSumSC <- sqDataSum %>%
  select(Deployment_Name, Height_cm, Num_Stems, OakDBH, Log.in.View, Squirrel_EDD_4S)

#Created scaled version of the covariates file
covSumSC2 <- covSumSC %>%
  mutate(Height_cm = scale(Height_cm),
         Num_Stems = scale(Num_Stems),
         OakDBH = scale(OakDBH),
         Squirrel_EDD_4S = scale(Squirrel_EDD_4S))

#Create input file for RPresence
pSumSC <- createPao(DHSumSC,unitcov = covSumSC,title="Gray Squirrel Summer",unitnames=sitenamesFull)
pSumSC2 <- createPao(DHSumSC,unitcov = covSumSC2,title="Gray Squirrel Summer",unitnames=sitenamesFull)

#Test to see if quadratic effect on height improves the height only model

modsHtest <- list(); 
modsHtest[[1]] <- occMod(model=list(psi ~ 1, p ~ Height_cm),data = pSumSC2,type = "so")
modsHtest[[2]] <- occMod(model=list(psi ~ 1, p ~ poly(Height_cm, 2, raw = T)),data = pSumSC2,type = "so")

resultsHtest <- createAicTable(modsHtest)
resultsHtest$table #there is no evidence of a quadratic effect on camera height here

#Test to see if there is an interaction between EDD and whether there is a log in view

modsEDDtest <- list(); 
modsEDDtest[[1]] <- occMod(model=list(psi ~ 1, p ~ Log.in.View + Squirrel_EDD_4S),data = pSumSC2,type = "so")
modsEDDtest[[2]] <- occMod(model=list(psi ~ 1, p ~ Log.in.View * Squirrel_EDD_4S),data = pSumSC2,type = "so")

resultsEDDtest <- createAicTable(modsEDDtest)
resultsEDDtest$table #there is no evidence of a interaction effect here

mods=list(); i=1
mods[[i]]=occMod(model=list(psi~1, p~1),    data=pSumSC2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~Height_cm),data=pSumSC2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~Num_Stems),data=pSumSC2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~OakDBH),data=pSumSC2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~Log.in.View),data=pSumSC2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~Squirrel_EDD_4S),data=pSumSC2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~Height_cm + Num_Stems),data=pSumSC2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~Height_cm + OakDBH),data=pSumSC2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~Height_cm + Log.in.View),data=pSumSC2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~Height_cm + Squirrel_EDD_4S),data=pSumSC2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~Num_Stems + OakDBH),data=pSumSC2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~Num_Stems + Log.in.View),data=pSumSC2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~Num_Stems + Squirrel_EDD_4S),data=pSumSC2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~OakDBH + Log.in.View),data=pSumSC2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~OakDBH + Squirrel_EDD_4S),data=pSumSC2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~Log.in.View + Squirrel_EDD_4S),data=pSumSC2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~Height_cm + Num_Stems + OakDBH),data=pSumSC2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~Height_cm + Num_Stems + Log.in.View),data=pSumSC2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~Height_cm + Num_Stems + Squirrel_EDD_4S),data=pSumSC2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~Height_cm + OakDBH + Log.in.View),data=pSumSC2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~Height_cm + OakDBH + Squirrel_EDD_4S),data=pSumSC2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~Height_cm + Log.in.View + Squirrel_EDD_4S),data=pSumSC2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~Num_Stems + OakDBH + Log.in.View),data=pSumSC2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~Num_Stems + OakDBH + Squirrel_EDD_4S),data=pSumSC2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~Num_Stems + Log.in.View + Squirrel_EDD_4S),data=pSumSC2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~OakDBH + Log.in.View + Squirrel_EDD_4S),data=pSumSC2,type="so");i=i+1


#     create AIC table of model results and print
resultsSumSC <- createAicTable(mods)
resultsSumSC$table

bestSumSC <- occMod(model = list(psi~1, p ~ OakDBH + Log.in.View + Squirrel_EDD_4S),data = pSumSC2,type = "so")
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

#Created scaled version of the covariates file
covFallSC2 <- covFallSC %>%
  mutate(Height_cm = scale(Height_cm),
         Num_Stems = scale(Num_Stems),
         OakDBH = scale(OakDBH),
         Squirrel_EDD_4S = scale(Squirrel_EDD_4S))

#Create input file for RPresence
pFallSC <- createPao(DHFallSC,unitcov = covFallSC,title="Gray Squirrel Fall",unitnames=sitenamesShort)
pFallSC2 <- createPao(DHFallSC,unitcov = covFallSC2,title="Gray Squirrel Fall - scaled",unitnames=sitenamesShort)

#Test to see if quadratic effect on height improves the height only model

modsHtest <- list(); 
modsHtest[[1]] <- occMod(model=list(psi ~ 1, p ~ Height_cm),data = pFallSC2,type = "so")
modsHtest[[2]] <- occMod(model=list(psi ~ 1, p ~ poly(Height_cm, 2, raw = T)),data = pFallSC2,type = "so")

resultsHtest <- createAicTable(modsHtest)
resultsHtest$table #there is no evidence of a quadratic effect on camera height here

#Test to see if there is an interaction between EDD and whether there is a log in view

modsEDDtest <- list(); 
modsEDDtest[[1]] <- occMod(model=list(psi ~ 1, p ~ Log.in.View + Squirrel_EDD_4S),data = pFallSC2,type = "so")
modsEDDtest[[2]] <- occMod(model=list(psi ~ 1, p ~ Log.in.View * Squirrel_EDD_4S),data = pFallSC2,type = "so")

resultsEDDtest <- createAicTable(modsEDDtest)
resultsEDDtest$table #there is evidence of a interaction effect here so it needs to be reflected in models below


mods=list(); i=1
mods[[i]]=occMod(model=list(psi~1, p~1),    data=pFallSC2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~Height_cm),data=pFallSC2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~Num_Stems),data=pFallSC2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~OakDBH),data=pFallSC2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~Log.in.View),data=pFallSC2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~Squirrel_EDD_4S),data=pFallSC2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~Height_cm + Num_Stems),data=pFallSC2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~Height_cm + OakDBH),data=pFallSC2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~Height_cm + Log.in.View),data=pFallSC2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~Height_cm + Squirrel_EDD_4S),data=pFallSC2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~Num_Stems + OakDBH),data=pFallSC2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~Num_Stems + Log.in.View),data=pFallSC2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~Num_Stems + Squirrel_EDD_4S),data=pFallSC2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~OakDBH + Log.in.View),data=pFallSC2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~OakDBH + Squirrel_EDD_4S),data=pFallSC2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~Log.in.View * Squirrel_EDD_4S),data=pFallSC2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~Height_cm + Num_Stems + OakDBH),data=pFallSC2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~Height_cm + Num_Stems + Log.in.View),data=pFallSC2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~Height_cm + Num_Stems + Squirrel_EDD_4S),data=pFallSC2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~Height_cm + OakDBH + Log.in.View),data=pFallSC2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~Height_cm + OakDBH + Squirrel_EDD_4S),data=pFallSC2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~Height_cm + Log.in.View * Squirrel_EDD_4S),data=pFallSC2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~Num_Stems + OakDBH + Log.in.View),data=pFallSC2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~Num_Stems + OakDBH + Squirrel_EDD_4S),data=pFallSC2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~Num_Stems + Log.in.View * Squirrel_EDD_4S),data=pFallSC2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~OakDBH + Log.in.View * Squirrel_EDD_4S),data=pFallSC2,type="so");i=i+1


#     create AIC table of model results and print
resultsFallSC <- createAicTable(mods)
resultsFallSC$table

bestFallSC <- occMod(model=list(psi~1, p~OakDBH + Log.in.View * Squirrel_EDD_4S),data=pFallSC2,type="so") 
bestFallSC$beta$p

#Sum up model weights for each of the 5 covariates
mnames=resultsFallSC$table$Model;
for (s in c("Height","Oak", "Stems", "Log", "EDD")) {
  i <- grep(s,mnames); #  get model names with each covariate name
  wgt=sum(resultsFallSC$table$wgt[i]); #  add up wgts of those models
  cat('CumWgt(',s,')=',wgt,'\n')     #  print sum of weights
}
