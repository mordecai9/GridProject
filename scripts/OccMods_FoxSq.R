#Assesement of optimal detection probability models for Fox Squirrel using occupancy models in Rpresence. Data comes from high resolution 27 camera grid is Posey Hollow. Detection histories were prepared using the script "PrepforDetectionHistory", where R objects of detection histories were created for each species and each season. These are loaded in here. The script also requires data files prepared for a separate GLM analysis, since these include the site covariates. 


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
  dplyr::select(Deployment_Name, Height_cm, Num_Stems, OakDBH, Log.in.View, Squirrel_EDD_WSp)

#Created scaled version of the covariates file
covSpSN2 <- covSpSN %>%
  mutate(Height_cm = scale(Height_cm),
         Num_Stems = scale(Num_Stems),
         OakDBH = scale(OakDBH),
         Squirrel_EDD_WSp = scale(Squirrel_EDD_WSp))

#Create input file for RPresence and scaled version
#pSpSN <- createPao(DHSpSN,unitcov = covSpSN,title="Fox Squirrel Spring",unitnames=sitenamesFull)
pSpSN2 <- createPao(DHSpSN,unitcov = covSpSN2,title="Fox Squirrel Spring",unitnames=sitenamesFull)

#Test to see if quadratic effect on height improves the height only model

modsHtest <- list(); 
  modsHtest[[1]] <- occMod(model=list(psi ~ 1, p ~ Height_cm),data = pSpSN2,type = "so")
  modsHtest[[2]] <- occMod(model=list(psi ~ 1, p ~ poly(Height_cm, 2, raw = T)),data = pSpSN2,type = "so")

resultsHtest <- createAicTable(modsHtest)
resultsHtest$table #there is evidence of a quadratic effect on camera height, so it should be in all models with camera height. 

modsHtest[[2]]$beta$p

#Test to see if there is an interaction between EDD and whether there is a log in view

modsEDDtest <- list(); 
modsEDDtest[[1]] <- occMod(model=list(psi ~ 1, p ~ Log.in.View + Squirrel_EDD_WSp),data = pSpSN2,type = "so")
modsEDDtest[[2]] <- occMod(model=list(psi ~ 1, p ~ Log.in.View * Squirrel_EDD_WSp),data = pSpSN2,type = "so")

resultsEDDtest <- createAicTable(modsEDDtest)
resultsEDDtest$table #there is evidence of a interaction effect here, so this needs to be in all models below with both of these parameters

#Full Model Comparison using only a max of 3 predictor variables at a time (though some have more parameters) 
  
mods=list(); i=1
  mods[[i]]=occMod(model=list(psi~1, p~1),    data=pSpSN2,type="so");i=i+1
  mods[[i]]=occMod(model=list(psi~1, p~poly(Height_cm, 2, raw = T)),data=pSpSN2,type="so");i=i+1
  mods[[i]]=occMod(model=list(psi~1, p~Num_Stems),data=pSpSN2,type="so");i=i+1
  mods[[i]]=occMod(model=list(psi~1, p~OakDBH),data=pSpSN2,type="so");i=i+1
  mods[[i]]=occMod(model=list(psi~1, p~Log.in.View),data=pSpSN2,type="so");i=i+1
  mods[[i]]=occMod(model=list(psi~1, p~Squirrel_EDD_WSp),data=pSpSN2,type="so");i=i+1
  mods[[i]]=occMod(model=list(psi~1, p~poly(Height_cm, 2, raw = T) + Num_Stems),data=pSpSN2,type="so");i=i+1
  mods[[i]]=occMod(model=list(psi~1, p~poly(Height_cm, 2, raw = T) + OakDBH),data=pSpSN2,type="so");i=i+1
  mods[[i]]=occMod(model=list(psi~1, p~poly(Height_cm, 2, raw = T) + Log.in.View),data=pSpSN2,type="so");i=i+1
  mods[[i]]=occMod(model=list(psi~1, p~poly(Height_cm, 2, raw = T) + Squirrel_EDD_WSp),data=pSpSN2,type="so");i=i+1
  mods[[i]]=occMod(model=list(psi~1, p~Num_Stems + OakDBH),data=pSpSN2,type="so");i=i+1
  mods[[i]]=occMod(model=list(psi~1, p~Num_Stems + Log.in.View),data=pSpSN2,type="so");i=i+1
  mods[[i]]=occMod(model=list(psi~1, p~Num_Stems + Squirrel_EDD_WSp),data=pSpSN2,type="so");i=i+1
  mods[[i]]=occMod(model=list(psi~1, p~OakDBH + Log.in.View),data=pSpSN2,type="so");i=i+1
  mods[[i]]=occMod(model=list(psi~1, p~OakDBH + Squirrel_EDD_WSp),data=pSpSN2,type="so");i=i+1
  mods[[i]]=occMod(model=list(psi~1, p~Log.in.View * Squirrel_EDD_WSp),data=pSpSN2,type="so");i=i+1
  mods[[i]]=occMod(model=list(psi~1, p~poly(Height_cm, 2, raw = T) + Num_Stems + OakDBH),data=pSpSN2,type="so");i=i+1
  mods[[i]]=occMod(model=list(psi~1, p~poly(Height_cm, 2, raw = T) + Num_Stems + Log.in.View),data=pSpSN2,type="so");i=i+1
  mods[[i]]=occMod(model=list(psi~1, p~poly(Height_cm, 2, raw = T) + Num_Stems + Squirrel_EDD_WSp),data=pSpSN2,type="so");i=i+1
  mods[[i]]=occMod(model=list(psi~1, p~poly(Height_cm, 2, raw = T) + OakDBH + Log.in.View),data=pSpSN2,type="so");i=i+1
  mods[[i]]=occMod(model=list(psi~1, p~poly(Height_cm, 2, raw = T) + OakDBH + Squirrel_EDD_WSp),data=pSpSN2,type="so");i=i+1
  mods[[i]]=occMod(model=list(psi~1, p~poly(Height_cm, 2, raw = T) + Log.in.View * Squirrel_EDD_WSp),data=pSpSN2,type="so");i=i+1
  mods[[i]]=occMod(model=list(psi~1, p~Num_Stems + OakDBH + Log.in.View),data=pSpSN2,type="so");i=i+1
  mods[[i]]=occMod(model=list(psi~1, p~Num_Stems + OakDBH + Squirrel_EDD_WSp),data=pSpSN2,type="so");i=i+1
  mods[[i]]=occMod(model=list(psi~1, p~Num_Stems + Log.in.View * Squirrel_EDD_WSp),data=pSpSN2,type="so");i=i+1
  mods[[i]]=occMod(model=list(psi~1, p~OakDBH + Log.in.View * Squirrel_EDD_WSp),data=pSpSN2,type="so");i=i+1
 
  
#     create AIC table of model results and print
resultsSpSN <- createAicTable(mods)
resultsSpSN$table

bestSpSN <- occMod(model=list(psi~1, p~poly(Height_cm, 2, raw = T)+ Squirrel_EDD_WSp),data=pSpSN2,type="so") 
unique(bestSpSN$real$p) 
unique(bestSpSN$real$psi)
bestSpSN$beta$p

#Average detection probability (Height, Squirrel EDD Wsp)
MeanpFoxSq <- predict(bestSpSN, newdata = data.frame(Height_cm = mean(covSpSN2$Height_cm), Squirrel_EDD_WSp = mean(covSpSN2$Squirrel_EDD_WSp)), param = "p", conf= 0.95)
row.names(MeanpFoxSq) <- "Fox Sq Spring"

#Sum up model weights for each of the 5 covariates
mnames=resultsSpSN$table$Model;
for (s in c("Height","Oak", "Stems", "Log", "EDD")) {
  i <- grep(s,mnames); #  get model names with each covariate name
  wgt=sum(resultsSpSN$table$wgt[i]); #  add up wgts of those models
  cat('CumWgt(',s,')=',wgt,'\n')     #  print sum of weights
}

# Response Curve - Height Fox Sq Spring -----------------------------------

seq.SpHgt <- seq(min(covSpSN2$Height_cm), max(covSpSN2$Height_cm), length.out = 100)

nd.seq.SpHgt <- data.frame(Height_cm = seq.SpHgt, Squirrel_EDD_WSp = median(sqDataSpr$Squirrel_EDD_WSp) )

pred.SpHgt <- predict(bestSpSN, newdata = nd.seq.SpHgt, param = "p", conf= 0.95)

predPlotFSpHgt <- data.frame(nd.seq.SpHgt, pred.SpHgt) 

load(file = "data/responseTheme.Rdata")

FSqSprHgtPlot <- ggplot(predPlotFSpHgt, aes(x=(Height_cm*sd(covSpSN$Height_cm))+mean(covSpSN$Height_cm), y=est)) +
  # Confidence region
  geom_ribbon(aes(ymin=lower_0.95, ymax=upper_0.95), alpha=0.25, show.legend = F) +
  # Prediction Lines
  geom_line() +
  xlab("Camera height (cm)") +   
  ylab("Estimated detection probability") +
  geom_text(aes(x = 34, y = .63, label = "C"), size = 8)+
  myTheme

#ggsave("results/FoxSqSprHgt_DP.tiff", width = 6.5, height = 4.0, units = "in" )

# Fox Squirrel Occupancy Models in Winter ---------------------------------


sitenamesShort <- row.names(as.data.frame(DHWinSN)) 
nsitesShort=nrow(DHWinSN) 
nSURVEYsWinSN=ncol(DHWinSN)  #  set number of sites,surveys from det. history data

#create covariates file

load("data/sqDataWin.RData") #I'm using general squirrel data here because the site covariates are the same for both species
covWinSN <- sqDataWin %>%
  dplyr::select(Deployment_Name, Height_cm, Num_Stems, OakDBH, Log.in.View, Squirrel_EDD_WSp)

#Created scaled version of the covariates file
covWinSN2 <- covWinSN %>%
  mutate(Height_cm = scale(Height_cm),
         Num_Stems = scale(Num_Stems),
         OakDBH = scale(OakDBH),
         Squirrel_EDD_WSp = scale(Squirrel_EDD_WSp))

#Create input file for RPresence
pWinSN <- createPao(DHWinSN,unitcov = covWinSN,title="Fox Squirrel Winter",unitnames=sitenamesShort)
pWinSN2 <- createPao(DHWinSN,unitcov = covWinSN2,title="Fox Squirrel Winter",unitnames=sitenamesShort)


#Test to see if quadratic effect on height improves the height only model

modsHtest <- list(); 
modsHtest[[1]] <- occMod(model=list(psi ~ 1, p ~ Height_cm),data = pWinSN2,type = "so")
modsHtest[[2]] <- occMod(model=list(psi ~ 1, p ~ poly(Height_cm, 2, raw = T)),data = pWinSN2,type = "so")

resultsHtest <- createAicTable(modsHtest)
resultsHtest$table #there is no evidence of a quadratic effect on camera height

#Test to see if there is an interaction between EDD and whether there is a log in view

modsEDDtest <- list(); 
modsEDDtest[[1]] <- occMod(model=list(psi ~ 1, p ~ Log.in.View + Squirrel_EDD_WSp),data = pWinSN2,type = "so")
modsEDDtest[[2]] <- occMod(model=list(psi ~ 1, p ~ Log.in.View * Squirrel_EDD_WSp),data = pWinSN2,type = "so")

resultsEDDtest <- createAicTable(modsEDDtest)
resultsEDDtest$table #there is evidence of a interaction effect here, so this needs to be in all models below with both of these parameters


mods=list(); i=1
mods[[i]]=occMod(model=list(psi~1, p~1),    data=pWinSN2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~Height_cm),data=pWinSN2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~Num_Stems),data=pWinSN2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~OakDBH),data=pWinSN2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~Log.in.View),data=pWinSN2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~Squirrel_EDD_WSp),data=pWinSN2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~Height_cm + Num_Stems),data=pWinSN2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~Height_cm + OakDBH),data=pWinSN2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~Height_cm + Log.in.View),data=pWinSN2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~Height_cm + Squirrel_EDD_WSp),data=pWinSN2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~Num_Stems + OakDBH),data=pWinSN2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~Num_Stems + Log.in.View),data=pWinSN2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~Num_Stems + Squirrel_EDD_WSp),data=pWinSN2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~OakDBH + Log.in.View),data=pWinSN2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~OakDBH + Squirrel_EDD_WSp),data=pWinSN2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~Log.in.View * Squirrel_EDD_WSp),data=pWinSN2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~Height_cm + Num_Stems + OakDBH),data=pWinSN2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~Height_cm + Num_Stems + Log.in.View),data=pWinSN2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~Height_cm + Num_Stems + Squirrel_EDD_WSp),data=pWinSN2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~Height_cm + OakDBH + Log.in.View),data=pWinSN2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~Height_cm + OakDBH + Squirrel_EDD_WSp),data=pWinSN2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~Height_cm + Log.in.View * Squirrel_EDD_WSp),data=pWinSN2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~Num_Stems + OakDBH + Log.in.View),data=pWinSN2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~Num_Stems + OakDBH + Squirrel_EDD_WSp),data=pWinSN2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~Num_Stems + Log.in.View * Squirrel_EDD_WSp),data=pWinSN2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~OakDBH + Log.in.View * Squirrel_EDD_WSp),data=pWinSN2,type="so");i=i+1


#     create AIC table of model results and print
resultsWinSN <- createAicTable(mods)
resultsWinSN$table

mods[[25]]$beta$p
mods[[9]]$beta$p

bestWinSN <- occMod(model=list(psi~1, p ~ Num_Stems + Log.in.View * Squirrel_EDD_WSp),data=pWinSN2,type="so")
bestWinSN$beta$p

MeanpFoxSq[2,] <- predict(bestWinSN, newdata = data.frame(Num_Stems = mean(covWinSN2$Num_Stems), Log.in.View = factor("NO", levels=c("NO", "YES")), Squirrel_EDD_WSp = mean(covWinSN2$Squirrel_EDD_WSp)), param = "p", conf= 0.95)
row.names(MeanpFoxSq)[2] <- "Fox Sq Winter"

#Sum up model weights for each of the 5 covariates
mnames=resultsWinSN$table$Model;
for (s in c("Height","Oak", "Stems", "Log", "EDD")) {
  i <- grep(s,mnames); #  get model names with each covariate name
  wgt=sum(resultsWinSN$table$wgt[i]); #  add up wgts of those models
  cat('CumWgt(',s,')=',wgt,'\n')     #  print sum of weights
}

#Response Curve for EDD log interaction - Fox Sq Winter----------------------------------------------------------------------- 
seq.WEdd <- seq(min(covWinSN2$Squirrel_EDD_WSp), max(covWinSN2$Squirrel_EDD_WSp), length.out = 100) 


nd.seq.WEdd <- data.frame(Squirrel_EDD_WSp = seq.WEdd, Log.in.View = factor("YES", levels=c("NO", "YES")), Num_Stems = median(covWinSN2$Num_Stems))
nd.seq.WEddNo <- data.frame(Squirrel_EDD_WSp = seq.WEdd, Log.in.View = factor("NO", levels=c("NO", "YES")), Num_Stems = median(covWinSN2$Num_Stems))

pred.WEdd <- predict(bestWinSN, newdata = nd.seq.WEdd, param = "p", conf= 0.95)
pred.WEddNo <- predict(bestWinSN, newdata = nd.seq.WEddNo, param = "p", conf= 0.95)


predVals <- rbind(as.data.frame(pred.WEdd), as.data.frame(pred.WEddNo))
allNewVals <- rbind(nd.seq.WEdd, nd.seq.WEddNo)

#I here make a dataframe combining the predictor values, with the estimated response values
predPlot <- data.frame(allNewVals, predVals) 

load(file = "data/responseTheme.Rdata")

FSqWinEDDPlot <- ggplot(predPlot, aes(x=Squirrel_EDD_WSp * sd(covWinSN$Squirrel_EDD_WSp) + mean(covWinSN$Squirrel_EDD_WSp), y=est)) +
  # Confidence region
  geom_ribbon(aes(ymin=lower_0.95, ymax=upper_0.95, fill = Log.in.View), alpha=0.25, show.legend = F) +
  # Prediction Lines
  geom_line(aes(color = Log.in.View), show.legend = F) +
  scale_colour_manual("",values=c("tomato","darkolivegreen"))+
  scale_fill_manual("",values=c("tomato","darkolivegreen"))+
  xlab("Effective Detection Distance (m)") +   
  ylab("Estimated detection probability") +
  geom_text(aes(x = 2.2, y = 1.0, label = "B"), size = 8)+
  myTheme

#ggsave("results/FSqWinEDDLog_DP.tiff", width = 6.5, height = 4.0, units = "in" ) #saves whatever last ggplot was made

# Fox Squirrel Occupancy Models in Summer ---------------------------------


nSURVEYsSumSN=ncol(DHSumSN)  #  set number of sites,surveys from det. history data

#create covariates file

load("data/sqDataSum.RData") #I'm using general squirrel data here because the site covariates are the same for both species
covSumSN <- sqDataSum %>%
  dplyr::select(Deployment_Name, Height_cm, Num_Stems, OakDBH, Log.in.View, Squirrel_EDD_4S)

#Created scaled version of the covariates file
covSumSN2 <- covSumSN %>%
  mutate(Height_cm = scale(Height_cm),
         Num_Stems = scale(Num_Stems),
         OakDBH = scale(OakDBH),
         Squirrel_EDD_4S = scale(Squirrel_EDD_4S))

#Create input file for RPresence as well as version with scaled covariates
pSumSN <- createPao(DHSumSN,unitcov = covSumSN,title="Fox Squirrel Summer",unitnames=sitenamesFull)
pSumSN2 <- createPao(DHSumSN,unitcov = covSumSN2,title="Fox Squirrel Summer - scaled",unitnames=sitenamesFull)

#Test to see if quadratic effect on height improves the height only model

modsHtest <- list(); 
modsHtest[[1]] <- occMod(model=list(psi ~ 1, p ~ Height_cm),data = pSumSN2,type = "so")
modsHtest[[2]] <- occMod(model=list(psi ~ 1, p ~ poly(Height_cm, 2, raw = T)),data = pSumSN2,type = "so")

resultsHtest <- createAicTable(modsHtest)
resultsHtest$table #there is no evidence of a quadratic effect on camera height here so need quadratic in all models for height

#Test to see if there is an interaction between EDD and whether there is a log in view

modsEDDtest <- list(); 
modsEDDtest[[1]] <- occMod(model=list(psi ~ 1, p ~ Log.in.View + Squirrel_EDD_4S),data = pSumSN,type = "so")
modsEDDtest[[2]] <- occMod(model=list(psi ~ 1, p ~ Log.in.View * Squirrel_EDD_4S),data = pSumSN,type = "so")

resultsEDDtest <- createAicTable(modsEDDtest)
resultsEDDtest$table #there is no evidence of a interaction effect here

mods=list(); i=1
mods[[i]]=occMod(model=list(psi~1, p~1),    data=pSumSN2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~Height_cm),data=pSumSN2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~Num_Stems),data=pSumSN2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~OakDBH),data=pSumSN2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~Log.in.View),data=pSumSN2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~Squirrel_EDD_4S),data=pSumSN2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~Height_cm + Num_Stems),data=pSumSN2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~Height_cm + OakDBH),data=pSumSN2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~Height_cm + Log.in.View),data=pSumSN2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~Height_cm + Squirrel_EDD_4S),data=pSumSN2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~Num_Stems + OakDBH),data=pSumSN2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~Num_Stems + Log.in.View),data=pSumSN2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~Num_Stems + Squirrel_EDD_4S),data=pSumSN2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~OakDBH + Log.in.View),data=pSumSN2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~OakDBH + Squirrel_EDD_4S),data=pSumSN2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~Log.in.View + Squirrel_EDD_4S),data=pSumSN2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~Height_cm + Num_Stems + OakDBH),data=pSumSN2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~Height_cm + Num_Stems + Log.in.View),data=pSumSN2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~Height_cm + Num_Stems + Squirrel_EDD_4S),data=pSumSN2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~Height_cm + OakDBH + Log.in.View),data=pSumSN2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~Height_cm + OakDBH + Squirrel_EDD_4S),data=pSumSN2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~Height_cm + Log.in.View + Squirrel_EDD_4S),data=pSumSN2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~Num_Stems + OakDBH + Log.in.View),data=pSumSN2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~Num_Stems + OakDBH + Squirrel_EDD_4S),data=pSumSN2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~Num_Stems + Log.in.View + Squirrel_EDD_4S),data=pSumSN2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~OakDBH + Log.in.View + Squirrel_EDD_4S),data=pSumSN2,type="so");i=i+1


#     create AIC table of model results and print
resultsSumSN <- createAicTable(mods)
resultsSumSN$table

bestSumSN <- occMod(model=list(psi~1, p~Height_cm + Log.in.View + Squirrel_EDD_4S),data=pSumSN2,type="so")
summary(bestSumSN)
unique(bestSumSN$real$p) 
unique(bestSumSN$real$psi)
bestSumSN$beta$p

MeanpFoxSq[3,] <- predict(bestSumSN, newdata = data.frame(Height_cm = mean(covSumSN2$Height_cm), Log.in.View = factor("NO", levels=c("NO", "YES")), Squirrel_EDD_4S = mean(covSumSN2$Squirrel_EDD_4S)), param = "p", conf= 0.95)
row.names(MeanpFoxSq)[3] <- "Fox Sq Summer"

#Sum up model weights for each of the 5 covariates
mnames=resultsSumSN$table$Model;
for (s in c("Height","Oak", "Stems", "Log", "EDD")) {
  i <- grep(s,mnames); #  get model names with each covariate name
  wgt=sum(resultsSumSN$table$wgt[i]); #  add up wgts of those models
  cat('CumWgt(',s,')=',wgt,'\n')     #  print sum of weights
}


# Response Curve - Height Fox Sq Summer -----------------------------------

seq.SHgt <- seq(min(covSumSN2$Height_cm), max(covSumSN2$Height_cm), length.out = 100)

nd.seq.SHgtL <- data.frame(Height_cm = seq.SHgt, Log.in.View = factor("YES", levels=c("NO", "YES")), Squirrel_EDD_4S = median(covSumSN2$Squirrel_EDD_4S) )

pred.SHgtL <- predict(bestSumSN, newdata = nd.seq.SHgtL, param = "p", conf= 0.95)

predPlotFSumHgt <- data.frame(nd.seq.SHgtL, pred.SHgtL) 

load(file = "data/responseTheme.Rdata")

FSqSumHgtPlot <- ggplot(predPlotFSumHgt, aes(x=(Height_cm*sd(covSumSN$Height_cm))+mean(covSumSN$Height_cm), y=est)) +
  # Confidence region
  geom_ribbon(aes(ymin=lower_0.95, ymax=upper_0.95), alpha=0.25, show.legend = F) +
  # Prediction Lines
  geom_line() +
  xlab("Camera Height (cm)") +   
  ylab("Estimated detection probability") +
  geom_text(aes(x = 34, y = .22, label = "B"), size = 8)+
  myTheme

# Fox Squirrel Occupancy Models in Fall ---------------------------------

sitenamesFall <- row.names(as.data.frame(DHFallSN)) 
nsitesFall=nrow(DHFallSN) 
nSURVEYsFallSN=ncol(DHFallSN)  #  set number of sites,surveys from det. history data

#create covariates file

load("data/sqDataFall.RData") #I'm using general squirrel data here because the site covariates are the same for both species
covFallSN <- sqDataFall %>%
  select(Deployment_Name, Height_cm, Num_Stems, OakDBH, Log.in.View, Squirrel_EDD_4S)

#Created scaled version of the covariates file
covFallSN2 <- covFallSN %>%
  mutate(Height_cm = scale(Height_cm),
         Num_Stems = scale(Num_Stems),
         OakDBH = scale(OakDBH),
         Squirrel_EDD_4S = scale(Squirrel_EDD_4S))

#Create input file for RPresence
pFallSN <- createPao(DHFallSN,unitcov = covFallSN,title="Fox Squirrel Fall",unitnames=sitenamesFall)
pFallSN2 <- createPao(DHFallSN,unitcov = covFallSN2,title="Fox Squirrel Fall - scaled",unitnames=sitenamesFall)

#Test to see if quadratic effect on height improves the height only model

modsHtest <- list(); 
modsHtest[[1]] <- occMod(model=list(psi ~ 1, p ~ Height_cm),data = pFallSN2,type = "so")
modsHtest[[2]] <- occMod(model=list(psi ~ 1, p ~ poly(Height_cm, 2, raw = T)),data = pFallSN2,type = "so")

resultsHtest <- createAicTable(modsHtest)
resultsHtest$table #there is no evidence of a quadratic effect on camera height here

#Test to see if there is an interaction between EDD and whether there is a log in view

modsEDDtest <- list(); 
modsEDDtest[[1]] <- occMod(model=list(psi ~ 1, p ~ Log.in.View + Squirrel_EDD_4S),data = pFallSN2,type = "so")
modsEDDtest[[2]] <- occMod(model=list(psi ~ 1, p ~ Log.in.View * Squirrel_EDD_4S),data = pFallSN2,type = "so")

resultsEDDtest <- createAicTable(modsEDDtest)
resultsEDDtest$table #there is no evidence of a interaction effect here


mods=list(); i=1
mods[[i]]=occMod(model=list(psi~1, p~1),    data=pFallSN2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~Height_cm),data=pFallSN2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~Num_Stems),data=pFallSN2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~OakDBH),data=pFallSN2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~Log.in.View),data=pFallSN2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~Squirrel_EDD_4S),data=pFallSN2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~Height_cm + Num_Stems),data=pFallSN2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~Height_cm + OakDBH),data=pFallSN2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~Height_cm + Log.in.View),data=pFallSN2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~Height_cm + Squirrel_EDD_4S),data=pFallSN2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~Num_Stems + OakDBH),data=pFallSN2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~Num_Stems + Log.in.View),data=pFallSN2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~Num_Stems + Squirrel_EDD_4S),data=pFallSN2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~OakDBH + Log.in.View),data=pFallSN2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~OakDBH + Squirrel_EDD_4S),data=pFallSN2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~Log.in.View + Squirrel_EDD_4S),data=pFallSN2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~Height_cm + Num_Stems + OakDBH),data=pFallSN2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~Height_cm + Num_Stems + Log.in.View),data=pFallSN2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~Height_cm + Num_Stems + Squirrel_EDD_4S),data=pFallSN2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~Height_cm + OakDBH + Log.in.View),data=pFallSN2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~Height_cm + OakDBH + Squirrel_EDD_4S),data=pFallSN2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~Height_cm + Log.in.View + Squirrel_EDD_4S),data=pFallSN2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~Num_Stems + OakDBH + Log.in.View),data=pFallSN2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~Num_Stems + OakDBH + Squirrel_EDD_4S),data=pFallSN2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~Num_Stems + Log.in.View + Squirrel_EDD_4S),data=pFallSN2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~OakDBH + Log.in.View + Squirrel_EDD_4S),data=pFallSN2,type="so");i=i+1


#     create AIC table of model results and print
resultsFallSN <- createAicTable(mods)
resultsFallSN$table

bestFallSN <- occMod(model=list(psi~1, p~Num_Stems + Log.in.View + Squirrel_EDD_4S),data=pFallSN2,type="so") 
bestFallSN$beta$p

MeanpFoxSq[4,] <- predict(bestFallSN, newdata = data.frame(Num_Stems = mean(covFallSN2$Num_Stems), Log.in.View = factor("NO", levels=c("NO", "YES")), Squirrel_EDD_4S = mean(covFallSN2$Squirrel_EDD_4S)), param = "p", conf= 0.95)
row.names(MeanpFoxSq)[4] <- "Fox Sq Fall"
save(MeanpFoxSq, file = "results/MeanPSeasons_FoxSq")

#Sum up model weights for each of the 5 covariates
mnames=resultsFallSN$table$Model;
for (s in c("Height","Oak", "Stems", "Log", "EDD")) {
  i <- grep(s,mnames); #  get model names with each covariate name
  wgt=sum(resultsFallSN$table$wgt[i]); #  add up wgts of those models
  cat('CumWgt(',s,')=',wgt,'\n')     #  print sum of weights
}
