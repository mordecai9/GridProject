#Assessment of optimal detection probalitity models for Pooled Squirrel data using occupancy models in Rpresence. Data comes from high resolution 27 camera grid is Posey Hollow. Detection histories were prepared using the script "PrepforDetectionHistory", where R objects of detection hustories were created for each species and each season. These are loaded in here. The script also requires data files prepared for a separate GLM analysis, since these include the site covariates. 


# Clear Environment and Load Packages -------------------------------------

rm(list = ls())
library(tidyverse)
library(RPresence)

# Load Detection Histories Gray Squirrel ------------------------------------------------

load("results/DetHistFall_AllSquirrel")

load("results/DetHistSum_AllSquirrel")

load("results/DetHistWin_AllSquirrel")

load("results/DetHistSpring_AllSquirrel")

# All Squirrel Occupancy Models in Spring ---------------------------------

sitenamesFull <- row.names(as.data.frame(DHSpSQ)) 
nsitesFull=nrow(DHSpSQ) 
nSURVEYsSpSC=ncol(DHSpSQ)  #  set number of sites,surveys from det. history data

#create covariates file

load("data/sqDataSpr.RData")
covSpSQ <- sqDataSpr %>%
  select(Deployment_Name, Height_cm, Num_Stems, OakDBH, Log.in.View, Squirrel_EDD_WSp)

#Created scaled version of the covariates file
covSpSQ2 <- covSpSQ %>%
  mutate(Height_cm = scale(Height_cm),
         Num_Stems = scale(Num_Stems),
         OakDBH = scale(OakDBH),
         Squirrel_EDD_WSp = scale(Squirrel_EDD_WSp))

#Create input file for RPresence
pSpSQ <- createPao(DHSpSQ,unitcov = covSpSQ,title="All Squirrel Spring",unitnames=sitenamesFull)
pSpSQ2 <- createPao(DHSpSQ,unitcov = covSpSQ2,title="All Squirrel Spring-scaled",unitnames=sitenamesFull)

#Test to see if quadratic effect on height improves the height only model

modsHtest <- list(); 
modsHtest[[1]] <- occMod(model=list(psi ~ 1, p ~ Height_cm),data = pSpSQ2,type = "so")
modsHtest[[2]] <- occMod(model=list(psi ~ 1, p ~ poly(Height_cm, 2, raw = T)),data = pSpSQ2,type = "so")

resultsHtest <- createAicTable(modsHtest)
resultsHtest$table #there is evidence of a quadratic effect on camera height, so it should be in all models with camera height

#Test to see if there is an interaction between EDD and whether there is a log in view

modsEDDtest <- list(); 
modsEDDtest[[1]] <- occMod(model=list(psi ~ 1, p ~ Log.in.View + Squirrel_EDD_WSp),data = pSpSQ2,type = "so")
modsEDDtest[[2]] <- occMod(model=list(psi ~ 1, p ~ Log.in.View * Squirrel_EDD_WSp),data = pSpSQ2,type = "so")

resultsEDDtest <- createAicTable(modsEDDtest)
resultsEDDtest$table #there is evidence of an interaction effect here, so this needs to be in all models below with both of these parameters


mods=list(); i=1
mods[[i]]=occMod(model=list(psi~1, p~1),    data=pSpSQ2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~poly(Height_cm, 2, raw = T)),data=pSpSQ2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~Num_Stems),data=pSpSQ2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~OakDBH),data=pSpSQ2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~Log.in.View),data=pSpSQ2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~Squirrel_EDD_WSp),data=pSpSQ2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~poly(Height_cm, 2, raw = T) + Num_Stems),data=pSpSQ2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~poly(Height_cm, 2, raw = T) + OakDBH),data=pSpSQ2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~poly(Height_cm, 2, raw = T) + Log.in.View),data=pSpSQ2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~poly(Height_cm, 2, raw = T) + Squirrel_EDD_WSp),data=pSpSQ2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~Num_Stems + OakDBH),data=pSpSQ2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~Num_Stems + Log.in.View),data=pSpSQ2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~Num_Stems + Squirrel_EDD_WSp),data=pSpSQ2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~OakDBH + Log.in.View),data=pSpSQ2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~OakDBH + Squirrel_EDD_WSp),data=pSpSQ2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~Log.in.View * Squirrel_EDD_WSp),data=pSpSQ2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~poly(Height_cm, 2, raw = T) + Num_Stems + OakDBH),data=pSpSQ2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~poly(Height_cm, 2, raw = T) + Num_Stems + Log.in.View),data=pSpSQ2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~poly(Height_cm, 2, raw = T) + Num_Stems + Squirrel_EDD_WSp),data=pSpSQ2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~poly(Height_cm, 2, raw = T) + OakDBH + Log.in.View),data=pSpSQ2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~poly(Height_cm, 2, raw = T) + OakDBH + Squirrel_EDD_WSp),data=pSpSQ2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~poly(Height_cm, 2, raw = T) + Log.in.View * Squirrel_EDD_WSp),data=pSpSQ2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~Num_Stems + OakDBH + Log.in.View),data=pSpSQ2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~Num_Stems + OakDBH + Squirrel_EDD_WSp),data=pSpSQ2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~Num_Stems + Log.in.View * Squirrel_EDD_WSp),data=pSpSQ2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~OakDBH + Log.in.View * Squirrel_EDD_WSp),data=pSpSQ2,type="so");i=i+1


#     create AIC table of model results and print
resultsSpSQ <- createAicTable(mods)
resultsSpSQ$table

bestSpSQ <- occMod(model=list(psi~1, p~Log.in.View * Squirrel_EDD_WSp),data=pSpSQ2,type="so") 
bestSpSQ$beta$p

#Average detection probability (Log, Squirrel EDD Wsp)
MeanpSq <- predict(bestSpSQ, newdata = data.frame(Log.in.View = factor("NO", levels=c("NO", "YES")), Squirrel_EDD_WSp = mean(covSpSQ2$Squirrel_EDD_WSp)), param = "p", conf= 0.95)
row.names(MeanpSq) <- "All Sq Spring"

#Sum up model weights for each of the 5 covariates
mnames=resultsSpSQ$table$Model;
for (s in c("Height","Oak", "Stems", "Log", "EDD")) {
  i <- grep(s,mnames); #  get model names with each covariate name
  wgt=sum(resultsSpSQ$table$wgt[i]); #  add up wgts of those models
  cat('CumWgt(',s,')=',wgt,'\n')     #  print sum of weights
}


# All Squirrel Occupancy Models in Winter ---------------------------------

sitenamesShort <- row.names(as.data.frame(DHWinSQ)) 
nsitesShort=nrow(DHWinSQ) 
nSURVEYsWinSC=ncol(DHWinSQ)  #  set number of sites,surveys from det. history data

#create covariates file

load("data/sqDataWin.RData") #I'm using general squirrel data here because the site covariates are the same for both species
covWinSQ <- sqDataWin %>%
  dplyr::select(Deployment_Name, Height_cm, Num_Stems, OakDBH, Log.in.View, Squirrel_EDD_WSp)

#Created scaled version of the covariates file
covWinSQ2 <- covWinSQ %>%
  mutate(Height_cm = scale(Height_cm),
         Num_Stems = scale(Num_Stems),
         OakDBH = scale(OakDBH),
         Squirrel_EDD_WSp = scale(Squirrel_EDD_WSp))


#Create input file for RPresence
pWinSQ <- createPao(DHWinSQ,unitcov = covWinSQ,title="All Squirrel Winter",unitnames=sitenamesShort)
pWinSQ2 <- createPao(DHWinSQ,unitcov = covWinSQ2,title="AllSquirrelWinter-scaled",unitnames=sitenamesShort)

#Test to see if quadratic effect on height improves the height only model

modsHtest <- list(); 
modsHtest[[1]] <- occMod(model=list(psi ~ 1, p ~ Height_cm),data = pWinSQ2,type = "so")
modsHtest[[2]] <- occMod(model=list(psi ~ 1, p ~ poly(Height_cm, 2, raw = T)),data = pWinSQ2,type = "so")

resultsHtest <- createAicTable(modsHtest)
resultsHtest$table #there is evidence of a quadratic effect on camera height, so this needs to be in models below
modsHtest[[2]]$beta$p

#Test to see if there is an interaction between EDD and whether there is a log in view

modsEDDtest <- list(); 
modsEDDtest[[1]] <- occMod(model=list(psi ~ 1, p ~ Log.in.View + Squirrel_EDD_WSp),data = pWinSQ2,type = "so")
modsEDDtest[[2]] <- occMod(model=list(psi ~ 1, p ~ Log.in.View * Squirrel_EDD_WSp),data = pWinSQ2,type = "so")

resultsEDDtest <- createAicTable(modsEDDtest)
resultsEDDtest$table #there is evidence of a interaction effect here, so this needs to be in all models below with both of these parameters


mods=list(); i=1
mods[[i]]=occMod(model=list(psi~1, p~1),    data=pWinSQ2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~poly(Height_cm, 2, raw = T)),data=pWinSQ2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~Num_Stems),data=pWinSQ2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~OakDBH),data=pWinSQ2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~Log.in.View),data=pWinSQ2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~Squirrel_EDD_WSp),data=pWinSQ2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~poly(Height_cm, 2, raw = T) + Num_Stems),data=pWinSQ2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~poly(Height_cm, 2, raw = T) + OakDBH),data=pWinSQ2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~poly(Height_cm, 2, raw = T) + Log.in.View),data=pWinSQ2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~poly(Height_cm, 2, raw = T) + Squirrel_EDD_WSp),data=pWinSQ2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~Num_Stems + OakDBH),data=pWinSQ2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~Num_Stems + Log.in.View),data=pWinSQ2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~Num_Stems + Squirrel_EDD_WSp),data=pWinSQ2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~OakDBH + Log.in.View),data=pWinSQ2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~OakDBH + Squirrel_EDD_WSp),data=pWinSQ2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~Log.in.View * Squirrel_EDD_WSp),data=pWinSQ2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~poly(Height_cm, 2, raw = T) + Num_Stems + OakDBH),data=pWinSQ2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~poly(Height_cm, 2, raw = T) + Num_Stems + Log.in.View),data=pWinSQ2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~poly(Height_cm, 2, raw = T) + Num_Stems + Squirrel_EDD_WSp),data=pWinSQ2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~poly(Height_cm, 2, raw = T) + OakDBH + Log.in.View),data=pWinSQ2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~poly(Height_cm, 2, raw = T) + OakDBH + Squirrel_EDD_WSp),data=pWinSQ2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~poly(Height_cm, 2, raw = T) + Log.in.View * Squirrel_EDD_WSp),data=pWinSQ2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~Num_Stems + OakDBH + Log.in.View),data=pWinSQ2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~Num_Stems + OakDBH + Squirrel_EDD_WSp),data=pWinSQ2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~Num_Stems + Log.in.View * Squirrel_EDD_WSp),data=pWinSQ2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~OakDBH + Log.in.View * Squirrel_EDD_WSp),data=pWinSQ2,type="so");i=i+1


#     create AIC table of model results and print
resultsWinSQ <- createAicTable(mods)
resultsWinSQ$table

bestWinSQ <- occMod(model=list(psi~1, p ~ Num_Stems + Log.in.View * Squirrel_EDD_WSp),data=pWinSQ2,type="so")
bestWinSQ$beta$p

#Average detection probability (Stems, Log, Squirrel EDD Wsp)
MeanpSq[2,] <- predict(bestWinSQ, newdata = data.frame(Num_Stems = mean(covWinSQ2$Num_Stems), Log.in.View = factor("NO", levels=c("NO", "YES")), Squirrel_EDD_WSp = mean(covWinSQ2$Squirrel_EDD_WSp)), param = "p", conf= 0.95)
row.names(MeanpSq)[2] <- "All Sq Winter"

#Sum up model weights for each of the 5 covariates
mnames=resultsWinSQ$table$Model;
for (s in c("Height","Oak", "Stems", "Log", "EDD")) {
  i <- grep(s,mnames); #  get model names with each covariate name
  wgt=sum(resultsWinSQ$table$wgt[i]); #  add up wgts of those models
  cat('CumWgt(',s,')=',wgt,'\n')     #  print sum of weights
}

#Response Curve for EDD log interaction - All Sq Winter----------------------------------------------------------------------- 
seq.WEdd <- seq(min(covWinSQ2$Squirrel_EDD_WSp), max(covWinSQ2$Squirrel_EDD_WSp), length.out = 100) 


nd.seq.WEdd <- data.frame(Squirrel_EDD_WSp = seq.WEdd, Log.in.View = factor("YES", levels=c("NO", "YES")), Num_Stems = median(covWinSQ2$Num_Stems))
nd.seq.WEddNo <- data.frame(Squirrel_EDD_WSp = seq.WEdd, Log.in.View = factor("NO", levels=c("NO", "YES")), Num_Stems = median(covWinSQ2$Num_Stems))

pred.WEdd <- predict(bestWinSQ, newdata = nd.seq.WEdd, param = "p", conf= 0.95)
pred.WEddNo <- predict(bestWinSQ, newdata = nd.seq.WEddNo, param = "p", conf= 0.95)


predVals <- rbind(as.data.frame(pred.WEdd), as.data.frame(pred.WEddNo))
allNewVals <- rbind(nd.seq.WEdd, nd.seq.WEddNo)

#I here make a dataframe combining the predictor values, with the estimated response values
predPlot <- data.frame(allNewVals, predVals) 

load(file = "data/responseTheme.Rdata")

SqWinEDDPlot <- ggplot(predPlot, aes(x=Squirrel_EDD_WSp * sd(covWinSQ$Squirrel_EDD_WSp) + mean(covWinSQ$Squirrel_EDD_WSp), y=est)) +
  # Confidence region
  geom_ribbon(aes(ymin=lower_0.95, ymax=upper_0.95, fill = Log.in.View), alpha=0.25, show.legend = F) +
  # Prediction Lines
  geom_line(aes(color = Log.in.View), show.legend = F) +
  scale_colour_manual("",values=c("tomato","darkolivegreen"))+
  scale_fill_manual("",values=c("tomato","darkolivegreen"))+
  xlab("Effective Detection Distance (m)") +   
  ylab("Est. det. probability") +
  myTheme
  
ggsave("results/SqWinEDDLog_DP_S2.tiff", width = 6.5, height = 4.0, units = "in", dpi = 300 ) #saves whatever last ggplot was made

# All Squirrel Occupancy Models in Summer ---------------------------------


nSURVEYsSumSC=ncol(DHSumSQ)  #  set number of sites,surveys from det. history data

#create covariates file

load("data/sqDataSum.RData") #I'm using general squirrel data here because the site covariates are the same for both species
covSumSQ <- sqDataSum %>%
  select(Deployment_Name, Height_cm, Num_Stems, OakDBH, Log.in.View, Squirrel_EDD_4S)

covSumSQ2 <- covSumSQ %>%
  mutate(Height_cm = scale(Height_cm),
         Num_Stems = scale(Num_Stems),
         OakDBH = scale(OakDBH),
         SSquirrel_EDD_4S = scale(Squirrel_EDD_4S))

#Create input file for RPresence
pSumSQ <- createPao(DHSumSQ,unitcov = covSumSQ,title="All Squirrel Summer",unitnames=sitenamesFull)
pSumSQ2 <- createPao(DHSumSQ,unitcov = covSumSQ2,title="AllSquirrelSummer-scaled",unitnames=sitenamesFull)

#Test to see if quadratic effect on height improves the height only model

modsHtest <- list(); 
modsHtest[[1]] <- occMod(model=list(psi ~ 1, p ~ Height_cm),data = pSumSQ2,type = "so")
modsHtest[[2]] <- occMod(model=list(psi ~ 1, p ~ poly(Height_cm, 2, raw = T)),data = pSumSQ2,type = "so")

resultsHtest <- createAicTable(modsHtest)
resultsHtest$table #there is no evidence of a quadratic effect on camera height 

#Test to see if there is an interaction between EDD and whether there is a log in view

modsEDDtest <- list(); 
modsEDDtest[[1]] <- occMod(model=list(psi ~ 1, p ~ Log.in.View + Squirrel_EDD_4S),data = pSumSQ2,type = "so")
modsEDDtest[[2]] <- occMod(model=list(psi ~ 1, p ~ Log.in.View * Squirrel_EDD_4S),data = pSumSQ2,type = "so")

resultsEDDtest <- createAicTable(modsEDDtest)
resultsEDDtest$table #there is no evidence of a interaction effect here

mods=list(); i=1
mods[[i]]=occMod(model=list(psi~1, p~1),    data=pSumSQ2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~Height_cm),data=pSumSQ2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~Num_Stems),data=pSumSQ2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~OakDBH),data=pSumSQ2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~Log.in.View),data=pSumSQ2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~Squirrel_EDD_4S),data=pSumSQ2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~Height_cm + Num_Stems),data=pSumSQ2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~Height_cm + OakDBH),data=pSumSQ2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~Height_cm + Log.in.View),data=pSumSQ2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~Height_cm + Squirrel_EDD_4S),data=pSumSQ2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~Num_Stems + OakDBH),data=pSumSQ2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~Num_Stems + Log.in.View),data=pSumSQ2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~Num_Stems + Squirrel_EDD_4S),data=pSumSQ2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~OakDBH + Log.in.View),data=pSumSQ2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~OakDBH + Squirrel_EDD_4S),data=pSumSQ2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~Log.in.View + Squirrel_EDD_4S),data=pSumSQ2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~Height_cm + Num_Stems + OakDBH),data=pSumSQ2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~Height_cm + Num_Stems + Log.in.View),data=pSumSQ2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~Height_cm + Num_Stems + Squirrel_EDD_4S),data=pSumSQ2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~Height_cm + OakDBH + Log.in.View),data=pSumSQ2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~Height_cm + OakDBH + Squirrel_EDD_4S),data=pSumSQ2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~Height_cm + Log.in.View + Squirrel_EDD_4S),data=pSumSQ2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~Num_Stems + OakDBH + Log.in.View),data=pSumSQ2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~Num_Stems + OakDBH + Squirrel_EDD_4S),data=pSumSQ2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~Num_Stems + Log.in.View + Squirrel_EDD_4S),data=pSumSQ2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~OakDBH + Log.in.View + Squirrel_EDD_4S),data=pSumSQ2,type="so");i=i+1


#     create AIC table of model results and print
resultsSumSQ <- createAicTable(mods)
resultsSumSQ$table

bestSumSQ <- occMod(model = list(psi~1, p ~ Height_cm + Log.in.View + Squirrel_EDD_4S),data = pSumSQ2,type = "so")
bestSumSQ$beta$p

#Average detection probability (Height, Log, Squirrel EDD 4S)
MeanpSq[3,] <- predict(bestSumSQ, newdata = data.frame(Height_cm = mean(covSumSQ2$Height_cm), Log.in.View = factor("NO", levels=c("NO", "YES")), Squirrel_EDD_4S = mean(covSumSQ2$Squirrel_EDD_4S)), param = "p", conf= 0.95)
row.names(MeanpSq)[3] <- "All Sq Summer"

#Sum up model weights for each of the 5 covariates
mnames=resultsSumSQ$table$Model;
for (s in c("Height","Oak", "Stems", "Log", "EDD")) {
  i <- grep(s,mnames); #  get model names with each covariate name
  wgt=sum(resultsSumSQ$table$wgt[i]); #  add up wgts of those models
  cat('CumWgt(',s,')=',wgt,'\n')     #  print sum of weights
}

# All Squirrel Occupancy Models in Fall ---------------------------------

nSURVEYsFallSC=ncol(DHFallSQ)  #  set number of sites,surveys from det. history data

#create covariates file

load("data/sqDataFall.RData") #I'm using general squirrel data here because the site covariates are the same for both species
covFallSQ <- sqDataFall %>%
  select(Deployment_Name, Height_cm, Num_Stems, OakDBH, Log.in.View, Squirrel_EDD_4S)

#Scaled version of site covariates
covFallSQ2 <- covFallSQ %>%
  mutate(Height_cm = scale(Height_cm),
         Num_Stems = scale(Num_Stems),
         OakDBH = scale(OakDBH),
         SSquirrel_EDD_4S = scale(Squirrel_EDD_4S))

#Create input file for RPresence
pFallSQ <- createPao(DHFallSQ,unitcov = covFallSQ,title="All Squirrel Fall",unitnames=sitenamesShort)
pFallSQ2 <- createPao(DHFallSQ,unitcov = covFallSQ2,title="AllSquirrelFall-scaled",unitnames=sitenamesShort)

#Test to see if quadratic effect on height improves the height only model

modsHtest <- list(); 
modsHtest[[1]] <- occMod(model=list(psi ~ 1, p ~ Height_cm),data = pFallSQ2,type = "so")
modsHtest[[2]] <- occMod(model=list(psi ~ 1, p ~ poly(Height_cm, 2, raw = T)),data = pFallSQ2,type = "so")

resultsHtest <- createAicTable(modsHtest)
resultsHtest$table #there is evidence of a quadratic effect on camera height here, so it needs to be in all models below

#Test to see if there is an interaction between EDD and whether there is a log in view

modsEDDtest <- list(); 
modsEDDtest[[1]] <- occMod(model=list(psi ~ 1, p ~ Log.in.View + Squirrel_EDD_4S),data = pFallSQ2,type = "so")
modsEDDtest[[2]] <- occMod(model=list(psi ~ 1, p ~ Log.in.View * Squirrel_EDD_4S),data = pFallSQ2,type = "so")

resultsEDDtest <- createAicTable(modsEDDtest)
resultsEDDtest$table #there is evidence of a interaction effect here so it needs to be reflected in models below


mods=list(); i=1
mods[[i]]=occMod(model=list(psi~1, p~1),    data=pFallSQ2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~poly(Height_cm, 2, raw = T)),data=pFallSQ2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~Num_Stems),data=pFallSQ2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~OakDBH),data=pFallSQ2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~Log.in.View),data=pFallSQ2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~Squirrel_EDD_4S),data=pFallSQ2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~poly(Height_cm, 2, raw = T) + Num_Stems),data=pFallSQ2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~poly(Height_cm, 2, raw = T) + OakDBH),data=pFallSQ2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~poly(Height_cm, 2, raw = T) + Log.in.View),data=pFallSQ2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~poly(Height_cm, 2, raw = T) + Squirrel_EDD_4S),data=pFallSQ2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~Num_Stems + OakDBH),data=pFallSQ2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~Num_Stems + Log.in.View),data=pFallSQ2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~Num_Stems + Squirrel_EDD_4S),data=pFallSQ2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~OakDBH + Log.in.View),data=pFallSQ2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~OakDBH + Squirrel_EDD_4S),data=pFallSQ2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~Log.in.View * Squirrel_EDD_4S),data=pFallSQ2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~poly(Height_cm, 2, raw = T) + Num_Stems + OakDBH),data=pFallSQ2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~poly(Height_cm, 2, raw = T) + Num_Stems + Log.in.View),data=pFallSQ2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~poly(Height_cm, 2, raw = T) + Num_Stems + Squirrel_EDD_4S),data=pFallSQ2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~poly(Height_cm, 2, raw = T) + OakDBH + Log.in.View),data=pFallSQ2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~poly(Height_cm, 2, raw = T) + OakDBH + Squirrel_EDD_4S),data=pFallSQ2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~poly(Height_cm, 2, raw = T) + Log.in.View * Squirrel_EDD_4S),data=pFallSQ2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~Num_Stems + OakDBH + Log.in.View),data=pFallSQ2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~Num_Stems + OakDBH + Squirrel_EDD_4S),data=pFallSQ2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~Num_Stems + Log.in.View * Squirrel_EDD_4S),data=pFallSQ2,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1, p~OakDBH + Log.in.View * Squirrel_EDD_4S),data=pFallSQ2,type="so");i=i+1


#     create AIC table of model results and print
resultsFallSQ <- createAicTable(mods)
resultsFallSQ$table

bestFallSQ <- occMod(model=list(psi~1, p ~ Log.in.View * Squirrel_EDD_4S),data=pFallSQ2,type="so") 
bestFallSQ$beta$p

#Average detection probability (Log, Squirrel EDD 4S)
MeanpSq[4,] <- predict(bestFallSQ, newdata = data.frame(Log.in.View = factor("NO", levels=c("NO", "YES")), Squirrel_EDD_4S = mean(covFallSQ2$Squirrel_EDD_4S)), param = "p", conf= 0.95)
row.names(MeanpSq)[4] <- "All Sq Fall"
save(MeanpSq, file = "results/MeanPSeasons_AllSq")

#Sum up model weights for each of the 5 covariates
mnames=resultsFallSQ$table$Model;
for (s in c("Height","Oak", "Stems", "Log", "EDD")) {
  i <- grep(s,mnames); #  get model names with each covariate name
  wgt=sum(resultsFallSQ$table$wgt[i]); #  add up wgts of those models
  cat('CumWgt(',s,')=',wgt,'\n')     #  print sum of weights
}