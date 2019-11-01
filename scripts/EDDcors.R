#Investigating whether cameras with higher detection distances for one species group have higher detections distances for other size groups.

rm(list = ls())

source("scripts/pairsPannelFunctions.r")
library(tidyverse)

#We do not have unique detection distance data for each species in each season, although that would be ideal. We had to pool across seasons for most, and in some cases use data from other species. Details are below for each species.

#Deer EDD Data. For white-tailed deer, it was necessary to pool Summer/Fall data to get EDD information for Summer and Fall. Here we had a minimum # of observations of 31 for each camera. We also pooled Winter and Spring data for White-tailed deer, and here unfortunately had 5 cameras with less than 30 observations (min 18)

Deer_EDD<-read.csv("data/DeerEDD_4S.csv")


#Squirrel EDD Data. For squirrels, we had adequate observations to estimate detection distance for squirrels in Winter/Spring when these seasons were pooled, but had to use pooled 4-season data for Summer and Fall analyses for squirrels. 


sqEDD <- read.csv("data/SquirrelEDD_4S.csv") 

#Raccoon EDD Data. For raccoons, when pooling across all 4 seasons, approximately half the cameras failed to meet the 30 observation threshold. For those cameras with adequate observations we calculated EDD. For the remaining, we used the EDDs from DEER, which were determined to correlate better with raccoon EDD than squirrel. 

racEDD <- read.csv("data/RaccoonEDD_4S.csv") 




boxplot(Deer_EDD[, 2:3]) #So this makes sense, overall the cameras can see further when the vegetation is not there. Though I though it would be a bigger difference
plot(Summer.Fall.EDD ~ Winter.Spring.EDD, data = Deer_EDD) # You would think this would show a much stronger correlation. The cameras with long detection distances in Summer/Fall should have longer detection distances in Winter/Spring. Generally it is true but the relationship is not strong.
lm1 <- lm(Summer.Fall.EDD ~ Winter.Spring.EDD, data = Deer_EDD)
summary(lm1)
abline(lm1)


# Combine all EDD data and look at correlations ---------------------------

eddAll <- Deer_EDD %>%
  rename(SumFallDeer = Summer.Fall.EDD) %>%
  rename(WingSprDeer =  Winter.Spring.EDD) %>%
  left_join(sqEDD, by = "Deployment") %>%
  left_join(racEDD, by = "Deployment")

pairs(eddAll[,2:6], diag.panel = panel.hist, lower.panel = panel.smooth, upper.panel = panel.cor)
