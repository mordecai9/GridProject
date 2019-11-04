#Jaccard Similarity Index for Species Detection History
#Uses detection histories created in "PrepforDetectionHistory.R" script. Trying to find out if detection histories at cameras near each other are more similar than those further apart, just as we can do with Moran's I and correlograms for capture rates.

#Galvez et al 2016 did this, but we do not have access to his script.

library(jaccard) #this doesn't seem to work because it has "qvalue" as a dependency and this is no longer available it seems

##Load species files for Spring
load("results/DetHistSpSciurus niger")
DHSpSN <- temp

load("results/DetHistSpSciurus carolinensis")
DHSpSC <- temp

load("results/DetHistSpOdocoileus virginianus")
DHSpOV <- temp

load("results/DetHistSpUrsus americanus")
DHSpUA <- temp

load("results/DetHistSpProcyon lotor")
DHSpPL <- temp

##Load species files for Summer
load("results/DetHistSumSciurus niger")
DHSumSN <- temp

load("results/DetHistSumSciurus carolinensis")
DHSumSC <- temp

load("results/DetHistSumOdocoileus virginianus")
DHSumOV <- temp

load("results/DetHistSumUrsus americanus")
DHSumUA <- temp

load("results/DetHistSumProcyon lotor")
DHSumPL <- temp

##Load species for Fall
load("results/DetHistFallSciurus niger")
DHFallSN <- temp

load("results/DetHistFallSciurus carolinensis")
DHFallSC <- temp

load("results/DetHistFallOdocoileus virginianus")
DHFallOV <- temp

load("results/DetHistFallUrsus americanus")
DHFallUA <- temp

load("results/DetHistFallProcyon lotor")
DHFallPL <- temp

##Load species for Winter
load("results/DetHistWinSciurus niger")
DHWinSN <- temp

load("results/DetHistWinSciurus carolinensis")
DHWinSC <- temp

load("results/DetHistWinOdocoileus virginianus")
DHWinOV <- temp

load("results/DetHistWinUrsus americanus")
DHWinUA <- temp

load("results/DetHistWinProcyon lotor")
DHWinPL <- temp

#Changing Temp Files to Data Frames With the Camera names as column names
DFDHFALLOV<-as.data.frame(t(DHFallOV))

##Jaccard Index between Cameras (same season and species)##

#First code attempt
#If my math is correct, this produces a wrong answer
library('clusteval')
cluster_similarity(DFDHFALLOV$`0101_F17` , DFDHFALLOV$`0102_F17`, similarity="jaccard", method="independence")
#Looks like this cluster similarity ignores position, so it doesn't matter when the 1s and 0s happen, I think. So we shouldn't use this.


# Creating Jaccard Index Function -----------------------------------------
#Taking jaccard function from Stack Overflow
##I believe "margin ==2" means columns and "margin ==1" is rows.
##The original code is from https://stats.stackexchange.com/questions/176613/jaccard-similarity-in-r

#Using the modification lower down in the thread that allows for all row combinations, which is what we need, when traps are as rows

# Jaccard Index
library(tidyverse)

# Test dataset
#df <- data.frame(t(data.frame(c1=rnorm(100),
#                              c2=rnorm(100),
#                              c3=rnorm(100),
#                              c4=rnorm(100),
#                              c5=rnorm(100),
#                              c6=rnorm(100))))

#df[df > 0] <- 1
#df[df <= 0] <- 0
#df

# Function returns the Jaccard index and Jaccard distance
# Parameters:
# 1. df, dataframe of interest
# 2. margin, axis in which the apply function is meant to move along
jaccard <- function(df, margin=1) {
  if (margin == 1 | margin == 2) {
    M_00 <- apply(df, margin, sum) == 0
    M_11 <- apply(df, margin, sum) == 2
    if (margin == 1) {
      df <- df[!M_00, ]
      JSim <- sum(M_11) / nrow(df)
    } else {
      df <- df[, !M_00]
      JSim <- sum(M_11) / length(df)
    }
    JDist <- 1 - JSim
    return(c(JSim = JSim, JDist = JDist))
  } else break
}

jaccard(df[1:2,], margin=2)


jaccard_per_row <- function(df, margin=1){
  require(magrittr)
  require(dplyr)
  key_pairs <- expand.grid(row.names(df), row.names(df))
  results <- t(apply(key_pairs, 1, function(row) jaccard(df[c(row[1], row[2]),], margin=margin)))
  key_pair <- key_pairs %>% mutate(pair = paste(Var1,"_",Var2,sep=""))
  results <- data.frame(results)
  row.names(results) <- key_pair$pair
  results
}

#jaccard_per_row(df, margin=2)


# Import Distance Matrix for All Cameras ----------------------------------
#Import distance matrix for all camera combinations. This was made in the script "GridAutoCorrelation.R"
load("data/distMatrix")



# Fox Squirrel Jaccard Calculations and Graph -----------------------------
#Need to remove cameras and ocassions that had NAs to make it work. First I pull out the first ocassion, because cameras were set on different days. Then I remove rows that have NAs.

#Spring Fox Squirrel
DHSpSN_b <- as.data.frame(DHSpSN) %>%
  dplyr::select(-1) 

DHSpSN_b <- DHSpSN_b[which(complete.cases(DHSpSN_b)),]

#Calculates jaccard similarity, and the inverse   
jacSpSN <- jaccard_per_row(DHSpSN_b, margin=2)

#Creating columns representing each of the two cameras
jacSpSN$cam1 <- gsub("_.*", "", row.names(jacSpSN))
jacSpSN$cam2 <- substr(row.names(jacSpSN), 11, 14)

#removes comparisons of same camera to same camera
jacSpSN <-jacSpSN[which(jacSpSN$cam1 != jacSpSN$cam2), ]

#Pull intercamera distances from distance matrix based on value of the two cameras
jacSpSN$dist <- NA

for (i in 1:nrow(jacSpSN)) {
  jacSpSN$dist[i] <- distMatS_named[as.character(as.numeric(jacSpSN$cam1[i])),as.character(as.numeric(jacSpSN$cam2[i]))]

}

#Graph the similarity index by distance

foxSqJacPlot_SP <- ggplot(data = jacSpSN, mapping = aes(x = dist, y = JSim)) +
  geom_point() +
  ylim(0, 1) +
  xlim(10,120)
foxSqJacPlot_SP

#Summer Fox Squirrel
DHSumSN_b <- as.data.frame(DHSumSN) %>%
  dplyr::select(-1) 

DHSumSN_b <- DHSumSN_b[which(complete.cases(DHSumSN_b)),]

#Calculates jaccard similarity, and the inverse   
jacSumSN <- jaccard_per_row(DHSumSN_b, margin=2)

#Creating columns representing each of the two cameras
jacSumSN$cam1 <- gsub("-.*", "", row.names(jacSumSN))
jacSumSN$cam2 <- substr(row.names(jacSumSN), 10, 13)

#removes comparisons of same camera to same camera
jacSumSN <-jacSumSN[which(jacSumSN$cam1 != jacSumSN$cam2), ]

#Pull intercamera distances from distance matrix based on value of the two cameras
jacSumSN$dist <- NA

for (i in 1:nrow(jacSumSN)) {
  jacSumSN$dist[i] <- distMatS_named[as.character(as.numeric(jacSumSN$cam1[i])),as.character(as.numeric(jacSumSN$cam2[i]))]
  
}

#Graph the similarity index by distance

foxSqJacPlot_SUM <- ggplot(data = jacSumSN, mapping = aes(x = dist, y = JSim)) +
  geom_point() +
  ylim(0, 1) +
  xlim(10,120)
foxSqJacPlot_SUM

#Winter Fox Squirrel
DHWinSN_b <- as.data.frame(DHWinSN) %>%
  dplyr::select(-1) 

DHWinSN_b <- DHWinSN_b[which(complete.cases(DHWinSN_b)),]

#Calculates jaccard similarity, and the inverse   
jacWinSN <- jaccard_per_row(DHWinSN_b, margin=2)

#Creating columns representing each of the two cameras
jacWinSN$cam1 <- gsub("_.*", "", row.names(jacWinSN))
jacWinSN$cam2 <- substr(row.names(jacWinSN), 10, 13)

#removes comparisons of same camera to same camera
jacWinSN <-jacWinSN[which(jacWinSN$cam1 != jacWinSN$cam2), ]

#Pull intercamera distances from distance matrix based on value of the two cameras
jacWinSN$dist <- NA

for (i in 1:nrow(jacWinSN)) {
  jacWinSN$dist[i] <- distMatS_named[as.character(as.numeric(jacWinSN$cam1[i])),as.character(as.numeric(jacWinSN$cam2[i]))]
  
}

#Graph the similarity index by distance

foxSqJacPlot_Win <- ggplot(data = jacWinSN, mapping = aes(x = dist, y = JSim)) +
  geom_point() +
  ylim(0, 1) +
  xlim(10,120)
foxSqJacPlot_Win

#Winter Fox Squirrel
DHFallSN_b <- as.data.frame(DHFallSN) %>%
  dplyr::select(-1) 

DHFallSN_b <- DHFallSN_b[which(complete.cases(DHFallSN_b)),]

#Calculates jaccard similarity, and the inverse   
jacFallSN <- jaccard_per_row(DHFallSN_b, margin=2)

#Creating columns representing each of the two cameras
jacFallSN$cam1 <- gsub("_.*", "", row.names(jacFallSN))
jacFallSN$cam2 <- substr(row.names(jacFallSN), 10, 13)

#removes comparisons of same camera to same camera
jacFallSN <-jacFallSN[which(jacFallSN$cam1 != jacFallSN$cam2), ]

#Pull intercamera distances from distance matrix based on value of the two cameras
jacFallSN$dist <- NA

for (i in 1:nrow(jacFallSN)) {
  jacFallSN$dist[i] <- distMatS_named[as.character(as.numeric(jacFallSN$cam1[i])),as.character(as.numeric(jacFallSN$cam2[i]))]
  
}

#Graph the similarity index by distance

foxSqJacPlot_Fall <- ggplot(data = jacFallSN, mapping = aes(x = dist, y = JSim)) +
  geom_point() +
  ylim(0, 1) +
  xlim(10,120)
foxSqJacPlot_Fall

# Deer Jaccard Calculations and Graph -----------------------------
#Need to remove cameras and ocassions that had NAs to make it work. First I pull out the first ocassion, because cameras were set on different days. Then I remove rows that have NAs.

#Spring Deer
DHSpOV_b <- as.data.frame(DHSpOV) %>%
  dplyr::select(-1) 

DHSpOV_b <- DHSpOV_b[which(complete.cases(DHSpOV_b)),]

#Calculates jaccard similarity, and the inverse   
jacSpOV <- jaccard_per_row(DHSpOV_b, margin=2)

#Creating columns representing each of the two cameras
jacSpOV$cam1 <- gsub("_.*", "", row.names(jacSpOV))
jacSpOV$cam2 <- substr(row.names(jacSpOV), 11, 14)

#removes comparisons of same camera to same camera
jacSpOV <-jacSpOV[which(jacSpOV$cam1 != jacSpOV$cam2), ]

#Pull intercamera distances from distance matrix based on value of the two cameras
jacSpOV$dist <- NA

for (i in 1:nrow(jacSpOV)) {
  jacSpOV$dist[i] <- distMatS_named[as.character(as.numeric(jacSpOV$cam1[i])),as.character(as.numeric(jacSpOV$cam2[i]))]
  
}

#Graph the similarity index by distance

DeerJacPlot_SP <- ggplot(data = jacSpOV, mapping = aes(x = dist, y = JSim)) +
  geom_point() +
  ylim(0, 1) +
  xlim(10,120)
DeerJacPlot_SP

#Summer Deer
DHSumOV_b <- as.data.frame(DHSumOV) %>%
  dplyr::select(-1) 

DHSumOV_b <- DHSumOV_b[which(complete.cases(DHSumOV_b)),]

#Calculates jaccard similarity, and the inverse   
jacSumOV <- jaccard_per_row(DHSumOV_b, margin=2)

#Creating columns representing each of the two cameras
jacSumOV$cam1 <- gsub("-.*", "", row.names(jacSumOV))
jacSumOV$cam2 <- substr(row.names(jacSumOV), 10, 13)

#removes comparisons of same camera to same camera
jacSumOV <-jacSumOV[which(jacSumOV$cam1 != jacSumOV$cam2), ]

#Pull intercamera distances from distance matrix based on value of the two cameras
jacSumOV$dist <- NA

for (i in 1:nrow(jacSumOV)) {
  jacSumOV$dist[i] <- distMatS_named[as.character(as.numeric(jacSumOV$cam1[i])),as.character(as.numeric(jacSumOV$cam2[i]))]
  
}

#Graph the similarity index by distance

DeerJacPlot_SUM <- ggplot(data = jacSumOV, mapping = aes(x = dist, y = JSim)) +
  geom_point() +
  ylim(0, 1) +
  xlim(10,120)
DeerJacPlot_SUM

#Winter Deer
DHWinOV_b <- as.data.frame(DHWinOV) %>%
  dplyr::select(-1) 

DHWinOV_b <- DHWinOV_b[which(complete.cases(DHWinOV_b)),]

#Calculates jaccard similarity, and the inverse   
jacWinOV <- jaccard_per_row(DHWinOV_b, margin=2)

#Creating columns representing each of the two cameras
jacWinOV$cam1 <- gsub("_.*", "", row.names(jacWinOV))
jacWinOV$cam2 <- substr(row.names(jacWinOV), 10, 13)

#removes comparisons of same camera to same camera
jacWinOV <-jacWinOV[which(jacWinOV$cam1 != jacWinOV$cam2), ]

#Pull intercamera distances from distance matrix based on value of the two cameras
jacWinOV$dist <- NA

for (i in 1:nrow(jacWinOV)) {
  jacWinOV$dist[i] <- distMatS_named[as.character(as.numeric(jacWinOV$cam1[i])),as.character(as.numeric(jacWinOV$cam2[i]))]
  
}

#Graph the similarity index by distance

DeerJacPlot_Win <- ggplot(data = jacWinOV, mapping = aes(x = dist, y = JSim)) +
  geom_point() +
  ylim(0, 1) +
  xlim(10,120)
DeerJacPlot_Win

#Fall Deer
DHFallOV_b <- as.data.frame(DHFallOV) %>%
  dplyr::select(-1) 

DHFallOV_b <- DHFallOV_b[which(complete.cases(DHFallOV_b)),]

#Calculates jaccard similarity, and the inverse   
jacFallOV <- jaccard_per_row(DHFallOV_b, margin=2)

#Creating columns representing each of the two cameras
jacFallOV$cam1 <- gsub("_.*", "", row.names(jacFallOV))
jacFallOV$cam2 <- substr(row.names(jacFallOV), 10, 13)

#removes comparisons of same camera to same camera
jacFallOV <-jacFallOV[which(jacFallOV$cam1 != jacFallOV$cam2), ]

#Pull intercamera distances from distance matrix based on value of the two cameras
jacFallOV$dist <- NA

for (i in 1:nrow(jacFallOV)) {
  jacFallOV$dist[i] <- distMatS_named[as.character(as.numeric(jacFallOV$cam1[i])),as.character(as.numeric(jacFallOV$cam2[i]))]
  
}

#Graph the similarity index by distance

DeerJacPlot_Fall <- ggplot(data = jacFallOV, mapping = aes(x = dist, y = JSim)) +
  geom_point() +
  ylim(0, 1) +
  xlim(10,120)
DeerJacPlot_Fall



# Raccoon Jaccard Calculations and Graph -----------------------------
#Need to remove cameras and ocassions that had NAs to make it work. First I pull out the first ocassion, because cameras were set on different days. Then I remove rows that have NAs.

#Spring Raccoon
DHSpPL_b <- as.data.frame(DHSpPL) %>%
  dplyr::select(-1) 

DHSpPL_b <- DHSpPL_b[which(complete.cases(DHSpPL_b)),]

#Calculates jaccard similarity, and the inverse   
jacSpPL <- jaccard_per_row(DHSpPL_b, margin=2)

#Creating columns representing each of the two cameras
jacSpPL$cam1 <- gsub("_.*", "", row.names(jacSpPL))
jacSpPL$cam2 <- substr(row.names(jacSpPL), 11, 14)

#removes comparisons of same camera to same camera
jacSpPL <-jacSpPL[which(jacSpPL$cam1 != jacSpPL$cam2), ]

#Pull intercamera distances from distance matrix based on value of the two cameras
jacSpPL$dist <- NA

for (i in 1:nrow(jacSpPL)) {
  jacSpPL$dist[i] <- distMatS_named[as.character(as.numeric(jacSpPL$cam1[i])),as.character(as.numeric(jacSpPL$cam2[i]))]
  
}

#Graph the similarity index by distance

RacJacPlot_SP <- ggplot(data = jacSpPL, mapping = aes(x = dist, y = JSim)) +
  geom_point() +
  ylim(0, 1) +
  xlim(10,120)
RacJacPlot_SP

#Summer Raccoon
DHSumPL_b <- as.data.frame(DHSumPL) %>%
  dplyr::select(-1) 

DHSumPL_b <- DHSumPL_b[which(complete.cases(DHSumPL_b)),]

#Calculates jaccard similarity, and the inverse   
jacSumPL <- jaccard_per_row(DHSumPL_b, margin=2)

#Creating columns representing each of the two cameras
jacSumPL$cam1 <- gsub("-.*", "", row.names(jacSumPL))
jacSumPL$cam2 <- substr(row.names(jacSumPL), 10, 13)

#removes comparisons of same camera to same camera
jacSumPL <-jacSumPL[which(jacSumPL$cam1 != jacSumPL$cam2), ]

#Pull intercamera distances from distance matrix based on value of the two cameras
jacSumPL$dist <- NA

for (i in 1:nrow(jacSumPL)) {
  jacSumPL$dist[i] <- distMatS_named[as.character(as.numeric(jacSumPL$cam1[i])),as.character(as.numeric(jacSumPL$cam2[i]))]
  
}

#Graph the similarity index by distance

RacJacPlot_SUM <- ggplot(data = jacSumPL, mapping = aes(x = dist, y = JSim)) +
  geom_point() +
  ylim(0, 1) +
  xlim(10,120)
RacJacPlot_SUM

#Winter Raccoon
DHWinPL_b <- as.data.frame(DHWinPL) %>%
  dplyr::select(-1) 

DHWinPL_b <- DHWinPL_b[which(complete.cases(DHWinPL_b)),]

#Calculates jaccard similarity, and the inverse   
jacWinPL <- jaccard_per_row(DHWinPL_b, margin=2)

#Creating columns representing each of the two cameras
jacWinPL$cam1 <- gsub("_.*", "", row.names(jacWinPL))
jacWinPL$cam2 <- substr(row.names(jacWinPL), 10, 13)

#removes comparisons of same camera to same camera
jacWinPL <-jacWinPL[which(jacWinPL$cam1 != jacWinPL$cam2), ]

#Pull intercamera distances from distance matrix based on value of the two cameras
jacWinPL$dist <- NA

for (i in 1:nrow(jacWinPL)) {
  jacWinPL$dist[i] <- distMatS_named[as.character(as.numeric(jacWinPL$cam1[i])),as.character(as.numeric(jacWinPL$cam2[i]))]
  
}

#Graph the similarity index by distance

RacJacPlot_Win <- ggplot(data = jacWinPL, mapping = aes(x = dist, y = JSim)) +
  geom_point() +
  ylim(0, 1) +
  xlim(10,120)
RacJacPlot_Win

#Fall Raccoon
DHFallPL_b <- as.data.frame(DHFallPL) %>%
  dplyr::select(-1) 

DHFallPL_b <- DHFallPL_b[which(complete.cases(DHFallPL_b)),]

#Calculates jaccard similarity, and the inverse   
jacFallPL <- jaccard_per_row(DHFallPL_b, margin=2)

#Creating columns representing each of the two cameras
jacFallPL$cam1 <- gsub("_.*", "", row.names(jacFallPL))
jacFallPL$cam2 <- substr(row.names(jacFallPL), 10, 13)

#removes comparisons of same camera to same camera
jacFallPL <-jacFallPL[which(jacFallPL$cam1 != jacFallPL$cam2), ]

#Pull intercamera distances from distance matrix based on value of the two cameras
jacFallPL$dist <- NA

for (i in 1:nrow(jacFallPL)) {
  jacFallPL$dist[i] <- distMatS_named[as.character(as.numeric(jacFallPL$cam1[i])),as.character(as.numeric(jacFallPL$cam2[i]))]
  
}

#Graph the similarity index by distance

RacJacPlot_Fall <- ggplot(data = jacFallPL, mapping = aes(x = dist, y = JSim)) +
  geom_point() +
  ylim(0, 1) +
  xlim(10,120)
RacJacPlot_Fall

# Bear Jaccard Calculations and Graph -----------------------------

#Summer Bear
DHSumUA_b <- as.data.frame(DHSumUA) %>%
  dplyr::select(-1) 

DHSumUA_b <- DHSumUA_b[which(complete.cases(DHSumUA_b)),]

#Calculates jaccard similarity, and the inverse   
jacSumUA <- jaccard_per_row(DHSumUA_b, margin=2)

#Creating columns representing each of the two cameras
jacSumUA$cam1 <- gsub("-.*", "", row.names(jacSumUA))
jacSumUA$cam2 <- substr(row.names(jacSumUA), 10, 13)

#removes comparisons of same camera to same camera
jacSumUA <-jacSumUA[which(jacSumUA$cam1 != jacSumUA$cam2), ]

#Pull intercamera distances from distance matrix based on value of the two cameras
jacSumUA$dist <- NA

for (i in 1:nrow(jacSumUA)) {
  jacSumUA$dist[i] <- distMatS_named[as.character(as.numeric(jacSumUA$cam1[i])),as.character(as.numeric(jacSumUA$cam2[i]))]
  
}

#Graph the similarity index by distance

BearJacPlot_SUM <- ggplot(data = jacSumUA, mapping = aes(x = dist, y = JSim)) +
  geom_point() +
  ylim(0, 1) +
  xlim(10,120)
BearJacPlot_SUM


#Bear Raccoon
DHFallUA_b <- as.data.frame(DHFallUA) %>%
  dplyr::select(-1) 

DHFallUA_b <- DHFallUA_b[which(complete.cases(DHFallUA_b)),]

#Calculates jaccard similarity, and the inverse   
jacFallUA <- jaccard_per_row(DHFallUA_b, margin=2)

#Creating columns representing each of the two cameras
jacFallUA$cam1 <- gsub("_.*", "", row.names(jacFallUA))
jacFallUA$cam2 <- substr(row.names(jacFallUA), 10, 13)

#removes comparisons of same camera to same camera
jacFallUA <-jacFallUA[which(jacFallUA$cam1 != jacFallUA$cam2), ]

#Pull intercamera distances from distance matrix based on value of the two cameras
jacFallUA$dist <- NA

for (i in 1:nrow(jacFallUA)) {
  jacFallUA$dist[i] <- distMatS_named[as.character(as.numeric(jacFallUA$cam1[i])),as.character(as.numeric(jacFallUA$cam2[i]))]
  
}

#Graph the similarity index by distance

BearJacPlot_Fall <- ggplot(data = jacFallUA, mapping = aes(x = dist, y = JSim)) +
  geom_point() +
  ylim(0, 1) +
  xlim(10,120)
BearJacPlot_Fall
