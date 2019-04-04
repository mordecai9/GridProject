#It looks like we are missing 0102_W17. Need to check to see where that is. 

library(tidyverse)

load("data/camdata_summary")

#Create species by season files for focal species. 

# Data Import and Exploration - WTD ---------------------------------------

load("data/deerDataSum.RData")
load("data/deerDataFall.RData")
load("data/deerDataWin.RData")
load("data/deerDataSpr.RData")

deerSumCR <- deerDataSum %>%
  dplyr::select(c("Deployment", "CR","NAD83_X", "NAD83_Y"))

deerWinCR <- deerDataWin %>%
  dplyr::select(c("Deployment", "CR","NAD83_X", "NAD83_Y"))

deerSprCR <- deerDataSpr %>%
  dplyr::select(c("Deployment", "CR","NAD83_X", "NAD83_Y"))

deerFallCR <- deerDataFall %>%
  dplyr::select(c("Deployment", "CR","NAD83_X", "NAD83_Y"))

# Data Import and Exploration - Bear ---------------------------------------
load("data/bearDataSum.RData")
load("data/bearDataFall.RData")
load("data/bearDataWin.RData")
load("data/bearDataSpr.RData")

bearSumCR <- bearDataSum %>%
  dplyr::select(c("Deployment", "CR","NAD83_X", "NAD83_Y"))

bearWinCR <- bearDataWin %>%
  dplyr::select(c("Deployment", "CR","NAD83_X", "NAD83_Y"))

bearSprCR <- bearDataSpr %>%
  dplyr::select(c("Deployment", "CR","NAD83_X", "NAD83_Y"))

bearFallCR <- bearDataFall %>%
  dplyr::select(c("Deployment", "CR","NAD83_X", "NAD83_Y"))

# Data Import and Exploration - Raccoon ---------------------------------------
load("data/racDataSum.RData")
load("data/racDataFall.RData")
load("data/racDataWin.RData")
load("data/racDataSpr.RData")

#raccoon
raccoonSumCR <- racDataSum %>%
  dplyr::select(c("Deployment", "CR","NAD83_X", "NAD83_Y"))

raccoonWinCR <- racDataWin %>%
  dplyr::select(c("Deployment", "CR","NAD83_X", "NAD83_Y"))

raccoonSprCR <- racDataSpr %>%
  dplyr::select(c("Deployment", "CR","NAD83_X", "NAD83_Y"))

raccoonFallCR <- racDataFall %>%
  dplyr::select(c("Deployment", "CR","NAD83_X", "NAD83_Y"))

# Data Import and Exploration - All Squirrels ---------------------------------------

load("data/sqDataSum.RData")
load("data/sqDataFall.RData")
load("data/sqDataWin.RData")
load("data/sqDataSpr.RData")

squirrelSumCR <- sqDataSum %>%
  dplyr::select(c("Deployment", "CR","NAD83_X", "NAD83_Y"))

squirrelWinCR <- sqDataWin %>%
  dplyr::select(c("Deployment", "CR","NAD83_X", "NAD83_Y"))

squirrelSprCR <- sqDataSpr %>%
  dplyr::select(c("Deployment", "CR","NAD83_X", "NAD83_Y"))

squirrelFallCR <- sqDataFall %>%
  dplyr::select(c("Deployment", "CR","NAD83_X", "NAD83_Y"))



# Correlograms ------------------------------------------------------------

library(ncf)
d <- 80 #max distance for correlograms

#White Tailed Deer. Will need to format as high res tiff if used in paper.
par.default <- par(no.readonly = T)
par(mfrow = c(2,2))

Correlog_DeerSum <- spline.correlog(x = deerDataSum[, "NAD83_X"],
                                    y = deerDataSum[, "NAD83_Y"],
                                    z = deerDataSum[, "CR"],
                                    xmax = d)
plot(Correlog_DeerSum)
abline(v = 20, lty = 2)
text(x = 1, y = 0.8, "A", cex = 2)

Correlog_DeerFall <- spline.correlog(
  x = deerDataFall[, "NAD83_X"],
  y = deerDataFall[, "NAD83_Y"],
  z = deerDataFall[, "CR"],
  xmax = d)
plot(Correlog_DeerFall)
abline(v = 20, lty = 2)
text(x = 1, y = 0.8, "B", cex = 2)

Correlog_DeerWin <- spline.correlog(
  x = deerDataWin[, "NAD83_X"],
  y = deerDataWin[, "NAD83_Y"],
  z = deerDataWin[, "CR"],
  xmax = d)
plot(Correlog_DeerWin)
abline(v = 20, lty = 2)
text(x = 1, y = 0.8, "C", cex = 2)

Correlog_DeerSpr <- spline.correlog(x = deerDataSpr[, "NAD83_X"],
                                    y = deerDataSpr[, "NAD83_Y"],
                                    z = deerDataSpr[, "CR"],
                                    xmax = d)
plot(Correlog_DeerSpr)
abline(v = 20, lty = 2)
text(x = 1, y = 0.8, "D", cex = 2)

par(par.default)

#Black Bear, ignoring Winter here

par.default <- par(no.readonly = T)
par(mfrow = c(2,2))

Correlog_BearSum <- spline.correlog(x = bearDataSum[, "NAD83_X"],
                                    y = bearDataSum[, "NAD83_Y"],
                                    z = bearDataSum[, "CR"],
                                    xmax = d)
plot(Correlog_BearSum, main = "Bear/Summer")
abline(v = 20, lty = 2)

Correlog_BearFall <- spline.correlog(
  x = bearDataFall[, "NAD83_X"],
  y = bearDataFall[, "NAD83_Y"],
  z = bearDataFall[, "CR"],
  xmax = d)
plot(Correlog_BearFall, main = "Bear/Fall")
abline(v = 20, lty = 2)



Correlog_BearSpr <- spline.correlog(x = bearDataSpr[, "NAD83_X"],
                                    y = bearDataSpr[, "NAD83_Y"],
                                    z = bearDataSpr[, "CR"],
                                    xmax = d)
plot(Correlog_BearSpr, main = "Bear/Spring")
abline(v = 20, lty = 2)

par(par.default)

#Raccoons
par.default <- par(no.readonly = T)
par(mfrow = c(2,2))

Correlog_RacSum <- spline.correlog(x = racDataSum[, "NAD83_X"],
                                    y = racDataSum[, "NAD83_Y"],
                                    z = racDataSum[, "CR"],
                                    xmax = d)
plot(Correlog_RacSum, main = "Raccoon/Summer")
abline(v = 20, lty = 2)

Correlog_RacFall <- spline.correlog(
  x = racDataFall[, "NAD83_X"],
  y = racDataFall[, "NAD83_Y"],
  z = racDataFall[, "CR"],
  xmax = d)
plot(Correlog_RacFall, main = "Raccoon/Fall")
abline(v = 20, lty = 2)

Correlog_RacWin <- spline.correlog(
  x = racDataWin[, "NAD83_X"],
  y = racDataWin[, "NAD83_Y"],
  z = racDataWin[, "CR"],
  xmax = d)
plot(Correlog_RacWin,   main = "Raccoon/Winter")
abline(v = 20, lty = 2)

Correlog_RacSpr <- spline.correlog(x = racDataSpr[, "NAD83_X"],
                                    y = racDataSpr[, "NAD83_Y"],
                                    z = racDataSpr[, "CR"],
                                    xmax = d)
plot(Correlog_RacSpr, main = "Raccoon/Spring")
abline(v = 20, lty = 2)

par(par.default)

#Squirrels
par.default <- par(no.readonly = T)
par(mfrow = c(2,2))

Correlog_SqSum <- spline.correlog(x = sqDataSum[, "NAD83_X"],
                                   y = sqDataSum[, "NAD83_Y"],
                                   z = sqDataSum[, "CR"],
                                   xmax = d)
plot(Correlog_SqSum, main = "Squirrel/Summer")
abline(v = 20, lty = 2)

Correlog_SqFall <- spline.correlog(
  x = sqDataFall[, "NAD83_X"],
  y = sqDataFall[, "NAD83_Y"],
  z = sqDataFall[, "CR"],
  xmax = d)
plot(Correlog_SqFall, main = "Squirrel/Fall")
abline(v = 20, lty = 2)

Correlog_SqWin <- spline.correlog(
  x = sqDataWin[, "NAD83_X"],
  y = sqDataWin[, "NAD83_Y"],
  z = sqDataWin[, "CR"],
  xmax = d)
plot(Correlog_SqWin,   main = "Squirrel/Winter")
abline(v = 20, lty = 2)

Correlog_SqSpr <- spline.correlog(x = sqDataSpr[, "NAD83_X"],
                                   y = sqDataSpr[, "NAD83_Y"],
                                   z = sqDataSpr[, "CR"],
                                   xmax = d)
plot(Correlog_SqSpr, main = "Squirrel/Spring")
abline(v = 20, lty = 2)

par(par.default)

#________________________________________________
#Calculating Moran's I - Spatial Autocorrelation
#________________________________________________

library(ape) #for Moran's I
library(geoR) #for variograms

#Generate distance matrix, one per season because of differing # of active points
#Then take inverse of distance matrix and fill diagonals with 0
distMatS<-as.matrix(dist(cbind(deerSumCR$NAD83_X, deerSumCR$NAD83_Y)))

#Checking min and max intercamera distance, which is 15.82 and 114.7292
diag(distMatS) <- NA
min(distMatS[!is.na(distMatS)])
max(distMatS[!is.na(distMatS)])

distMatS<-as.matrix(dist(cbind(deerSumCR$NAD83_X, deerSumCR$NAD83_Y))) #return to normal

invDistMatS <- 1/distMatS
diag(invDistMatS)<-0

distMatF<-as.matrix(dist(cbind(deerFallCR$NAD83_X, deerFallCR$NAD83_Y)))
invDistMatF <- 1/distMatF
diag(invDistMatF)<-0

distMatW<-as.matrix(dist(cbind(deerWinCR$NAD83_X, deerWinCR$NAD83_Y)))
invDistMatW <- 1/distMatW
diag(invDistMatW)<-0

distMatSp<-as.matrix(dist(cbind(deerSprCR$NAD83_X, deerSprCR$NAD83_Y)))
invDistMatSp <- 1/distMatSp
diag(invDistMatSp)<-0


#Use coords to create XY matrix
XYmatrixS<-as.matrix(deerSumCR[,c(3,4)])
XYmatrixF <- as.matrix(deerFallCR[,c(3,4)])
XYmatrixSp <- as.matrix(deerSprCR[,c(3,4)])
XYmatrixW <- as.matrix(deerWinCR[,c(3,4)])

#_____________________________________
#Deer Analysis for Autocorrelation####
#_____________________________________
#Deer_MI_inv[1:5, 1:5] #Need to ask Josey why she did this

#Calculate Moran's I for Deer for all Seasons
deerSumMoran<-Moran.I(deerSumCR$CR, invDistMatS)
deerFallMoran<-Moran.I(deerFallCR$CR, invDistMatF)
deerWinMoran<-Moran.I(deerWinCR$CR, invDistMatW) #almost sig at 0.078
deerSprMoran<-Moran.I(deerSprCR$CR, invDistMatSp)

#Deer Variograms for 4 Seasons. Need to remove the cameras with no data from Fall and Winter to have it work
#Need to read more about setting the right bin size and how far we can interpret teh variogram (perhaps only half the field size?)
#The scaled option doesn't seem to work, as everythign just becomes zero for some reason
#A typical semivariogram would start near zero, with there being very little variance between points that are close together. Then the graph would increase until a sill or plateau, reaching some sort of background baseline variance, indicating no spatial autocorrelation. Here I'm not seeing that pattern. This indicates either 1) the points are all spatially correlated at the same level and the cameras aren't far enough to diminish this or 2) 20meters is enough to have independent detection rates of these species. We essentially above the sill (which may be the data variance). The latter is supported by the fact that Moran's I is not significant for almost all of these species with some exceptions.
par(mfrow = c(2,2), mar = c(5,5,5,5))
DvarioS<-variog(coords = XYmatrixS, data = deerSumCR$CR, breaks = seq(0, 150,10), 
               option = "bin")
plot(DvarioS, xlim = c(0,120), ylim = c(-2000,2000), main = "Summer", pts.range.cex = c(1,3))

DvarioW<-variog(coords = XYmatrixW, data = deerWinCR$CR, breaks = seq(0, 150,15), 
                option = "bin")
plot(DvarioW, xlim = c(0,120), ylim = c(-1000,1000), main = "Winter", scaled = T, var.lines = T)

DvarioSp<-variog(coords = XYmatrixSp, data = deerSprCR$CR, breaks = seq(0, 150,15), 
                option = "bin")
plot(DvarioSp, xlim = c(0,120), ylim = c(-1000,1000), main = "Spring", var.lines = T, scaled = T, pts.range.cex = c(1,3))

DvarioF<-variog(coords = XYmatrixF, data = deerFallCR$CR, breaks = seq(0, 150,15), 
                option = "bin", bin.cloud = T)
plot(DvarioF, xlim = c(0,120), ylim = c(0,20000), main = "Fall", var.lines = T)

#another plot format option
#plot(DvarioF, xlim = c(0,120), ylim = c(0,15000), main = "Fall", var.lines = T, #scaled = T, bin.cloud = T)



#################################
#Raccoon Analysis####
#Calculate Moran's I for bear for all Seasons. None are even close to significant
raccoonSumMoran<-Moran.I(raccoonSumCR$CR, invDistMatS)
raccoonFallMoran<-Moran.I(raccoonFallCR$CR, invDistMatF)
raccoonWinMoran<-Moran.I(raccoonWinCR$CR, invDistMatW)
raccoonSprMoran<-Moran.I(raccoonSprCR$CR, invDistMatSp)

#raccoon Variograms for 4 Seasons. Need to remove the cameras with no data from Fall and Winter to have it work
#Need to read more about setting the right bin size and how far we can interpret teh variogram (perhaps only half the field size?)
par(mfrow = c(2,2), mar = c(5,5,5,5))
RvarioS<-variog(coords = XYmatrixS, data = raccoonSumCR$CR, breaks = seq(0, 150,20), 
                option = "bin")
plot(RvarioS, xlim = c(0,120), ylim = c(0,400), main = "Summer", var.lines =T)

RvarioW<-variog(coords = XYmatrixW, data = raccoonWinCR$CR, breaks = seq(0, 150,20), 
                option = "bin")
plot(RvarioW, xlim = c(0,120), ylim = c(0,2000), main = "Winter", var.lines = T)

RvarioSp<-variog(coords = XYmatrixSp, data = raccoonSprCR$CR, breaks = seq(0, 150,20), 
                 option = "bin")
plot(RvarioSp, xlim = c(0,120), ylim = c(0,400), main = "Spring", var.lines = T)

RvarioF<-variog(coords = XYmatrixF, data = raccoonFallCR$CR, breaks = seq(0, 150,20), 
                option = "bin")
plot(RvarioF, xlim = c(0,120), ylim = c(0,400), main = "Fall", var.lines = T)

#################################
#Bear Analysis####
#Calculate Moran's I for Raccoon for all Seasons. There may not be enough captures in winter or spring to warrant running those tests.
bearSumMoran<-Moran.I(bearSumCR$CR, invDistMatS)
bearFallMoran<-Moran.I(bearFallCR$CR, invDistMatF)
bearSprMoran<-Moran.I(bearSprCR$CR, invDistMatSp) #nearly sig 0.094 

#bear Variograms for 4 Seasons. Need to remove the cameras with no data from Fall and Winter to have it work
#Need to read more about setting the right bin size and how far we can interpret teh variogram (perhaps only half the field size?)
par(mfrow = c(2,2), mar = c(5,5,5,5))
BvarioS<-variog(coords = XYmatrixS, data = bearSumCR$CR, breaks = seq(0, 150,20), 
                option = "bin")
plot(BvarioS, xlim = c(0,120), ylim = c(0,150), main = "Summer")

BvarioW<-variog(coords = XYmatrixW, data = bearWinCR$CR, breaks = seq(0, 150,20), 
                option = "bin")
plot(BvarioW, xlim = c(0,120), ylim = c(0,150), main = "Winter")

BvarioSp<-variog(coords = XYmatrixSp, data = bearSprCR$CR, breaks = seq(0, 150,20), 
                 option = "bin")
plot(BvarioSp, xlim = c(0,120), ylim = c(0,150), main = "Spring")

BvarioF<-variog(coords = XYmatrixF, data = bearFallCR$CR, breaks = seq(0, 150,20), 
                option = "bin")
plot(BvarioF, xlim = c(0,120), ylim = c(0,150), main = "Fall")

#################################
#Squirrel Analysis####
#Calculate Moran's I for squirrel for all Seasons. 

sqSumMoran<-Moran.I(squirrelSumCR$CR, invDistMatS)
sqFallMoran<-Moran.I(squirrelFallCR$CR, invDistMatF)
sqWinMoran<-Moran.I(squirrelWinCR$CR, invDistMatW)
sqSprMoran<-Moran.I(squirrelSprCR$CR, invDistMatSp)

#########################################
#Plotting Covariance
#########################################
#Calculating Deer Covariance
#Dcr = Deer_MI$Capture_Rate
#Ddist<-as.data.frame(Deervario)
#cov(Dcr, Ddist)