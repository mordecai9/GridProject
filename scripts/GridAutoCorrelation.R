#It looks like we are missing 0102_W17. Need to check to see where that is. 

library(tidyverse)

load("data/camdata_summary")

#Create species by season files for focal species. Think about creating these in RegressionPrep file and using them here? May have to create them again there, and I need to pool squirrels here, like I do there.
#Deer
deerSumCR <- camdata_summary %>%
  filter(Species == "Odocoileus virginianus" & Season == "Summer 2017") %>%
  dplyr::select(c("Deployment_Name2", "CR","NAD83_X", "NAD83_Y"))

deerWinCR <- camdata_summary %>%
  filter(Species == "Odocoileus virginianus" & Season == "Winter 2017") %>%
  dplyr::select(c("Deployment_Name2", "CR","NAD83_X", "NAD83_Y"))

deerSprCR <- camdata_summary %>%
  filter(Species == "Odocoileus virginianus" & Season == "Spring 2018") %>%
  dplyr::select(c("Deployment_Name2", "CR","NAD83_X", "NAD83_Y"))

deerFallCR <- camdata_summary %>%
  filter(Species == "Odocoileus virginianus" & Season == "Fall 2017") %>%
  dplyr::select(c("Deployment_Name2", "CR","NAD83_X", "NAD83_Y"))
#Bear
bearSumCR <- camdata_summary %>%
  filter(Species == "Ursus americanus" & Season == "Summer 2017") %>%
  dplyr::select(c("Deployment_Name2", "CR","NAD83_X", "NAD83_Y"))

bearWinCR <- camdata_summary %>%
  filter(Species == "Ursus americanus" & Season == "Winter 2017") %>%
  dplyr::select(c("Deployment_Name2", "CR","NAD83_X", "NAD83_Y"))

bearSprCR <- camdata_summary %>%
  filter(Species == "Ursus americanus" & Season == "Spring 2018") %>%
  dplyr::select(c("Deployment_Name2", "CR","NAD83_X", "NAD83_Y"))

bearFallCR <- camdata_summary %>%
  filter(Species == "Ursus americanus" & Season == "Fall 2017") %>%
  dplyr::select(c("Deployment_Name2", "CR","NAD83_X", "NAD83_Y"))

#raccoon
raccoonSumCR <- camdata_summary %>%
  filter(Species == "Procyon lotor" & Season == "Summer 2017") %>%
  dplyr::select(c("Deployment_Name2", "CR","NAD83_X", "NAD83_Y"))

raccoonWinCR <- camdata_summary %>%
  filter(Species == "Procyon lotor" & Season == "Winter 2017") %>%
  dplyr::select(c("Deployment_Name2", "CR","NAD83_X", "NAD83_Y"))

raccoonSprCR <- camdata_summary %>%
  filter(Species == "Procyon lotor" & Season == "Spring 2018") %>%
  dplyr::select(c("Deployment_Name2", "CR","NAD83_X", "NAD83_Y"))

raccoonFallCR <- camdata_summary %>%
  filter(Species == "Procyon lotor" & Season == "Fall 2017") %>%
  dplyr::select(c("Deployment_Name2", "CR","NAD83_X", "NAD83_Y"))

#________________________________________________
#Calculating Moran's I - Spatial Autocorrelation
#________________________________________________

library(ape) #for Moran's I
library(geoR) #for variograms

#Generate distance matrix, one per season because of differing # of active points
#Then take inverse of distance matrix and fill diagonals with 0
distMatS<-as.matrix(dist(cbind(deerSumCR$NAD83_X, deerSumCR$NAD83_Y)))
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
bearWinMoran<-Moran.I(bearWinCR$CR, invDistMatW) #significant
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



#########################################
#Plotting Covariance
#########################################
#Calculating Deer Covariance
#Dcr = Deer_MI$Capture_Rate
#Ddist<-as.data.frame(Deervario)
#cov(Dcr, Ddist)