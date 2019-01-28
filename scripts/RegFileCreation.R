#Create Individual Species/Seasons Data frames with SeqCount and Covariates to be used in Poisson GLM regression analysis. This file starts working with the object "camdata_summary" which is created and cleaned in previous scripts including CamDataPrep.R


require(tidyverse)

load("data/camdata_summary")
camdata <- camdata_summary
names(camdata)[1] <- "Deployment" #Facilitates merges later
names(camdata)[4] <- "nSeqs" #making shorter name here since this will be response variable

#create individual species data frames with species name, capture rate, and coordinates, In reality we will only need squirrels, deer, bear and raccoon, but the others are here if needed.
deerData<-subset(camdata, Species == "Odocoileus virginianus")
bearData<-subset(camdata, Species == "Ursus americanus")
coyoteData<-subset(camdata, Species == "Canis latrans")
foxsqrlData<-subset(camdata, Species == "Sciurus niger")
redfoxData<-subset(camdata, Species == "Vulpes vulpes")
bobcatData<-subset(camdata, Species == "Lynx rufus")
grsqrlData<-subset(camdata, Species == "Sciurus carolinensis")
raccoonData<-subset(camdata, Species == "Procyon lotor")
opossumData<-subset(camdata, Species =="Didelphis virginiana")
unknownsqrlData<-subset(camdata, Species == "Unknown Squirrel")


#For calculating Effective Detection Distance(EDD) all squirrel species were merged.
#Since EDD will be a variable, the squirrel species will be merged now

sqrlData<- foxsqrlData %>%
  dplyr::select(Deployment, Deployment_Name, Species, nSeqs, CR) %>%
  left_join(grsqrlData, by = "Deployment_Name") %>%
  left_join(unknownsqrlData, by = "Deployment_Name" )

sqrlData$allSqSeqs <- (sqrlData$nSeqs.x + sqrlData$nSeqs.y + sqrlData$nSeqs)
sqrlData$allSqCR <- sqrlData$allSqSeqs/sqrlData$Deploy.Duration.x *100
sqrlData <- sqrlData[, -c(3:8,11:12, 15:25)]
sqrlData$Species <- "All Squirrels"
names(sqrlData) <- c("Deployment","Deployment_Name", "Season", "Deploy.Duration", "NAD83_X", "NAD83_Y", "nSeqs", "CR", "Species") 

####Bring in the different variables####

#Variable - Camera Height####
#Bring in Camera Height Data
#CAMERA HEIGHT is the distance in centimeters from the ground to the camera lens#

camHeight<-read.csv("data/Camera_heights.csv")
sqrlData <- merge(sqrlData, camHeight, by = "Deployment")
deerData <- merge(deerData, camHeight, by = "Deployment")
bearData <- merge(bearData, camHeight, by = "Deployment")
raccoonData <- merge(raccoonData, camHeight, by = "Deployment")


#Variable - Number of Trees in camera sight####
#Working with SIGEO tree grid information to try to associate with capture rates from our high resolution camera grid. Coordinates should be UTM Zone 17S

library(rgeos)
library(rgdal)
library(sp)
library(maptools)
library(raster)
library(grid)

#Bring in geo-reference tree data from entire SIGEO grid
trees<-read.csv("data/scbi.full2_vegdata.csv")

#Change data frame into a Spatial Points Data frame using coordinates of the cameras
coordinates(trees)<- c("NAD83_X", "NAD83_Y")
class(trees)
#plot(SIGEOtrees) #this is a big file, and includes all trees in SIGEO plot, not just in our small grid

#plot the camera trap coordinates to show a point where the camera traps are in the grid. note that this plots each point many times. 
plot(camdata$NAD83_X,
     camdata$NAD83_Y,
     xlim = c(747420, 747560),
     ylim = c(4308900,4309040))

#Convert this grid coordinate information into a spatialpoints object so it is not plotted on a graph
#First need to make the x and y coordinates a separate matrix
trapxy <- camdata %>%
  dplyr::select(NAD83_X, NAD83_Y) 
  

trapxySP <- SpatialPointsDataFrame(coords = trapxy, data = camdata,
                                         proj4string = CRS(proj4string(trees)))
plot(trapxySP)

#Create a clipping polygon to reduce the size of the SIGEO grid to just the area of interest
#I'm setting the extent as 50m around the extreme trap coordinates
c<-50
clp <- as(extent(min(trapxy$NAD83_X)-c, 
                     max(trapxy$NAD83_X)+c,
                     min(trapxy$NAD83_Y)-c,
                     max(trapxy$NAD83_Y)+c),
              "SpatialPolygons")

#Assign the coordinate reference system of "trees" to the new clipping polygon         
proj4string(clp) <- CRS(proj4string(trees))
plot(clp)

#You could also use gIntersect below but it does not preserve the original attribute data
treesSmall <- intersect(trees, clp)

#plot grid with tree and cameras
par(mar = c(3, 9, 3, 2.1), xpd = T) #xpd = T allows text to exist outside the plot area
plot(treesSmall, col = "darkgreen", pch = 3,cex.main = 4)
plot(trapxySP, pch = 19, col = "red", add = T)
legend(747340,4308970, legend = c("Tree", "Camera"), col = c("darkgreen", "red"), pch = c(3,19), cex =1.0, bty = "n")

#Add scale, try to move to put under legend
scale.len <- 20
x1 <- c(747338.5,747338.5+scale.len)
y1<- c(4308890, 4308890)
lines(x1,y1,lwd = 2)
text(747372.9, 4308890, '20m', cex = 1.0)

#Add Deployment label to each camera
pointLabel(unique(coordinates(trapxySP)[,1]),unique(coordinates(trapxySP)[,2]),labels=as.character(unique(trapxySP@data$Deployment)), cex = 0.7, allowSmallOverlap = T)

#Create 4 point polygon to represent camera view of each camera
#Create data frame of the 4 points per camera
#These numbers were chosen so the polygon cone extends 20 meters north of each camera point and covers a 40 degree view angle from the camera trap
camview <- camdata[, c(1,9,10)]
camview$X1<-(camview$NAD83_X + 6.84)
camview$Y1<-(camview$NAD83_Y + 18.79)
camview$X2<-(camview$NAD83_X)
camview$Y2<-(camview$NAD83_Y + 20)
camview$X3<-(camview$NAD83_X - 6.84)
camview$Y3<-(camview$NAD83_Y + 18.79)

#Remove repeated lines
camview1<- dplyr::distinct(camview)

#All XY coordinates need to be in the same X and Y column (three columns total).
#This code will move all the points to the original NAD83 X and Y columns based on deloyment name
#Each camera will have 4 rows and each row will have a new XY coordinate point

#Make a seperate df for each set of X,Y coords
cam_poly_points<-camview1 [,c(1:3)]
cam_poly_points1<-camview1 [,c(1,4:5)]
names(cam_poly_points1) <- c("Deployment", "NAD83_X","NAD83_Y")
cam_poly_points2<-camview1 [,c(1,6:7)]
names(cam_poly_points2) <- c("Deployment", "NAD83_X","NAD83_Y")
cam_poly_points3<-camview1 [,c(1,8:9)]
names(cam_poly_points3) <- c("Deployment", "NAD83_X","NAD83_Y")

#Rbind all df of coords together
#Each camera will have 4 rows and each row will have a new XY coordinate pair
#I think this can be done in 1 line just listing all objects to rbind.
camview2<-rbind(cam_poly_points,cam_poly_points1)
camview3<-rbind(camview2,cam_poly_points2)
camview4<-rbind(camview3,cam_poly_points3)

camview_list<-split(camview4, camview1$Deployment)
camview5<-lapply(camview_list, function(x) {x["Deployment"]<- NULL; x})

#create sp object and convert coords to polygon to prepare for plotting
cvpp <- lapply(camview5, Polygon)

#add id variable
#this will help group individuals later
cvp<-lapply(seq_along(cvpp), function(i) Polygons(list(cvpp[[i]]),ID = names(camview_list)[i]))

#create sp object
camview_spo<-SpatialPolygons(cvp, proj4string = CRS(proj4string(trees)))

#Create spdf with IDs (one unique ID per poly) and plot polygons
camview_spo.df<-SpatialPolygonsDataFrame(camview_spo,data.frame(id = unique(camview1$Deployment),row.names = unique(camview1$Deployment)))
plot(camview_spo.df, add = T)

#Cut out tree data from within polygons
clip_polys<-intersect(trees,camview_spo.df)
plot(clip_polys)
plot(camview_spo.df, add = T)
cvtrees<-as.data.frame(clip_polys)

#Pull and total the # of trees per deployment and change column names
cvtreecount<-cvtrees[,c(4,28)]
cvtreecount1<-aggregate(cvtreecount[,1], by = list(cvtreecount$d),sum)
colnames(cvtreecount1)[2]<-"Num_Stems"
colnames(cvtreecount1)[1]<-"Deployment"

#Merge this tree data into each species file here
sqrlData <- merge(sqrlData, cvtreecount1, by = "Deployment")
deerData <- merge(deerData, cvtreecount1, by = "Deployment")
bearData <- merge(bearData, cvtreecount1, by = "Deployment")
raccoonData <- merge(raccoonData, cvtreecount1, by = "Deployment")



#____________________________
#Variable - Oak Trees####
#____________________________
#Pull Oak Tree Data from grid
Oaks <-subset(treesSmall,sp %in% c('qual','quru','quco','qufa','qupr','quve','qusp','qumi'))
plot(Oaks, pch = 19)
#plot camera locations in red
plot(trapxySP, pch = 1, col = "red", add = T)

#add column to study site tree info that divides trees into 5 color size cateories
Oaks$Size_Category[Oaks$dbh <150]<-'461' #turqoise
Oaks$Size_Category[Oaks$dbh >150]<-'68' #dark blue
Oaks$Size_Category[Oaks$dbh >300]<-'47' #yellow
Oaks$Size_Category[Oaks$dbh >600]<-'139' #green
Oaks$Size_Category[Oaks$dbh >900]<-'8' #gray
Oaks$Size_Category[Oaks$dbh >1200]<-'550' #pink

#plot Oak Tree sizes by color
par(mar=c(5,17,4,2))
plot(Oaks,pch = 19, col = Oaks$Size_Category)
plot(trapxySP, pch = 1, col = "red", add = T)
#Legend matching color to size (not working for some reason)
legend(747285,4309044, legend = c("< 15cm","> 15cm","> 30cm","> 60cm","> 90cm","> 120cm"), col = c("461", "68", "47","139", "8", "550"), pch = 19, title = "DBH of Oak Trees", bty = 'n')


#Cut out oak tree data from within the cones
#library(rowr) #what is this for?
pardefault <- par(no.readonly = T)
par(mar = c(5.1,4.1,4.1,2.1))
polyoaktrees<-intersect(Oaks, camview_spo.df)
plot(polyoaktrees)
plot(camview_spo.df, add = T)
polyoaktreesdf<-as.data.frame(polyoaktrees)
par(pardefault)

#Pull # of oaks out of each deployment and rename columns to prepare for merge
oakcount<-polyoaktreesdf[,c(4,29)]
oakcount1<-aggregate(oakcount[,1],by=list(oakcount$d), sum)
colnames(oakcount1)[2]<-"Num_Oaks"
colnames(oakcount1)[1]<-"Deployment"

#Pull DBH of oaks from each deployment and add total
Oak_DBH<-polyoaktreesdf[,c(11,29)]
Oak_DBH1<-aggregate(Oak_DBH[,1],by=list(Oak_DBH$d), sum)
colnames(Oak_DBH1)[1]<-"Deployment"
colnames(Oak_DBH1)[2]<-"OakDBH"

#Merge oak information into each species file

sqrlData <- merge(sqrlData, Oak_DBH1, by = "Deployment", all = T)
sqrlData[is.na(sqrlData)] <- 0
deerData <- merge(deerData, Oak_DBH1, by = "Deployment", all = T)
deerData[is.na(deerData)] <- 0
bearData <- merge(bearData, Oak_DBH1, by = "Deployment", all = T)
bearData[is.na(bearData)] <- 0
raccoonData <- merge(raccoonData, Oak_DBH1, by = "Deployment", all = T)
raccoonData[is.na(raccoonData)] <- 0


# Variable - Log in View Y/N ----------------------------------------------

camOpS<-read.csv("data/CameraOperation_SI17.csv")

logView <- camOpS %>%
  dplyr::select(c(Log.in.View, Deployment_Name2)) %>%
  rename(Deployment = Deployment_Name2)

sqrlData <- merge(sqrlData, logView, by = "Deployment")
deerData <- merge(deerData, logView, by = "Deployment")
bearData <- merge(bearData, logView, by = "Deployment")
raccoonData <- merge(raccoonData, logView, by = "Deployment")

# Variable - Effective Detection Distance ---------------------------------
#Bring in Detection Distance Data for Each Species
#We do not have unique detection distance data for each species in each season, although that would be ideal. We had to pool across seasons for most, and in some cases use data from other species. Details are below for each species.

#Deer EDD Data. For white-tailed deer, it was necessary to pool Summer/Fall data to get EDD information for Summer and Fall. Here we had a minimum # of observations of 31 for each camera. We also pooled Winter and Spring data for White-tailed deer, and here unfortunately had 5 cameras with less than 30 observations (min 18)

Deer_EDD<-read.csv("data/DeerEDD_4S.csv")
#Might be interesting here to plot the EDDs for the two season blocks to compare them. At first glance it looks like there is consistent trend from SummFall to WintSpring
deerData <- merge(deerData, Deer_EDD, by = "Deployment")

boxplot(Deer_EDD[, 2:3]) #So this makes sense, overall the cameras can see further when the vegetation is not there. Though I though it would be a bigger difference
plot(Summer.Fall.EDD ~ Winter.Spring.EDD, data = Deer_EDD) # You would think this would show a much stronger correlation. The cameras with long detection distances in Summer/Fall should have longer detection distances in Winter/Spring. Generally it is true but the relationship is not strong.
lm1 <- lm(Summer.Fall.EDD ~ Winter.Spring.EDD, data = Deer_EDD)
summary(lm1)
abline(lm1)

#Create new dataframes specific to each season, removing the EDD data from the opposite season
deerDataSum <- deerData %>%
  filter(Season == "Summer 2017") %>%
  dplyr::select(-Winter.Spring.EDD)
save(deerDataSum, file = "data/deerDataSum.RData")

deerDataFall <- deerData %>%
  filter(Season == "Fall 2017") %>%
  dplyr::select(-Winter.Spring.EDD)
save(deerDataFall, file = "data/deerDataFall.RData")

deerDataWin <- deerData %>%
  filter(Season == "Winter 2017") %>%
  dplyr::select(-Summer.Fall.EDD)
save(deerDataWin, file = "data/deerDataWin.RData")

deerDataSpr <- deerData %>%
  filter(Season == "Spring 2018") %>%
  dplyr::select(-Summer.Fall.EDD)
save(deerDataSpr, file = "data/deerDataSpr.RData")


#Bear EDD Data. For Black bear, we used EDD values for white-tailed deer. There were not sufficient data to estimate EDD for black bear in any season.

bearData <- merge(bearData, Deer_EDD, by = "Deployment")


bearDataSum <- bearData %>%
  filter(Season == "Summer 2017") %>%
  dplyr::select(-Winter.Spring.EDD)
save(bearDataSum, file = "data/bearDataSum.RData")

bearDataFall <- bearData %>%
  filter(Season == "Fall 2017") %>%
  dplyr::select(-Winter.Spring.EDD)
save(bearDataFall, file = "data/bearDataFall.RData")

bearDataWin <- bearData %>%
  filter(Season == "Winter 2017") %>%
  dplyr::select(-Summer.Fall.EDD)
save(bearDataWin, file = "data/bearDataWin.RData")

bearDataSpr <- bearData %>%
  filter(Season == "Spring 2018") %>%
  dplyr::select(-Summer.Fall.EDD)
save(bearDataSpr, file = "data/bearDataSpr.RData")

#Squirrel EDD Data. For squirrels, we had adequate observations to estimate detection distance for squirrels in Winter/Spring when these seasons were pooled, but had to use pooled 4-season data for Summer and Fall analyses for squirrels. 
#Waiting for clarification from Josey, because there is only one column of EDD data in the squirrel EDD csv file.

sqEDD <- read.csv("data/SquirrelEDD_4S.csv") 

sqrlData <- merge(sqrlData, sqEDD, by = "Deployment")

sqDataSum <- sqrlData %>%
  filter(Season == "Summer 2017") %>%
  dplyr::select(-Squirrel_EDD_WSp)
save(sqDataSum, file = "data/sqDataSum.RData")

sqDataFall <- sqrlData %>%
  filter(Season == "Fall 2017") %>%
  dplyr::select(-Squirrel_EDD_WSp)
save(sqDataFall, file = "data/sqDataFall.RData")

sqDataWin <- sqrlData %>%
  filter(Season == "Winter 2017") %>%
  dplyr::select(-Squirrel_EDD_4S)
save(sqDataWin, file = "data/sqDataWin.RData")

sqDataSpr <- sqrlData %>%
  filter(Season == "Spring 2018") %>%
  dplyr::select(-Squirrel_EDD_4S)
save(sqDataSpr, file = "data/sqDataSpr.RData")

#Raccoon EDD Data. For raccoons, when pooling across all 4 seasons, approximately half the cameras failed to meet the 30 observation threshold. For those cameras with adequate observations we calculated EDD. For the remaining, we used the EDDs from DEER, which were determined to correlate better with raccoon EDD than squirrel. 

racEDD <- read.csv("data/RaccoonEDD_4S.csv") 

raccoonData <- merge(raccoonData, racEDD, by = "Deployment")

racDataSum <- raccoonData %>%
  filter(Season == "Summer 2017") 
save(racDataSum, file = "data/racDataSum.RData")

racDataFall <- raccoonData %>%
  filter(Season == "Fall 2017") 
save(racDataFall, file = "data/racDataFall.RData")

racDataWin <- raccoonData %>%
  filter(Season == "Winter 2017") 
save(racDataWin, file = "data/racDataWin.RData")

racDataSpr <- raccoonData %>%
  filter(Season == "Spring 2018") 
save(racDataSpr, file = "data/racDataSpr.RData")