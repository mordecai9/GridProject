#Create Individual Species/Seasons Data frames with SeqCount and Covariates to be used in Poisson GLM regression analysis. I think here I can use "camdatasummary", not CR3, since it has the raw number of sequences and the camera effort, as well as season and both kinds of deployment names

load("data/camdata_summary")
camdata <- camdata_summary
names(camdata)[1] <- "Deployment" #Facilitates merges later

#create individual species data frames with species name, capture rate, and coordinates
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
  dplyr::select(Deployment, Deployment_Name, Species, Number_of_Sequences, CR) %>%
  left_join(grsqrlData, by = "Deployment_Name") %>%
  left_join(unknownsqrlData, by = "Deployment_Name" )

sqrlData$allSqSeqs <- (sqrlData$Number_of_Sequences.x + sqrlData$Number_of_Sequences.y + sqrlData$Number_of_Sequences)
sqrlData$allSqCR <- sqrlData$allSqSeqs/sqrlData$Deploy.Duration.x *100
sqrlData <- sqrlData[, -c(3:8,11:12, 15:25)]
sqrlData$Species <- "All Squirrels"

####Bring in the different variables####
#We can merge in these new variables into each of our 4-5 focal species files individually, once we have the code working for all of them. Below I am doing it for the full camdata file. I'm not sure which is best really. The height, and both tree variables will not change per species or per season. The only thing season and species specific is the EDD, so we obviously can't easily merge the EDD data into the full camdata file.


#Variable - Camera Height####
#Bring in Camera Height Data
#CAMERA HEIGHT is the distance in centimeters from the ground to the camera lens#

camHeight<-read.csv("data/Camera_heights.csv")
camdata <- merge(camdata, camHeight, by = "Deployment")


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

#plot the camera trap coordinates to show a point where the camera traps are in the grid
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
plot(treesSmall, col = "darkgreen", pch = 3,cex.main = 4)
plot(trapxySP, pch = 19, col = "red", add = T)

#Add a legend (this needs to be tweaked, maybe we need to change plotting margins with "par"?)
#par(font = 2)
#legend(747300,4308970, legend = c("Tree", "Camera"), col = c("darkgreen", "red"), pch = c(3,19), cex =1.5, bty = "n")

#Add scale (this doesn't appear to be working)
#scale.len <- 20
#x1 <- c(747308.5,747308.5+scale.len)
#y1<- c(4308890, 4308890)
#lines(x1,y1,lwd = 2)
#text(747347.9, 4308890, '20m', cex = 1.5)

#Add Deployment label to each camera
pointLabel(unique(coordinates(trapxySP)[,1]),unique(coordinates(trapxySP)[,2]),labels=as.character(unique(trapxySP@data$Deployment)), cex = 0.7, allowSmallOverlap = T)

#Josey, we will have to review this code below. I looked through it but did not understand it all. Will need to update names of files for sure. camdata_coordinates_SI17 is now camdata, but it has many more rows than before, because it has all 4 seasons' data.

#Create 4 point polygon to represent camera view
#Create data frame of the 4 points per camera
camview_SI17 <- camdata_coordinates_SI17[, c(2,3,5)]
camview_SI17$X1<-(camview_SI17$NAD83_X + 6.84)
camview_SI17$Y1<-(camview_SI17$NAD83_Y + 18.79)
camview_SI17$X2<-(camview_SI17$NAD83_X)
camview_SI17$Y2<-(camview_SI17$NAD83_Y + 20)
camview_SI17$X3<-(camview_SI17$NAD83_X - 6.84)
camview_SI17$Y3<-(camview_SI17$NAD83_Y + 18.79)

camview1_SI17<- camdata_coordinates_SI17 [,c(2,3,5)]
camview1_SI17[28:54,]<-(camview_SI17[1:27, c(4:5,3)]) #Could not figure out why these rows were being called out in this way.
camview1_SI17[55:81,]<-(camview_SI17[1:27, c(6:7,3)])
camview1_SI17[82:108,]<-(camview_SI17[1:27, c(8:9,3)])

camview_list_SI17<-split(camview1_SI17, camview1_SI17$Deployment)
camview_list_SI17<-lapply(camview_list_SI17, function(x) {x["Deployment"]<- NULL; x})

#create sp object and convert coords to polygon to prepare for 
cvpp_SI17 <- lapply(camview_list_SI17, Polygon)

#add id variable
cvp_SI17<-lapply(seq_along(cvpp_SI17), function(i) Polygons(list(cvpp_SI17[[i]]),ID = names(camview_list_SI17)[i]))

#create spobject
camview_spo_SI17<-SpatialPolygons(cvp_SI17, proj4string = CRS(proj4string(SIGEOtrees_SI17)))

#Create spdf with IDs (one unique ID per poly) and plot polygons
camview_spo.df_SI17<-SpatialPolygonsDataFrame(camview_spo_SI17,data.frame(id = unique(camview1_SI17$Deployment),row.names = unique(camview1_SI17$Deployment)))
plot(camview_spo.df_SI17, add = T)

#Cut out tree data from within polygons
clip_polys_SI17<-intersect(SIGEOsmall_SI17,camview_spo.df_SI17)
plot(clip_polys_SI17)
cvtrees_SI17<-as.data.frame(clip_polys_SI17)

#Pull and total the # of trees per deployment and change column names
cvtreecount_SI17<-cvtrees_SI17[,c(4,28)]
cvtreecount1_SI17<-aggregate(cvtreecount_SI17[,1], by = list(cvtreecount_SI17$d),sum)
colnames(cvtreecount1_SI17)[2]<-"Number_of_Trees"
colnames(cvtreecount1_SI17)[1]<-"Deployment"

#____________________________
#Variable - Oak Trees####
#____________________________
#Pull Oak Tree Data from grid
Oaks <-subset(treesSmall,sp %in% c('qual','quru','quco','qufa','qupr','quve','qusp','qumi'))
plot(Oaks, pch = 19)
#plot camera locations in red
plot(trapxySP, pch = 22, col = "red", add = T)

#add column to study site tree info that divides trees into 5 color size cateories
Oaks$Size_Category[Oaks$dbh <150]<-'461' #turqoise
Oaks$Size_Category[Oaks$dbh >150]<-'68' #dark blue
Oaks$Size_Category[Oaks$dbh >300]<-'47' #yellow
Oaks$Size_Category[Oaks$dbh >600]<-'139' #green
Oaks$Size_Category[Oaks$dbh >900]<-'8' #gray
Oaks$Size_Category[Oaks$dbh >1200]<-'550' #pink

#plot Oak Tree sizes by color
par(mar=c(5,17,4,2))
plot(Oaks,pch = 19, col = Oaks$Size_Category, add = T)

#Legend matching color to size (not working for some reason)
legend(747285,4309044, legend = c("< 15cm","> 15cm","> 30cm","> 60cm","> 90cm","> 120cm"), col = c("461", "68", "47","139", "8", "550"), pch = 19, title = "DBH of Oak Trees", bty = 'n')

#Cut out oak tree data from within the cones
library(rowr) #what is this for?
polyoaktrees_SI17<-intersect(Oak_Trees_SI17, camview_spo.df_SI17)
plot(polyoaktrees_SI17)
polyoaktreesdf_SI17<-as.data.frame(polyoaktrees_SI17)

#Pull # of oaks out of each deployment and rename columns to prepare for merge
oakcount_SI17<-polyoaktreesdf_SI17[,c(4,29)]
oakcount1_SI17<-aggregate(oakcount_SI17[,1],by=list(oakcount_SI17$d), sum)
colnames(oakcount1_SI17)[2]<-"Num_Oaks"
colnames(oakcount1_SI17)[1]<-"Deployment"

#Pull DBH of oaks from each deployment and add total
Oak_DBH_SI17<-polyoaktreesdf_SI17[,c(11,29)]
Oak_DBH1_SI17<-aggregate(Oak_DBH_SI17[,1],by=list(Oak_DBH_SI17$d), sum)
colnames(Oak_DBH1_SI17)[1]<-"Deployment"
colnames(Oak_DBH1_SI17)[2]<-"DBH"

#########################################
#Variable - Estimated Detection Distance
#########################################
#Bring in Detection Distance Data for Each Species
setwd("C:/Users/josey/Documents/CT Grid/Summer2017")
#Deer EDD Data
Deer_EDD_Summer2017<-read.csv("Deer_EDD_S17.csv")
Deer_EDD_Summer2017<-as.data.frame(Deer_EDD_Summer2017)

#Bear EDD Data
#Squirrel EDD Data
#Raccoon EDD Data

############################################################
#Deer Data Frame with XY coords
Deercamdata_SI17<-deercamdata_coords_SI17[-c(2,5)]

#Add Variables to Data Frame
Deercamdata_SI17$Cam_Nights<-camnightdata_SI17$Deploy.Duration
Deercamdata_SI17$Log<-camnightdata_SI17$Log.in.View
Deercamdata_SI17$EDD<-Deer_EDD_Summer2017$ESW.EDR
Deercamdata_SI17$Cam_Height<-Cam_heights_SI17$Camera_Height
Deercamdata_SI17$Num_Stems<-cvtreecount1_SI17$Number_of_Trees


#Merge because not all cameras have oak trees
Deercamdata_SI17<-merge(Deercamdata_SI17, Oak_DBH1_SI17, all = TRUE)
Deercamdata_SI17[is.na(Deercamdata_SI17)] <- 0


#######################################################
#Bear Data Frame with XY coords
Bearcamdata_SI17<-bearcamdata_coords_SI17[-c(2,5)]

#Add Variables to Data Frame
Bearcamdata_SI17$Cam_Nights<-camnightdata_SI17$Deploy.Duration
Bearcamdata_SI17$Log<-camnightdata_SI17$Log.in.View
#Bearcamdata_SI17$EDD<-Bear_EDD_Summer2017$ESW_EDR
Bearcamdata_SI17$Cam_Height<-Cam_heights_SI17$Camera_Height
Bearcamdata_SI17$Num_Stems<-cvtreecount1_SI17$Number_of_Trees

#Merge because not all cameras have oak trees
Bearcamdata_SI17<-merge(Bearcamdata_SI17, Oak_DBH1_SI17, all = TRUE)
Bearcamdata_SI17[is.na(Bearcamdata_SI17)] <- 0
########################################################
#Squirrel Data Frame with XY coords
Sqrlcamdata_SI17<-sqrlcamdata_coords_SI17[-c(2,5:19)]

#Add Variables to Data Frame
Sqrlcamdata_SI17$Cam_Nights<-camnightdata_SI17$Deploy.Duration
Sqrlcamdata_SI17$Log<-camnightdata_SI17$Log.in.View
#Sqrlcamdata_SI17$EDD<-Sqrl_EDD_Summer2017$ESW_EDR
Sqrlcamdata_SI17$Cam_Height<-Cam_heights_SI17$Camera_Height
Sqrlcamdata_SI17$Num_Stems<-cvtreecount1_SI17$Number_of_Trees

#Merge because not all cameras have oak trees
Sqrlcamdata_SI17<-merge(Sqrlcamdata_SI17, Oak_DBH1_SI17 , all = TRUE)
Sqrlcamdata_SI17[is.na(Sqrlcamdata_SI17)] <- 0

#########################################################
#Raccoon Data Frame with XY coords
Raccooncamdata_SI17<-raccooncamdata_coords_SI17[-c(2,5)]

#Add Variables to Data Frame
Raccooncamdata_SI17$Cam_Nights<-camnightdata_SI17$Deploy.Duration
Raccooncamdata_SI17$Log<-camnightdata_SI17$Log.in.View
#Raccooncamdata_SI17$EDD<-Raccoon_EDD_Summer2017$ESW_EDR
Raccooncamdata_SI17$Cam_Height<-Cam_heights_SI17$Camera_Height
Raccooncamdata_SI17$Num_Stems<-cvtreecount1_SI17$Number_of_Trees

#Merge because not all cameras have oak trees
Raccooncamdata_SI17<-merge(Raccooncamdata_SI17, Oak_DBH1_SI17, all = TRUE)
Raccooncamdata_SI17[is.na(Raccooncamdata_SI17)] <- 0

########################################################
#Save species data frames to bring in for Analysis later
########################################################
#Change working directory to general Camera Grid Folder
setwd("C:/Users/josey/Documents/CT Grid")
write.csv(Deercamdata_SI17, file = "Deer_DF.csv")
write.csv(Bearrcamdata_SI17, file = "Bear_DF.csv")
write.csv(Sqrlcamdata_SI17, file = "Sqrl_DF.csv")
write.csv(Raccooncamdata_SI17, file = "Raccoon_DF.csv")