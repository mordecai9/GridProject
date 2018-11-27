#Script for formatting and analysing camera trap data from high resolution
#camera grid established in Posey Hollow SCBI end of May 2017. 
#Twenty-seven Reconyx cameras approximately 20 meters apart.
#Cameras run in 2 month seasonal periods and distance markers at 2.5 meters apart up to 10 meters.
#Photos uploaded to eMammal and data downloaded from eMammal as csv file.
#Data currently stored at X:\1 eMammal\eMammal Projects\SCBICameraTestGrid
##
############################################################################
#Calculate Number of Species Sequences for Each Deployment (Camera Station)
############################################################################
library(dplyr)
library(reshape2)

#Bringing in emammal raw data from summer 2017 deployments
camdata_SI17<-read.csv("data/rawGridData_4S_final.csv", stringsAsFactors = FALSE)

#Convert Begin.Time and End.Time columns to date types
#replace "T" in time columns with a space first
camdata_SI17$Begin.Time<-sub(pattern = 'T',' ', camdata_SI17$Begin.Time)
camdata_SI17$End.Time<-sub(pattern = 'T', ' ', camdata_SI17$End.Time)
camdata_SI17$Begin.Time<-as.POSIXct(as.character(camdata_SI17$Begin.Time, "%Y-%b-%d %H:%M:%S", tz = "EST5EDT"))
camdata_SI17$End.Time<-as.POSIXct(as.character(camdata_SI17$End.Time, "%Y-%b-%d %H:%M:%S", tz = "EST5EDT"))


#Substitute Deployment named '0104-SI17-2' to '0104-SI17'
#This was the result of uploading this deployment on two seperate days to eMammal
camdata_SI17$Deployment.Name<-sub(pattern = '0104-SI17-2','0104-SI17',camdata_SI17$Deployment.Name)

#bringing in deployment length information
#Maps_CamChecks
camnightdata_SI17<-read.csv("data/CameraOperation_SI17.csv", stringsAsFactors = FALSE)

#Convert time in and out date columns to read as dates
#Replace "T" in the time columns with a space first
camnightdata_SI17$Date.out<-sub(pattern = 'T',' ', camnightdata_SI17$Date.out)
camnightdata_SI17$Date.in<-sub(pattern = 'T', ' ', camnightdata_SI17$Date.in)
camnightdata_SI17$Shift.begin.1<-sub(pattern = 'T', ' ', camnightdata_SI17$Shift.begin.1)
camnightdata_SI17$Shift.end.1<-sub(pattern = 'T', ' ', camnightdata_SI17$Shift.end.1)
camnightdata_SI17$Shift.begin.2<-sub(pattern = 'T', ' ', camnightdata_SI17$Shift.begin.2)
camnightdata_SI17$Shift.end.2<-sub(pattern = 'T', ' ', camnightdata_SI17$Shift.end.2)
camnightdata_SI17$Date.out<-as.POSIXct(as.character(camnightdata_SI17$Date.out,"%Y-%b-%d %H:%M:%S", tz = "EST5EDT"))
camnightdata_SI17$Date.in<-as.POSIXct(as.character(camnightdata_SI17$Date.in,"%Y-%b-%d %H:%M:%S", tz = "EST5EDT"))
camnightdata_SI17$Shift.begin.1<-as.POSIXct(as.character(camnightdata_SI17$Shift.begin.1, "%Y-%b-%d %H:%M:%S", tz = "EST5EDT"))
camnightdata_SI17$Shift.end.1<-as.POSIXct(as.character(camnightdata_SI17$Shift.end.1, "%Y-%b-%d %H:%M:%S", tz = "EST5EDT"))
camnightdata_SI17$Shift.begin.2<-as.POSIXct(as.character(camnightdata_SI17$Shift.begin.2, "%Y-%b-%d %H:%M:%S", tz = "EST5EDT"))
camnightdata_SI17$Shift.end.2<-as.POSIXct(as.character(camnightdata_SI17$Shift.end.2, "%Y-%b-%d %H:%M:%S", tz = "EST5EDT"))

#Calculate Deployment Duration by subtracting shift length from deployment length
#Shift means a period of days where the camera was shifted and is unusable data
camnightdata_SI17$Deploy.Length<-difftime(camnightdata_SI17$Date.in ,camnightdata_SI17$Date.out , units = c("days"))
camnightdata_SI17$Deploy.Shifts.1<-difftime(camnightdata_SI17$Shift.end.1 ,camnightdata_SI17$Shift.begin.1 , units = c("days"))
camnightdata_SI17$Deploy.Shifts.2<-difftime(camnightdata_SI17$Shift.end.2 ,camnightdata_SI17$Shift.begin.2 , units = c("days"))
camnightdata_SI17$Deploy.Duration<-camnightdata_SI17$Deploy.Length-camnightdata_SI17$Deploy.Shifts.1-camnightdata_SI17$Deploy.Shifts.2

#Substitute periods in column titles to underscores (if needed)
names(camdata_SI17) <- sub(pattern = 'Deployment.', 'Deployment_',names(camdata_SI17))

#Merge camdata and number of camera nights by deployment name
#Make sure all "Scent" deployments from camnightdata are gone after merge
camdata_summary_SI17<-merge(camnightdata_SI17, camdata_SI17, by = "Deployment_Name")

#remove rows that occur after end dates,before dates, and shift dates for deployments
camdata_summaryd_SI17<-subset(camdata_summary_SI17, Begin.Time>=Date.out&Begin.Time<=Date.in)
camdata_summaryd_SI17<-subset(camdata_summaryd_SI17, End.Time<Shift.begin.1|Begin.Time>=Shift.end.1)
camdata_summaryd_SI17<-subset(camdata_summaryd_SI17, End.Time<Shift.begin.2|Begin.Time>=Shift.end.2)

setwd("C:/Users/josey/Documents/CT Grid")
write.csv(camdata_summaryd_SI17, file = "Camdata_SI17.csv")

#########################################################################
#Calculate Average, Minimum, and Maximum Number of Sequences Across Grid
#########################################################################
#create summary table of frequencies of sequences and convert it to a data frame
camdata_summarym_SI17<-table(camdata_summaryd_SI17$Deployment_Name, 
                        camdata_summaryd_SI17$Species.Name)
camdata_summarym_SI17<-as.data.frame(camdata_summarym_SI17)

#change column names
colnames(camdata_summarym_SI17)<-c("Deployment_Name", "Species", "Number_of_Sequences")


#calculating average number of sequences for each species across the grid
#This can be seen as a baseline sequence number against which we can compare sequence numbers at individual cameras
baselines_SI17<- tapply(camdata_summarym_SI17$Number_of_Sequences, camdata_summarym_SI17$Species, mean)
baselines_SI17<- as.data.frame(baselines_SI17)
names(baselines_SI17)[1] <- "Mean_Seq"

#SD of average number of sequences across the grid for each species
baselines_sd_SI17 <- tapply(camdata_summarym_SI17$Number_of_Sequences, camdata_summarym_SI17$Species, sd)
baselines_SI17$sd<- baselines_sd_SI17
names(baselines_SI17)[2] <- "StDev_Seq"

#number of sequences by species by camera, then transposing to have species as rows
CR_bD_bS_SI17 <- tapply(camdata_summarym_SI17$Number_of_Sequences, list(camdata_summarym_SI17$Species, camdata_summarym_SI17$Deployment_Name), mean)
CR_bD_bS_SI17 <- as.data.frame(CR_bD_bS_SI17)

#Min and max of number of sequences for each species across the grid
baselines_SI17$min <- apply(CR_bD_bS_SI17,1,min)
names(baselines_SI17)[3]<- "Min_Seq"
baselines_SI17$max <- apply(CR_bD_bS_SI17,1,max)
names(baselines_SI17)[4]<-"Max_Seq"

#Add columns to CR_bD_bS for standard deviation, mean, max and min.
CR_bD_bS_SI17$Mean_Seq <- baselines_SI17$Mean_Seq
CR_bD_bS_SI17$StDev_Seq<-baselines_SI17$StDev_Seq
CR_bD_bS_SI17$Min_Seq<-baselines_SI17$Min_Seq
CR_bD_bS_SI17$Max_Seq<-baselines_SI17$Max_Seq

#remove humans,birds,and misfires
CR_bD_bS1_SI17<-CR_bD_bS_SI17[-c(1,2,7,9,11,16:21), ]


#############################################################
#Calculate Proportion of Species Captures for Each Camera
#############################################################
#Create replicate data frame and remove mean,sd,min,max columns
CR_bD_bS2_SI17<-CR_bD_bS1_SI17
CR_bD_bS2_SI17<-CR_bD_bS2_SI17[,-c(28:31)]

#convert all values >0 to 1
CR_bD_bS3_SI17<-as.data.frame((ifelse(CR_bD_bS2_SI17==0,0,1)))

#create column of camera capture proportion per species for whole grid to original data frame
CR_bD_bS1_SI17$Prop_Grid<-rowSums(CR_bD_bS3_SI17)/length(unique(camdata_SI17$Deployment_Name))

##############################################################################################
#Create Boxplot of Number of Sequences and Proportion of sequences per Species per Camera
##############################################################################################

#Transpose headers
CR_SI17<-t(CR_bD_bS2_SI17)

#reshape data so that each row is per deployment, per species, per capture rate
CR1_SI17<-melt(CR_SI17)
colnames(CR1_SI17)<-c ("Deployment", "Species", "Number_of_Sequences")

#Calculates proportion of cameras that captured each species
CR2_SI17<-subset(CR1_SI17, Number_of_Sequences !=0)
prop_per_cam_SI17<-as.data.frame(table(CR2_SI17$Species))
colnames(prop_per_cam_SI17)<-c("Species","Prop")
####################################################
#Edits stop here
#########################################
#Merge proportion dataframe with rest of data
CR3_SI17<-merge(CR1_SI17, prop_per_cam_SI17, by = "Species")

CR3_SI17["tot"]<-length(unique(CR1_SI17$Deployment))

#Adds column label which outputs the desired text labels for the boxplot
CR3_SI17["label"]<-paste(CR3_SI17$Species,"|", CR3_SI17$Prop, "/", CR3_SI17$tot)

#Boxplot of the sequence number per species
par(mar=c(9,17,4,2))
plot<-boxplot(CR3_SI17$Number_of_Sequences~CR3_SI17$label,
              cex.main = 2.5,
              main = "Summer 2017",
              cex.lab = 1.7,
              cex.axis = 1.3,
              horizontal = T, las = 2, cex.axis = 1,
              names.arg = CR3_SI17$label, at=rank(tapply(CR3_SI17$Number_of_Sequences,CR3_SI17$Species,mean)))
mtext(expression(bold("Species")), side = 2, line = 15, cex = 1.7)
mtext(expression(bold("Number_of_Sequences (per 100 camera-nights)")), side = 1 , line = 4, cex = 1.7)

################################################################################################################
#Create Individual Species Data frames with Deployment, #Seq, Detection Distance, Cam Height, #Trees, duration
##############################################################################################################
#Subset number of sequences by species and deployment
CR3_Sequences_SI17<-CR3_SI17[,1:3]

#camdata_summaryd_SI17$Caprate <- (deercamdata_coords$Number_of_Sequences/camdata_summaryd_SI17$Deploy.Duration) *100

#Bring in csv with correct coordinates
setwd("C:/Users/josey/Documents/CT Grid/Summer2017")
camdata_coordinates_SI17 <- read.csv("Grid_Coordinates_SI17.csv")

#Merge capture rates per species with coordinates of Deployments by the Deployment column
camdata_coords_SI17<-merge(camdata_coordinates_SI17, CR3_Sequences_SI17,  by = "Deployment")

#create individual species data frames with species name, capture rate, and coordinates
deercamdata_coords_SI17<-subset(camdata_coords_SI17, Species == "Odocoileus virginianus")
bearcamdata_coords_SI17<-subset(camdata_coords_SI17, Species == "Ursus americanus")
coyotecamdata_coords_SI17<-subset(camdata_coords_SI17, Species == "Canis latrans")
foxsqrlcamdata_coords_SI17<-subset(camdata_coords_SI17, Species == "Sciurus niger")
redfoxcamdata_coords_SI17<-subset(camdata_coords_SI17, Species == "Vulpes vulpes")
bobcatcamdata_coords_SI17<-subset(camdata_coords_SI17, Species == "Lynx rufus")
grsqrlcamdata_coords_SI17<-subset(camdata_coords_SI17, Species == "Sciurus carolinensis")
raccooncamdata_coords_SI17<-subset(camdata_coords_SI17, Species == "Procyon lotor")
opossumcamdata_coords_SI17<-subset(camdata_coords_SI17, Species =="Didelphis virginiana")
unknownsqrlcamdata_coords_SI17<-subset(camdata_coords_SI17, Species == "Unknown Squirrel")

###Find Capture Rate###
#deercamdata_coords$Deploy.Duration<-camnightdata_SI17$Deploy.Duration
#deercamdata_coords$Deploy.Duration<-as.numeric(deercamdata_coords$Deploy.Duration)
#deercamdata_coords$Capture_Rate <- (deercamdata_coords$Number_of_Sequences/deercamdata_coords$Deploy.Duration) *100


#Merge all squirrel species
sqrlcamdata_coords_SI17<-merge(foxsqrlcamdata_coords_SI17, grsqrlcamdata_coords_SI17, by = "Deployment")
sqrlcamdata_coords_SI17<-merge(sqrlcamdata_coords_SI17, unknownsqrlcamdata_coords_SI17, by = "Deployment")
sqrlcamdata_coords_SI17$Number_of_Sequences_Tot<-sqrlcamdata_coords_SI17$Number_of_Sequences.x + 
  sqrlcamdata_coords_SI17$Number_of_Sequences.y + sqrlcamdata_coords_SI17$Number_of_Sequences


##################################
#Bring in the different variables
#Variable - Camera Height
################################
#Bring in Camera Height Data
Cam_heights_SI17<-read.csv("Camera_Heights_SI17.csv")
Cam_heights_SI17<-as.data.frame(Cam_heights_SI17)

#Rename columns to match for merge
names(Cam_heights_SI17)[1]<- "Deployment"
names(Cam_heights_SI17)[3]<- "Camera_Height"

#CAMERA HEIGHT is the distance in centimeters from the ground to the camera lens#

#Check that 'Camera_Height' is an integer variable
str(Cam_heights_SI17)

############################################
#Variable - Number of Trees in camera sight
###########################################
#Working with SIGEO tree grid information to try to associate with capture rates from our
#high resolution camera grid. Grid established summer 2017, running through summer 2018.
#Coordinates should be UTM Zone 17S

library(rgeos)
library(rgdal)
library(sp)
library(maptools)
library(raster)
library(grid)

#Bring in geo-reference tree data from entire SIGEO grid
setwd("C:/Users/josey/Documents/CT Grid")
list.files()
SIGEOtrees_SI17<-read.csv("scbi.full2_UTM_lat_long_12012017.csv")

#Change data frame into a Spatial Points Data frame
head(SIGEOtrees_SI17)
coordinates(SIGEOtrees_SI17)<- c("NAD83_X", "NAD83_Y")
class(SIGEOtrees_SI17)
#plot(SIGEOtrees)

#plot the coordinates
plot(camdata_coordinates_SI17$NAD83_X,
     camdata_coordinates_SI17$NAD83_Y,
     xlim = c(747420, 747560),
     ylim = c(4308900,4309040))

#Convert this trap coordinate information into a spatialpoints object
#First need to have the xy coordinates as a separate matrix
trapxy_SI17 <- camdata_coordinates_SI17[, c(2,3)]
trapxy_sp_SI17 <- SpatialPointsDataFrame(coords = trapxy_SI17, data = camdata_coordinates_SI17,
                                    proj4string = CRS(proj4string(SIGEOtrees_SI17)))
plot(trapxy_SI17)

#Create a clipping polygon to reduce the size of the SIGEO grid to just the area of interest
#I'm setting the extent as 50m around the extreme trap coordinates
c<-50
CP_SI17 <- as(extent(min(trapxy_SI17$NAD83_X)-c, 
                max(trapxy_SI17$NAD83_X)+c,
                min(trapxy_SI17$NAD83_Y)-c,
                max(trapxy_SI17$NAD83_Y)+c),
         "SpatialPolygons")

#Assign the coordinate reference system of SIGEOtrees to the new clipping polygon         
proj4string(CP_SI17) <- CRS(proj4string(SIGEOtrees_SI17))
plot(CP_SI17)

#You could also use gIntersect below but it does not preserve the original attribute data
SIGEOsmall_SI17 <- intersect(SIGEOtrees_SI17, CP_SI17)

#plot grid with tree and cameras
plot(SIGEOsmall_SI17, col = "darkgreen", pch = 3,cex.main = 4)
plot(trapxy_sp_SI17, pch = 19, col = "red", add = T)

#Add a legend
par(font = 2)
legend(747300,4308970, legend = c("Tree", "Camera"), col = c("darkgreen", "red"), 
       pch = c(3,19), cex =1.5, bty = "n")

#Add scale
scale.len <- 20
x <- c(747308.5,747308.5+scale.len)
y<- c(4308890, 4308890)
lines(x,y,lwd = 2)
text(747347.9, 4308890, '20m', cex = 1.5)

#Add Deployment label to each camera
#pointLabel(coordinates(trapxy_sp),labels=trapxy_sp@data$Deployment, cex = 0.7, allowSmallOverlap = T)

#########################################################
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
camview1_SI17[28:54,]<-(camview_SI17[1:27, c(4:5,3)])
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

######################
#Variable - Oak Trees
######################
#Pull Oak Tree Data from grid
Oak_Trees_SI17<-subset(SIGEOsmall_SI17,sp %in% c('qual','quru','quco','qufa','qupr','quve','qusp','qumi'))
plot(Oak_Trees_SI17, pch = 19)
#plot camera locations in red
plot(trapxy_sp_SI17, pch = 22, col = "red", add = T)

#add column to study site tree info that divides trees into 5 color size cateories
Oak_Trees_SI17$Size_Category[Oak_Trees_SI17$dbh <150]<-'461' #turqoise
Oak_Trees_SI17$Size_Category[Oak_Trees_SI17$dbh >150]<-'68' #dark blue
Oak_Trees_SI17$Size_Category[Oak_Trees_SI17$dbh >300]<-'47' #yellow
Oak_Trees_SI17$Size_Category[Oak_Trees_SI17$dbh >600]<-'139' #green
Oak_Trees_SI17$Size_Category[Oak_Trees_SI17$dbh >900]<-'8' #gray
Oak_Trees_SI17$Size_Category[Oak_Trees_SI17$dbh >1200]<-'550' #pink

#plot Oak Tree sizes by color
par(mar=c(5,17,4,2))
plot(Oak_Trees_SI17,pch = 19, col = Oak_Trees_SI17$Size_Category, add = T)

#Legend matching color to size 
legend(747285,4309044, legend = c("< 15cm","> 15cm","> 30cm","> 60cm","> 90cm","> 120cm"), col = c("461", "68", "47","139", "8", "550"), pch = 19, title = "DBH of Oak Trees", bty = 'n')

#Cut out oak tree data from within the cones
library(rowr)
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
