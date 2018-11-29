

###############################################################################
#Regression Analysis of Deer Capture Rate by Deer Estimated Detection Distance
###############################################################################
#create data frame from Deer EDD csv
setwd("C:/Users/josey/Documents/CT Grid/Summer2017")
Deer_EDD_Summer2017<-read.csv("Deer_EDD_S17.csv")
Deer_EDD_Summer2017<-as.data.frame(Deer_EDD_Summer2017)

CR4<-subset(CR3, Species =="Odocoileus virginianus")

#Rename columns to match camdata data frame
colnames(Deer_EDD_Summer2017)[6]<-"Deployment"
colnames(Deer_EDD_Summer2017)[19]<-"Number_of_Detections"
colnames(Deer_EDD_Summer2017)[14]<-"ESW_EDR"
Deer_EDD_Summer2017$CapRate<-deercamdata_coords$Capture_Rate

#Remove deployments with <20 detections
Deer_EDD<- Deer_EDD_Summer2017[which(Deer_EDD_Summer2017$CapRate >= 20),]

#Merge deer caprate and EDD by Deployment and Species
Deer_caprate_EDD<-merge(Deer_EDD, CR4, by = "Deployment")

#boxplot of merged data to identify outliers
par(mar=c(5,5,4,2))
boxplot(Deer_caprate_EDD$ESW_EDR)
boxplot.stats(Deer_caprate_EDD$ESW_EDR)$out

#remove outliers
Deer_caprate_EDD1<-Deer_caprate_EDD[which(Deer_caprate_EDD$ESW_EDR <=10),]

#correlate capture rate with EDD and show regression line
par(mar=c(5,6,4,2))
with(Deer_caprate_EDD1, plot(CapRate ~ ESW_EDR,
                             pch = 19,
                             xlab = expression(bold("Effective Detection Distance (m)")),
                             ylab = expression(bold("Deer Capture Rate")),
                             main = "Regression Analysis of Deer Capture Rate and Effective Detection Distance",
                             cex.main = 2,
                             cex.axis = 1.3,
                             cex.lab = 1.7))
lm.outdedd = lm(CapRate ~ ESW_EDR, data = Deer_caprate_EDD1)
abline(lm.outdedd, col="blue")
summary(lm.outdedd)

#Add Rsquared value expression
#Need to create an object with just rsquared value first
Rsquared<-summary(lm.outdedd)$r.squared
text(5.866514,191.3044, as.expression(substitute(italic(R)^2 == r,list(r=round(Rsquared,3)))), cex = 1.5, font = 2)


#########################################################################################################
#Regression Analysis of Deer Estimated Detection Distance with Bear, Gray Squirrel, and Racoon Capture Rate
#########################################################################################################
#subset capture rate per deployment of bear
bearcam_caprates<-bearcamdata_coords[,c(1,6,7)]

#rename column title
names(bearcam_caprates)[3]<-"Bear_Capture_Rate"

#Merge Deer EDD and bear capture rate info by Deployment
bear_dEDD<-merge(Deer_EDD, bearcam_caprates, by = "Deployment")

#plot results
boxplot(bear_dEDD$ESW_EDR, main = "Effective Detection Distance for Black Bear",
        cex.main = 1.7)
boxplot.stats(bear_dEDD$ESW_EDR)$out

#remove outliers
bear_dEDD1<-bear_dEDD[which(bear_dEDD$ESW_EDR <=10),]

#plot regression analysis of Bear Capture Rate by Deer EDD without outliers
par(mar=c(6,6,4,6))
with(bear_dEDD1, plot(Bear_Capture_Rate ~ ESW_EDR, main = "Regression Analysis of Bear Capture Rate and Deer Effective Detection Distance",
                      cex.main = 2.2,
                      xlab = expression(bold("Deer Effective Detection Distance")),
                      ylab = expression(bold("Bear Capture Rate")),
                      cex.axis = 1.3,
                      cex.lab = 1.5))
lm.out_db = lm(Bear_Capture_Rate ~ ESW_EDR, data = bear_dEDD1)
abline(lm.out_db, col="blue")
summary(lm.out_db)

##################################################
#subset gray squirrel capture rate per deployment
grsqrlcam_caprates<-grsqrlcamdata_coords[,c(1,6,7)]

#rename column title
names(grsqrlcam_caprates)[3]<-"Gray_Squirrel_Capture_Rate"

#Merge gray squirrel info and Deer EDD
grsqrl_dEDD<-merge(Deer_EDD, grsqrlcam_caprates, by = "Deployment")

#plot results
boxplot(grsqrl_dEDD$ESW_EDR)
boxplot.stats(grsqrl_dEDD$ESW_EDR)$out

#remove outliers
grsqrl_dEDD1<-grsqrl_dEDD[which(grsqrl_dEDD$ESW_EDR <=10),]

#plot regression analysis of Gray Squirrel Capture Rate by Deer EDD without outliers
par(mar=c(11,6,4,6))
with(grsqrl_dEDD1, plot(Gray_Squirrel_Capture_Rate ~ ESW_EDR, main = "Regression Analysis of Gray Squirrel Capture Rate and Deer Effective Detection Distance",
                        cex.main = 1.9,
                        xlab = expression(bold("Deer Effective Detection Distance")),
                        ylab = expression(bold("Gray Squirrel Capture Rate")),
                        cex.axis = 1,
                        cex.lab = 1.5))
lm.out_gd = lm(Gray_Squirrel_Capture_Rate ~ ESW_EDR, data = grsqrl_dEDD1)
abline(lm.out_gd, col="blue")
summary(lm.out_gd)

############################################
#subset raccoon capture rates per deployment
raccooncam_caprates<-raccooncamdata_coords [,c(1,6,7)]

#rename column
names(raccooncam_caprates)[3]<-"Raccoon_Capture_Rate"

#merge raccoon capture rate and Deer EDD by Deployment
raccoon_dEDD<-merge(Deer_EDD, raccooncam_caprates, by = "Deployment")

#plot results
boxplot(raccoon_dEDD$ESW_EDR)
boxplot.stats(raccoon_dEDD$ESW_EDR)$out

#remove outliers
raccoon_dEDD1<-raccoon_dEDD[which(raccoon_dEDD$ESW_EDR <=10),]

#plot regression analysis of Raccoon Capture Rate by Deer EDD without outliers
par(mar=c(11,6,4,6))
with(raccoon_dEDD1, plot(Raccoon_Capture_Rate ~ ESW_EDR, main = "Regression Analysis of Raccoon Capture Rate and Deer Effective Detection Distance",
                         cex.main = 1.9,
                         xlab = expression(bold("Deer Effective Detection Distance")),
                         ylab = expression(bold("Raccoon Capture Rate")),
                         cex.axis = 1,
                         cex.lab = 1.5))
lm.out_rd = lm(Raccoon_Capture_Rate ~ ESW_EDR, data = raccoon_dEDD1)
abline(lm.out_rd, col="blue")
summary(lm.out_rd)

##################################################################################################
#Regression Analysis of the Deer Capture Rate to Capture Rate of Bear, Gray Squirrel, and Raccoon
#################################################################################################
# DEER AND BEAR CAPTURE RATE
#Subset Species, Deployment and caprate from deer and bear data frames
deercam_caprate<-deercamdata_coords[,c(1,6,7)]
names(deercam_caprate)[3]<-"Deer_Capture_Rate"
bearcam_caprate<-bearcamdata_coords[,c(1,6,7)]
names(bearcam_caprate)[3]<-"Bear_Capture_Rate"
deerbearcam<-merge(bearcam_caprate, deercam_caprate, by = "Deployment")

#boxplot to identify outliers
boxplot(deercam_caprate$Deer_Capture_Rate,bearcam_caprate$Bear_Capture_Rate)
boxplot.stats(deercam_caprate$Deer_Capture_Rate)$out
boxplot.stats(bearcam_caprate$Bear_Capture_Rate)$out

#remove outliers
deercam_caprate1<-deercam_caprate[which(deercam_caprate$Deer_Capture_Rate <=166),]
bearcam_caprate1<-bearcam_caprate[which(bearcam_caprate$Bear_Capture_Rate <=17),]

#Plot Regression Analysis of capture rate of deer and bear
deerbearcam<-merge(bearcam_caprate1, deercam_caprate1, by = "Deployment")
with(deerbearcam, plot(Deer_Capture_Rate ~ Bear_Capture_Rate, main = " Capture Rate Correlation Between Deer and Black Bear",
                       cex.main = 2.7,
                       font.lab = 2.3,
                       cex.axis = 1.5,
                       cex.lab = 1.7))
lm.outbd = lm(Deer_Capture_Rate ~ Bear_Capture_Rate, data = deerbearcam)
abline(lm.outbd, col="blue")
summary(lm.outbd)

#################################################
#DEER AND GRAY SQUIRREL CAPTURE RATE
#Subset Species, Deployment and caprate from individual species data frames
deercam_caprate<-deercamdata_coords[,c(1,6,7)]
names(deercam_caprate)[3]<-"Deer_Capture_Rate"
grsqrlcam_caprate<-grsqrlcamdata_coords[,c(1,6,7)]
names(grsqrlcam_caprate)[3]<-"GrSqrl_Capture_Rate"
deergrsqrlcam<-merge(grsqrlcam_caprate, deercam_caprate, by = "Deployment")

#boxplot to identify outliers
boxplot(deercam_caprate$Deer_Capture_Rate,grsqrlcam_caprate$GrSqrl_Capture_Rate)
boxplot.stats(deercam_caprate$Deer_Capture_Rate)$out
boxplot.stats(grsqrlcam_caprate$GrSqrl_Capture_Rate )$out

#remove outliers
deercam_caprate<-deercam_caprate[which(deercam_caprate$Deer_Capture_Rate <=166),]
grsqrlcam_caprate<-grsqrlcam_caprate[which(grsqrlcam_caprate$GrSqrl_Capture_Rate <=53),]

#Plot Regression Analysis of capture rate of Deer and Gray Squirrel
deersqrlcam<-merge(grsqrlcam_caprate, deercam_caprate, by = "Deployment")
with(deersqrlcam, plot(Deer_Capture_Rate ~ GrSqrl_Capture_Rate, main = " Capture Rate Correlation Between Deer and Gray Squirrel",
                       cex.main = 2.6,
                       font.lab = 2.3,
                       cex.axis = 1.5,
                       cex.lab = 1.7))
lm.outgd = lm(Deer_Capture_Rate ~ GrSqrl_Capture_Rate, data = deersqrlcam)
abline(lm.outgd, col="blue")
summary(lm.outgd)

####################
#DEER AND RACOON CAPTURE RATE
#Subset Species, Deployment and caprate from individual species data frames
deercam_caprate<-deercamdata_coords[,c(1,6,7)]
names(deercam_caprate)[3]<-"Deer_Capture_Rate"
raccooncam_caprate<-raccooncamdata_coords[,c(1,6,7)]
names(raccooncam_caprate)[3]<-"Raccoon_Capture_Rate"
deerraccooncam<-merge(raccooncam_caprate, deercam_caprate, by = "Deployment")

#boxplot to identify outliers
boxplot(deercam_caprate$Deer_Capture_Rate,raccooncam_caprate$Raccoon_Capture_Rate)
boxplot.stats(deercam_caprate$Deer_Capture_Rate)$out
boxplot.stats(raccooncam_caprate$Raccoon_Capture_Rate)$out

#remove outliers
deercam_caprate<-deercam_caprate[which(deercam_caprate$Deer_Capture_Rate <=166),]
raccooncam_caprate<-raccooncam_caprate[which(raccooncam_caprate$Raccoon_Capture_Rate <=12),]

#Plot Regression Analysis of capture rate of Deer and Raccoon
deerraccooncam<-merge(raccooncam_caprate, deercam_caprate, by = "Deployment")
with(deerraccooncam, plot(Deer_Capture_Rate ~ Raccoon_Capture_Rate, main = " Capture Rate Correlation Between Deer and Raccoon",
                          cex.main = 2.6,
                          font.lab = 2.3,
                          cex.axis = 1.5,
                          cex.lab = 1.7))
lm.outrd = lm(Deer_Capture_Rate ~ Raccoon_Capture_Rate, data = deerraccooncam)
abline(lm.outrd, col="blue")
summary(lm.outrd)

#############################################################################################
#Correlation between Camera Height and Capture Rate of Deer, Bear, Gray Squirrel, and Raccoon
#############################################################################################
#Bring in camera height data
Cam_heights<-read.csv("Camera_Heights.csv")
Cam_heights<-as.data.frame(Cam_heights)

#Rename columns to match for merge
names(Cam_heights)[1]<- "Deployment"
names(Cam_heights)[3]<- "Camera_Height"

#CAMERA HEIGHT is the distance in centimeters from the ground to the camera lens#

#Check variable types in data frame and change the camera height to a number variable
str(Cam_heights)


#Merge the species' capture rate per deployment and each deployment's camera height, by "Deployment" column
Deer_Cr_Ch<-merge(deercamdata_coords, Cam_heights,  by = "Deployment")
Grsqrl_Cr_Ch<-merge(grsqrlcamdata_coords, Cam_heights, by = "Deployment")
Bear_Cr_Ch<-merge(bearcamdata_coords, Cam_heights, by = "Deployment")
Raccoon_Cr_Ch<-merge(raccooncamdata_coords, Cam_heights, by = "Deployment")

############################################################
#Regression Analysis of Deer capture rate by camera's height
par(mar=c(5,5,4,2))
with(Deer_Cr_Ch, plot(Capture_Rate ~ Camera_Height , main = "Regression Analysis of Deer Capture Rate and Camera Height",
                      cex.main = 2.2,
                      xlab = expression(bold("Camera Height")),
                      ylab = expression(bold("Capture Rate")),
                      cex.axis = 1.3,
                      cex.lab = 1.6))
lm.outdcrch = lm(Capture_Rate ~ Camera_Height, data = Deer_Cr_Ch)
abline(lm.outdcrch, col="blue")
summary(lm.outdcrch)

###################################################################
#Regression Analysis of Black Bear capture rate by camera's height
par(mar=c(5,5,4,2))
with(Bear_Cr_Ch, plot(Capture_Rate ~ Camera_Height , main = "Regression Analysis of Bear Capture Rate and Camera Height",
                      cex.main = 2.2,
                      xlab = expression(bold("Camera Height")),
                      ylab = expression(bold("Capture Rate")),
                      cex.axis = 1.3,
                      cex.lab = 1.6))
lm.outbcrch = lm(Capture_Rate ~ Camera_Height, data = Bear_Cr_Ch)
abline(lm.outbcrch, col="blue")
summary(lm.outbcrch)

#######################################################################
#Regression Analysis of Gray Squirrel capture rate by camera's height
par(mar=c(5,5,4,2))
with(Grsqrl_Cr_Ch, plot(Capture_Rate ~ Camera_Height , main = "Regression Analysis of Gray Squirrel Capture Rate and Camera Height",
                        cex.main = 2.2,
                        xlab = expression(bold("Camera Height")),
                        ylab = expression(bold("Capture Rate")),
                        cex.axis = 1.3,
                        cex.lab = 1.6))
lm.outgscrch = lm(Capture_Rate ~ Camera_Height, data = Grsqrl_Cr_Ch)
abline(lm.outgscrch, col="blue")
summary(lm.outgscrch)

###################################
#Boxplot Raccoon capture rate and remove outliers
boxplot(Raccoon_Cr_Ch$Capture_Rate)
boxplot.stats(Raccoon_Cr_Ch$Capture_Rate)$out
Raccoon_Cr_Ch<-Raccoon_Cr_Ch[which(Raccoon_Cr_Ch$Capture_Rate <=12),]

#Regression Analysis of Raccoon capture rate by camera's height
par(mar=c(5,5,4,2))
with(Raccoon_Cr_Ch, plot(Capture_Rate ~ Camera_Height , main = "Regression Analysis of Raccoon Capture Rate and Camera Height",
                         cex.main = 2.2,
                         xlab = expression(bold("Camera Height")),
                         ylab = expression(bold("Capture Rate")),
                         cex.axis = 1.3,
                         cex.lab = 1.6))
lm.outrcrch = lm(Capture_Rate ~ Camera_Height, data = Raccoon_Cr_Ch)
abline(lm.outrcrch, col="blue")
summary(lm.outrcrch)

################################################
#Correlation Between Camera Height and Deer EDD
################################################
#Merge Deer EDD with Camera Heights, by Deployment Column.
#Deer EDD data frame does not include deployments with <20 detections of deer
DEDD_CM<-merge(Deer_EDD, Cam_heights, by = "Deployment")

#Regression Analysis of Camera height by Deer EDD
boxplot(DEDD_CM$ESW_EDR)
boxplot.stats(DEDD_CM$ESW_EDR)$out
Raccoon_Cr_Ch<-Raccoon_Cr_Ch[which(Raccoon_Cr_Ch$Capture_Rate <=10),]

par(mar=c(5,5,4,2))
with(DEDD_CM, plot(ESW_EDR ~ Camera_Height , main = "Regression Analysis of Camera Height to Deer EDD",
                   cex.main = 2.2,
                   xlab = expression(bold("Camera Height")),
                   ylab = expression(bold("Deer EDD")),
                   cex.axis = 1.3,
                   cex.lab = 1.6))
lm.outcmedd = lm(ESW_EDR ~ Camera_Height, data = DEDD_CM)
abline(lm.outcmedd, col="blue")
summary(lm.outcmedd)

##############################################################################################
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
SIGEOtrees<-read.csv("scbi.full2_UTM_lat_long_12012017.csv")

#Change data frame into a Spatial Points Data frame
head(SIGEOtrees)
coordinates(SIGEOtrees)<- c("NAD83_X", "NAD83_Y")
class(SIGEOtrees)
#plot(SIGEOtrees)

#Bring in csv with correct coordinates
setwd("C:/Users/josey/Documents/CT Grid/Summer2017")
camdata_coordinates <- read.csv("Grid_Coordinates.csv")
camdata_coordinates <- as.data.frame(camdata_coordinates)

#plot the coordinates
plot(camdata_coordinates$NAD83_X,
     camdata_coordinates$NAD83_Y,
     xlim = c(747420, 747560),
     ylim = c(4308900,4309040))

#Convert this trap coordinate information into a spatialpoints object
#First need to have the xy coordinates as a separate matrix
trapxy <- camdata_coordinates[, c(2,3)]
trapxy_sp <- SpatialPointsDataFrame(coords = trapxy, data = camdata_coordinates,
                                    proj4string = CRS(proj4string(SIGEOtrees)))
plot(trapxy)

#Create a clipping polygon to reduce the size of the SIGEO grid to just the area of interest
#I'm setting the extent as 50m around the extreme trap coordinates
c<-50
CP <- as(extent(min(trapxy$NAD83_X)-c, 
                max(trapxy$NAD83_X)+c,
                min(trapxy$NAD83_Y)-c,
                max(trapxy$NAD83_Y)+c),
         "SpatialPolygons")

#Assign the coordinate reference system of SIGEOtrees to the new clipping polygon         
proj4string(CP) <- CRS(proj4string(SIGEOtrees))
plot(CP)

#You could also use gIntersect below but it does not preserve the original attribute data
SIGEOsmall <- intersect(SIGEOtrees, CP)

#plot grid with tree and cameras
plot(SIGEOsmall, col = "darkgreen", pch = 3,cex.main = 4)
plot(trapxy_sp, pch = 19, col = "red", add = T)

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

##################################################################
#10m Buffer zones and related code
##################################################################
#buffer each camera by 10m. Maybe use actual max detection distance for each camera instead?
cams10m <- gBuffer(trapxy_sp, width=10, byid=TRUE, quadsegs = 4)
plot(cams10m, add = T)

#Cut out tree data from within the 10m buffers
trees10m <- intersect(SIGEOtrees, cams10m)
plot(trees10m)
gridtrees10m<-as.data.frame(trees10m)

#Check if trees are listed twice in buffer zone overlaps
doubletree<-gridtrees10m[,c(2,32)]

#Pull and total the # of trees per deployment and change column names 
treecount<-gridtrees10m[,c(4,32)]
treecount1<-aggregate(treecount[,1],by=list(treecount$Deployment), sum)
colnames(treecount1)[2]<-"Number_of_Trees"
colnames(treecount1)[1]<-"Deployment"

#Merge number of trees per deployment with deer capture rate per deployment by Deployment
Trees_Per_Dcr<-merge(deercamdata_coords, treecount1,  by = "Deployment")

#Boxplot of Number of Trees and remove outliers
boxplot(Trees_Per_Dcr$Number_of_Trees)
boxplot.stats(Trees_Per_Dcr$Number_of_Trees)$out
Trees_pDCR<-Trees_Per_Dcr[which(Trees_Per_Dcr$Number_of_Trees <=34),]

boxplot(deercamdata_coords$Capture_Rate)
boxplot.stats(deercamdata_coords$Capture_Rate)$out
deercam_caprate<-deercam_caprate[which(deercam_caprate$Deer_Capture_Rate <=166),]
#Plot Regression Analysis of # of trees to deer capture rate
par(mar=c(5,5,4,2))
with(Trees_pDCR, plot(Capture_Rate ~ Number_of_Trees, main = "Regression Analysis of Deer Capture Rate and Number of Trees per Deployment",
                      cex.main = 2.2,
                      xlab = expression(bold("Number of Trees Per Deployment")),
                      ylab = expression(bold("Capture Rate")),
                      cex.axis = 1.3,
                      cex.lab = 1.6))
lm.outdct = lm(Capture_Rate ~ Number_of_Trees, data = Trees_Per_Dcr)
abline(lm.outdct, col="blue")
summary(lm.outdct)


#########################################################
#Create 4 point polygon to represent camera view
###########################################################
#Create data frame of the 4 points per camera
camview <- camdata_coordinates[, c(2,3,5)]
camview$X1<-(camview$NAD83_X + 6.84)
camview$Y1<-(camview$NAD83_Y + 18.79)
camview$X2<-(camview$NAD83_X)
camview$Y2<-(camview$NAD83_Y + 20)
camview$X3<-(camview$NAD83_X - 6.84)
camview$Y3<-(camview$NAD83_Y + 18.79)

camview1<- camdata_coordinates [,c(2,3,5)]
camview1[28:54,]<-(camview[1:27, c(4:5,3)])
camview1[55:81,]<-(camview[1:27, c(6:7,3)])
camview1[82:108,]<-(camview[1:27, c(8:9,3)])

camview_list<-split(camview1, camview1$Deployment)
camview_list<-lapply(camview_list, function(x) {x["Deployment"]<- NULL; x})

#create sp object and convert coords to polygon to prepare for 
cvpp <- lapply(camview_list, Polygon)

#add id variable
cvp<-lapply(seq_along(cvpp), function(i) Polygons(list(cvpp[[i]]),ID = names(camview_list)[i]))

#create sp object
camview_spo<-SpatialPolygons(cvp, proj4string = CRS(proj4string(SIGEOtrees)))

#Create spdf with IDs (one unique ID per poly) and plot polygons
camview_spo.df<-SpatialPolygonsDataFrame(camview_spo,data.frame(id = unique(camview1$Deployment),row.names = unique(camview1$Deployment)))
plot(camview_spo.df, add = T)

#Cut out tree data from within polygons
clip_polys<-intersect(SIGEOsmall,camview_spo.df)
plot(clip_polys)
cvtrees<-as.data.frame(clip_polys)

#Pull and total the # of trees per deployment and change column names
cvtreecount<-cvtrees[,c(4,28)]
cvtreecount1<-aggregate(cvtreecount[,1], by = list(cvtreecount$d),sum)
colnames(cvtreecount1)[2]<-"Number_of_Trees"
colnames(cvtreecount1)[1]<-"Deployment"

###################################################################################
#Analyse relationship between # of Trees in cameras view with Species Capture Rate
###################################################################################
#Tree count vs Deer cr
#Merge tree count per deployment with deer capture rate per deployment by Deployment
cvTrees_per_Dcr<-merge(deercamdata_coords, cvtreecount1, by = "Deployment")

#Boxplot and remove outliers
boxplot.stats(cvTrees_per_Dcr$Number_of_Trees)$out
boxplot.stats(cvTrees_per_Dcr$Capture_Rate)$out

#Remove outliers
cvTrees_T<-cvTrees_per_Dcr[which(cvTrees_per_Dcr$Number_of_Trees <=38),]
cvTrees_Dcr<-cvTrees_per_Dcr[which(cvTrees_per_Dcr$Capture_Rate <=166),]
cvTrees_pDCR<-merge(cvTrees_Dcr, cvTrees_T, by = "Deployment")

#Plot Regression Analysis of # of trees to deer capture rate
par(mar=c(5,5,4,2))
with(cvTrees_pDCR, plot(Capture_Rate.x ~ Number_of_Trees.x, main = "Regression Analysis of Deer Capture Rate and Number of Trees per Deployment",
                        cex.main = 2.1,
                        xlab = expression(bold("Trees Per Deployment")),
                        ylab = expression(bold("Capture Rate")),
                        cex.axis = 1.3,
                        cex.lab = 1.6))
lm.outdctc = lm(Capture_Rate ~ Number_of_Trees, data = cvTrees_pDCR)
abline(lm.outdctc, col="blue")
summary(lm.outdctc)

#Post Rsqrd value on plot
Rsquared<-summary(lm.outdctc)$r.squared
text(13.30551,126.1788, as.expression(substitute(italic(R)^2 == r,list(r=round(Rsquared,3)))), cex = 1.5)

#######################################
#Tree count vs Bear cr
#plot boxplot to identify outliers
boxplot(cvtreecount1$Number_of_Trees)
boxplot.stats(cvtreecount1$Number_of_Trees)$out

#Remove outliers
cvTrees_T<-cvtreecount1[which(cvtreecount1$Number_of_Trees <=38),]

cvTrees_pBCR<-merge(cvTrees_T, bearcamdata_coords, by = "Deployment")

#Plot # of trees to bear cap rate
par(mar=c(5,5,4,2))
with(cvTrees_pBCR, plot(Capture_Rate ~ Number_of_Trees, main = "Regression Analysis of Bear Capture Rate and Number of Trees per Deployment",
                        cex.main = 2.1,
                        xlab = expression(bold("Trees Per Deployment")),
                        ylab = expression(bold("Capture Rate")),
                        cex.axis = 1.3,
                        cex.lab = 1.6))
lm.outbctc = lm(Capture_Rate ~ Number_of_Trees, data = cvTrees_pBCR)
abline(lm.outbctc, col="blue")
summary(lm.outbctc)

#Post Rsqrd value on plot
Rsquared<-summary(lm.outbctc)$r.squared
text(12.37023,16.05302, as.expression(substitute(italic(R)^2 == r,list(r=round(Rsquared,3)))), cex = 1.5)

########################################
#Tree count vs Squirrel CR
#Merge # of trees and Squirrel CR
cvTrees_per_Scr<-merge(grsqrlcamdata_coords, cvtreecount1, by = "Deployment")

#Boxplot and remove outliers
boxplot(cvTrees_per_Scr$Number_of_Trees)
boxplot.stats(cvTrees_per_Scr$Number_of_Trees)$out

#Remove outliers
cvTrees_pSCR<-cvTrees_per_Scr[which(cvTrees_per_Scr$Number_of_Trees <=38),]


#Plot # of trees to bear cap rate
par(mar=c(5,5,4,2))
with(cvTrees_pSCR, plot(Capture_Rate ~ Number_of_Trees, main = "Regression Analysis of Squirrel Capture Rate and Number of Trees per Deployment",
                        cex.main = 2,
                        xlab = expression(bold("Trees Per Deployment")),
                        ylab = expression(bold("Capture Rate")),
                        cex.axis = 1.3,
                        cex.lab = 1.6))
lm.outsctc = lm(Capture_Rate ~ Number_of_Trees, data = cvTrees_pSCR)
abline(lm.outsctc, col="blue")
summary(lm.outsctc)

#Post Rsqrd value on plot
Rsquared<-summary(lm.outsctc)$r.squared
text(12.37023,97.61948, as.expression(substitute(italic(R)^2 == r,list(r=round(Rsquared,3)))), cex = 1.5)

########################################
#Tree count vs Raccoon cr
#Merge # of trees and Raccoon CR
cvTrees_per_Rcr<-merge(raccooncamdata_coords, cvtreecount1, by = "Deployment")

#Boxplot and remove outliers
boxplot(cvTrees_per_Rcr$Number_of_Trees)
boxplot.stats(cvTrees_per_Rcr$Number_of_Trees)$out
cvTrees_pRCR<-cvTrees_per_Rcr[which(cvTrees_per_Rcr$Number_of_Trees <=38),]

#Plot # of trees to bear cap rate
par(mar=c(5,5,4,2))
with(cvTrees_pRCR, plot(Capture_Rate ~ Number_of_Trees, main = "Regression Analysis of Raccoon Capture Rate and Number of Trees per Deployment",
                        cex.main = 2,
                        xlab = expression(bold("Trees Per Deployment")),
                        ylab = expression(bold("Capture Rate")),
                        cex.axis = 1.3,
                        cex.lab = 1.6))
lm.outrctc = lm(Capture_Rate ~ Number_of_Trees, data = cvTrees_pRCR)
abline(lm.outrctc, col="blue")
summary(lm.outrctc)

#Post Rsqrd value on plot
Rsquared<-summary(lm.outrctc)$r.squared
text(12.37023,11.21435, as.expression(substitute(italic(R)^2 == r,list(r=round(Rsquared,3)))), cex = 1.5)

##########################################
#Deer EDD vs Tree Count
#Merge Deer EDD and tree count
DEDD_TC<-merge(Deer_EDD, cvtreecount1, by = "Deployment")

boxplot(DEDD_TC$Number_of_Trees)
boxplot.stats(DEDD_TC$Number_of_Trees)$out
DEDD_TC1<-DEDD_TC[which(DEDD_TC$Number_of_Trees <=38),]

#Plot # of trees to Deer Estimated Detection Distance
par(mar=c(5,5,4,2))
with(DEDD_TC1, plot(ESW_EDR ~ Number_of_Trees, main = "Analysis of Trees per Deployment on Deer Detection Distance",
                    cex.main = 2,
                    xlab = expression(bold("Trees Per Deployment")),
                    ylab = expression(bold("Estimated Detection Distance")),
                    cex.axis = 1.3,
                    cex.lab = 1.6))
lm.outdeddt = lm(ESW_EDR ~ Number_of_Trees, data = DEDD_TC1)
abline(lm.outdeddt, col="blue")
summary(lm.outdeddt)

#Post Rsqrd value on plot
Rsquared<-summary(lm.outdeddt)$r.squared
text(12.37023,12.21435, as.expression(substitute(italic(R)^2 == r,list(r=round(Rsquared,3)))), cex = 1.5)

#######################################
#Oak Tree Data 
#######################################
#Pull Oak Tree Data from grid
Oak_Trees<-subset(SIGEOsmall,sp %in% c('qual','quru','quco','qufa','qupr','quve','qusp','qumi'))
plot(Oak_Trees, pch = 19)
#plot camera locations in red
plot(trapxy_sp, pch = 22, col = "red", add = T)

#add column to study site tree info that divides trees into 5 color size cateories
Oak_Trees$Size_Category[Oak_Trees$dbh <150]<-'461' #turqoise
Oak_Trees$Size_Category[Oak_Trees$dbh >150]<-'68' #dark blue
Oak_Trees$Size_Category[Oak_Trees$dbh >300]<-'47' #yellow
Oak_Trees$Size_Category[Oak_Trees$dbh >600]<-'139' #green
Oak_Trees$Size_Category[Oak_Trees$dbh >900]<-'8' #gray
Oak_Trees$Size_Category[Oak_Trees$dbh >1200]<-'550' #pink

#plot Oak Tree sizes by color
par(mar=c(5,17,4,2))
plot(Oak_Trees ,pch = 19, col = Oak_Trees$Size_Category, add = T)

#Legend matching color to size 
legend(747285,4309044, legend = c("< 15cm","> 15cm","> 30cm","> 60cm","> 90cm","> 120cm"), col = c("461", "68", "47","139", "8", "550"), pch = 19, title = "DBH of Oak Trees", bty = 'n')


########################################################################
#Regression Analysis of Oak Trees per Deployment and Deer Capture Rate
########################################################################
#Cut out oak tree data from within the cones
polyoaktrees<- intersect(Oak_Trees, clip_polys)
plot(polyoaktrees)
polyoaktreesdf<-as.data.frame(polyoaktrees)

#Pull # of oaks out of each deployment and rename columns to prepare for merge
oakcount<-polyoaktreesdf[,c(4,33)]
oakcount1<-aggregate(oakcount[,1],by=list(oakcount$Deployment), sum)
colnames(oakcount1)[2]<-"Num_Oak_Trees"
colnames(oakcount1)[1]<-"Deployment"


#Merge number of oak trees within buffers with deer capture rate
Oaks_Per_Dcr<-merge(deercamdata_coords, oakcount1,  by = "Deployment", all.x = TRUE)
Oaks_Per_Dcr$Num_Oak_Trees[is.na (Oaks_Per_Dcr$Num_Oak_Trees)] = 0

#Boxplot of Number of oak trees and remove outliers
boxplot(Oaks_Per_Dcr$Capture_Rate)
boxplot.stats(Oaks_Per_Dcr$Capture_Rate)$out
Oaks_Per_Dcr1<-Oaks_Per_Dcr[which(Oaks_Per_Dcr$Capture_Rate <166),]

#Plot Regression Analysis of # of oak trees to deer capture rate
par(mar=c(5,5,4,2))
with(Oaks_Per_Dcr1, plot(Capture_Rate ~ Num_Oak_Trees, main = "Regression Analysis of Deer Capture Rate and Number of Oak Trees per Deployment",
                         cex.main = 2.2,
                         xlab = expression(bold("Number of Oaks Per Deployment")),
                         ylab = expression(bold("Capture Rate")),
                         cex.axis = 1.3,
                         cex.lab = 1.6))
lm.outodc = lm(Capture_Rate ~ Num_Oak_Trees, data = Oaks_Per_Dcr)
abline(lm.outodc, col="blue")
summary(lm.outodc)

##############################################################################
#Regression Analysis of Oak Trees Per Deployment and Gray Squirrel Capture Rate
##############################################################################
#Merge number of oak trees within buffers with gray squirrel capture rate
Oaks_Per_GrSqCR<-merge(grsqrlcamdata_coords, oakcount1,  by = "Deployment", all.x = TRUE)
Oaks_Per_GrSqCR$Num_Oak_Trees[is.na (Oaks_Per_GrSqCR$Num_Oak_Trees)] = 0

#Boxplot of Number of oak trees and remove outliers
boxplot(Oaks_Per_GrSqCR$Capture_Rate)
boxplot.stats(Oaks_Per_GrSqCR$Capture_Rate)$out
Oaks_Per_GrSqCR1<-Oaks_Per_GrSqCR[which(Oaks_Per_GrSqCR$Capture_Rate <52),]

#Plot Regression Analysis of # of oak trees to gray squirrel capture rate
par(mar=c(5,5,4,2))
with(Oaks_Per_GrSqCR1, plot(Capture_Rate ~ Num_Oak_Trees, main = "Regression Analysis of Gray Squirrel Capture Rate and Number of Oak Trees per Deployment",
                            cex.main = 2.2,
                            xlab = expression(bold("Number of Oaks Per Camera")),
                            ylab = expression(bold("Capture Rate")),
                            cex.axis = 1.3,
                            cex.lab = 1.6))
lm.outogsc = lm(Capture_Rate ~ Num_Oak_Trees, data = Oaks_Per_GrSqCR)
abline(lm.outogsc, col="blue")
summary(lm.outogsc)

################################################
#Calculating Moran's I - Spatial Autocorrelation
################################################
library(ape)
#Deer MI
#Pull deer capture rate and coord data per deployment
Deer_MI<-deercamdata_coords[,c(1,7,4,3)]

#Generate distance matrix
Deer_MI1<-as.matrix(dist(cbind(Deer_MI$NAD83_X, Deer_MI$NAD83_Y)))

#take inverse of matrix values and replace the diagonal entries with zero
Deer_MI_inv<-1/Deer_MI1
diag(Deer_MI_inv)<-0

Deer_MI_inv[1:5, 1:5]

#Calculate Moran's I 
DeerMoran<-Moran.I(Deer_MI$Capture_Rate, Deer_MI_inv)


#################################
#Bear MI
#Pull bear capture rate and coord data per deployment
Bear_MI<-bearcamdata_coords[,c(1,7,3,4)]

#Generate distance matrix
Bear_MI1<-as.matrix(dist(cbind(Bear_MI$NAD83_X, Bear_MI$NAD83_Y)))

#take inverse of matrix values and replace the diagonal entries with zero
Bear_MI_inv<-1/Bear_MI1
diag(Bear_MI_inv)<-0

Bear_MI_inv[1:5, 1:5]

#Calculate Moran's I 
BearMoran<-Moran.I(Bear_MI$Capture_Rate, Bear_MI_inv)


#################################
#Squirrel MI
#Pull squirrel capture rate and coord data per deployment
Sqrl_MI<-grsqrlcamdata_coords[,c(1,7,3,4)]

#Generate distance matrix
Sqrl_MI1<-as.matrix(dist(cbind(Sqrl_MI$NAD83_X, Sqrl_MI$NAD83_Y)))

#take inverse of matrix values and replace the diagonal entries with zero
Sqrl_MI_inv<-1/Sqrl_MI1
diag(Sqrl_MI_inv)<-0

Sqrl_MI_inv[1:5, 1:5]

#Calculate Moran's I 
SqrlMoran<-Moran.I(Sqrl_MI$Capture_Rate, Sqrl_MI_inv)


#################################
#Raccoon MI
#Pull raccoon capture rate and coord data per deployment
Rac_MI<-raccooncamdata_coords[,c(1,7,3,4)]

#Generate distance matrix
Rac_MI1<-as.matrix(dist(cbind(Rac_MI$NAD83_X, Rac_MI$NAD83_Y)))

#take inverse of matrix values and replace the diagonal entries with zero
Rac_MI_inv<-1/Rac_MI1
diag(Rac_MI_inv)<-0

Rac_MI_inv[1:5, 1:5]

#Calculate Moran's I 
RacMoran<-Moran.I(Rac_MI$Capture_Rate, Rac_MI_inv)


######################################################
#Variograms for Deer, Bear, Gray Squirrel, and Raccoon
######################################################
#Deer Variogram
library(geoR)
#Use pulled deer coords to create number matrix
Deervario<-as.matrix(Deer_MI[,c(4,3)])

#Plot variogram based on distance bins
Dvario<-variog(coords = Deervario, data = Deer_MI$Capture_Rate, breaks = seq(0, 150,7), 
               option = "bin")
plot(Dvario)

#################################
#variogram of Gray Squirrel

#Use pulled coords to create matrix
Squirrelvario<-as.matrix(Sqrl_MI[,c(3,4)])

#Plot variogram based on distance bins
Sqrlvario<-variog(coords = Squirrelvario, data = Sqrl_MI$Capture_Rate, breaks = seq(0, 150,7), 
                  option = "bin")
plot(Sqrlvario)


##################################
#Black Bear Variogram
#Pulled bear coords
Bbearvario<-as.matrix(Bear_MI[,c(3,4)])

#Plot variogram based on distance bins
Bearvario<-variog(coords = Bbearvario, data = Bear_MI$Capture_Rate, breaks = seq(0, 150,7), 
                  option = "bin")
plot(Bearvario)

##################################
#Raccoon Variogram

#Pull raccoon coords
Araccoonvario<-as.matrix(Raccoon_MI[,c(3,4)])

#Plot variogram based on distance bins
Raccoonvario<-variog(coords = Araccoonvario, data = Raccoon_MI$Capture_Rate, breaks = seq(0, 150,7), 
                     option = "bin")
plot(Raccoonvario)

#########################################
#Plotting Covariance
#########################################
#Calculating Deer Covariance
Dcr = Deer_MI$Capture_Rate
Ddist<-as.data.frame(Deervario)
cov(Dcr, Ddist)