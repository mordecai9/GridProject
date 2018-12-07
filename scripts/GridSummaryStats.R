#This file calculates various summary statistics and figures for the SCBI Camera Grid project. It should call in the clean data object called "camdataAllClean" from the data folder in the git repository (which brings in a df called "camdataMSeq". This object is generated in the script 'CamDataPrep.R". Note that sequences here have already been defined with a 10min independence threshold.

library(reshape)
library(dplyr)
library(tidyverse)
load("data/CamdataAllClean.RData")

alldata <- camdataMSeq #shorter name for ease of coding

#Here I remove unwanted species and categories right away.
badSP <- c("Animal Not on List", "Camera Misfire", "Camera Trapper", "Corvus brachyrhynchos", "Corvus corax", "Cyanocitta cristata", "Homo sapiens", "Meleagris gallopavo", "No Animal", "Other Bird species", "Owl species", "Raptor Species", "Time Lapse", "Unknown Animal", "Unknown Bird", "Unknown Canid", "Unknown Skunk_Badger")
alldata <- subset(alldata, !alldata$Species.Name %in% badSP)
levels(as.factor(alldata$Species.Name))

#Summarize Sequence Counts and Capture Rates for Whole Grid, All Seasons --------

#create summary table of frequencies of sequences and convert it to a data frame
camdata_summary<-as.data.frame(table(alldata$Deployment_Name, 
                             alldata$Species.Name))

#change column names
colnames(camdata_summary)<-c("Deployment_Name", "Species", "Number_of_Sequences")

#Pull out season IDs for each Deployment Name
seasonID <- alldata %>%
  select(Deployment_Name, Season) %>%
  distinct

#get unique camera effort values from sequence data for each deployment. This immediately excludes cameras that had no sequences recorded at all, which is OK.
effort <- alldata %>%
  select(Deployment_Name, Deploy.Duration) %>%
  distinct
effort <- merge (effort, seasonID, by = "Deployment_Name")
totEffort <- sum(as.numeric(effort$Deploy.Duration)) #Total camera nights over all 4 seasons
seasonEffort <- as.data.frame(tapply(effort$Deploy.Duration, effort$Season, sum)) #Camera effort within each season
names(seasonEffort)[1] <- "camNights"

#Overall Capture Rate per species, across everything (NOT THE MEAN)
CRAll <- as.data.frame(tapply(camdata_summary$Number_of_Sequences, camdata_summary$Species, sum))
names(CRAll)[1] <- "TotalSeqs"
CRAll$TotCR <- (CRAll$TotalSeqs/totEffort)*100

#Overall CapRates by Season (not MEANS)
camdata_summary <- merge(camdata_summary, seasonID, by = "Deployment_Name" )

CRSeasonTots <- as.data.frame(tapply(camdata_summary$Number_of_Sequences, list(camdata_summary$Species, camdata_summary$Season), sum))
names(CRSeasonTots) <- c("FallSeqs", "SpringSeqs", "SummerSeqs", "WinterSeqs")
CRSeasonTots$FallCR <- CRSeasonTots$FallSeqs/seasonEffort$camNights[1]*100
CRSeasonTots$SpringCR <- CRSeasonTots$SpringSeqs/seasonEffort$camNights[2]*100
CRSeasonTots$SummerCR <- CRSeasonTots$SummerSeqs/seasonEffort$camNights[3]*100
CRSeasonTots$WinterCR <- CRSeasonTots$WinterSeqs/seasonEffort$camNights[4]*100

CRSummary <- merge(CRSeasonTots, CRAll, by = "row.names") #Final Season by Species summary of capture rates including overall pooled sequences capture rate over whole study

#Merge in the camera effort and calculate capture rate for each species and each deployment, as rate of sequences per 100 camera nights
camdata_summary <- merge(camdata_summary, effort, by = "Deployment_Name" )
camdata_summary$Deploy.Duration <- as.numeric(camdata_summary$Deploy.Duration)
camdata_summary$CR <- (camdata_summary$Number_of_Sequences/camdata_summary$Deploy.Duration) *100
camdata_summary$Season.y <- NULL
names(camdata_summary)[4] <- "Season"

#Shift table so species are along top and values are # of sequences, fill in blanks with zeros.Not sure if I need this format, but it is here just in case.
camdata_summary_seqs<- spread(
  data = camdata_summary,
  key = Species,
  value = Number_of_Sequences,
  fill = 0
)
camdata_summary_seqs <- dplyr::select(camdata_summary_seqs, -("CR")) #unwanted column

#Shift table so species are along top and values are capture rates, fill in blanks with zeros. Not sure if I need this format, but it is here just in case.
camdata_summary_caprate<- spread(
  data = camdata_summary,
  key = Species,
  value = CR,
  fill = 0
)
camdata_summary_caprate <- dplyr::select(camdata_summary_caprate, -("Number_of_Sequences")) #unwanted column


#____________________________________________________
# Summary Stats - CapRates Overall and By Season ------------------------
#____________________________________________________

#calculating average capture rate for each species across the grid. This can be seen as a baseline number against which we can compare numbers at individual cameras. Note that in this case, this is across all 4 seasons, so we are averaging over more than 100 data points. Can probably do this using summarize like I did for seasonal calcs

baselinesAll<- as.data.frame(tapply(camdata_summary$CR, camdata_summary$Species, mean))
names(baselinesAll)[1] <- "Mean_CR"

#SD of average number of sequences across the grid for each species
baselinesAll_sd <- tapply(camdata_summary$CR, camdata_summary$Species, sd)
baselinesAll$sd<- baselinesAll_sd
names(baselinesAll)[2] <- "StDev_CR"

#Capture rate by species by camera, here with species as rows
CR_bD_bS <- as.data.frame(tapply(camdata_summary$CR, list(camdata_summary$Species, camdata_summary$Deployment_Name), mean))

#Min and max CR for each species across the grid
baselinesAll$minCR <- apply(CR_bD_bS,1,min)
names(baselinesAll)[3]<- "Min_CR"
baselinesAll$maxCR <- apply(CR_bD_bS,1,max)
names(baselinesAll)[4]<-"Max_CR"

#Add columns to CR_bD_bS for standard deviation, mean, max and min.
CR_bD_bS$Mean_CR <- baselinesAll$Mean_CR
CR_bD_bS$StDev_CR<-baselinesAll$StDev_CR
CR_bD_bS$Min_CR<-baselinesAll$Min_CR
CR_bD_bS$Max_CR<-baselinesAll$Max_CR

baselinesAll <- merge(baselinesAll, CRSummary, by.x = "row.names", by.y = "Row.names")

#Capture Rates and Proportion of Detections, by Season and Species####

#Note that the total count for cameras in Fall was 25, winter was 25, summer and spring were both 27, so need to take this into account when lookign at proportion of cameras with detections. 

SeasonsCR <- camdata_summary %>%
  group_by(.dots = c("Season", "Species")) %>%
  summarize(meanCR = mean(CR),
            sdCR = sd(CR),
            minCR = min(CR),
            maxCR = max(CR),
            seCR = sd(CR)/sqrt(27),
            lowCI = mean(CR)-(1.96*(sd(CR)/sqrt(27))),
            highCI = mean(CR)+(1.96*(sd(CR)/sqrt(27))),            
            DetCount = length(which(Number_of_Sequences >0)))
SeasonsCR$lowCI <- ifelse(SeasonsCR$lowCI >= 0, SeasonsCR$lowCI, 0)

#_______________________________
#Calculate Proportion of Camera Deployments Confirming Each Species####
#_______________________________

#Create replicate data frame and remove mean,sd,min,max columns
CR_bD_bS2<- dplyr::select(CR_bD_bS,-c("Mean_CR", "StDev_CR", "Min_CR", "Max_CR"))

#convert all values >0 to 1
CR_bD_bS2<-as.data.frame((ifelse(CR_bD_bS2==0,0,1)))

#create column of camera capture proportion per species for whole grid to original data frame. SHOULD THINK ABOUT ADDING THIS INFO TO CR SUMMARY NOT CR_bD_bS
CR_bD_bS$Prop_Grid<-rowSums(CR_bD_bS2)/length(unique(alldata$Deployment_Name))
CR_bD_bS$DetCount <- rowSums(CR_bD_bS2)
CR_bD_bS$TotDepl <- length(unique(camdata_summary$Deployment_Name))

#At this point I just have baseline summary values for capture rates by species, across all seasons, in a few difffert formats for various uses later.

save(camdata_summary, file = "data/camdata_summary")

#_________________________________________________________
#Boxplot of CapRates and Proportions of Deployments----
#_________________________________________________________

#Transpose headers
CRt<-t(CR_bD_bS2)

#reshape data so that each row is per deployment, per species, per capture rate
#FIGURE OUT IF I ALREADY HAVE ALL THIS. MIGHT BE REDUNDANT
CR1<-melt(CRt)
colnames(CR1)<-c ("Deployment", "Species", "Detected")

#Calculates proportion of cameras that captured each species
CR2<-subset(CR1, Detected !=0)
NumPres<-as.data.frame(table(CR2$Species))
colnames(NumPres)<-c("Species","NumPres")

#Same but for each season
labelsFall <- SeasonsCR %>%
  select(Species, Season, DetCount) %>%
  filter(Season == "Fall 2017") %>%
  mutate(FallDeployCt = 25) %>%
  mutate(FallLabel = paste(Species,"|", DetCount, "/", FallDeployCt) ) %>%
  rename(FallDetCount = DetCount)  

labelsSum <- SeasonsCR %>%
  select(Species, Season, DetCount) %>%
  filter(Season == "Summer 2017") %>%
  mutate(SumDeployCt = 27) %>%
  mutate(SumLabel = paste(Species,"|", DetCount, "/", SumDeployCt) ) %>%
  rename(SumDetCount = DetCount) 

labelsWinter <- SeasonsCR %>%
  select(Species, Season, DetCount) %>%
  filter(Season == "Winter 2017") %>%
  mutate(WinterDeployCt = 25) %>%
  mutate(WinterLabel = paste(Species,"|", DetCount, "/", WinterDeployCt) ) %>%
  rename(WinDetCount = DetCount) 

labelsSpring <- SeasonsCR %>%
  select(Species, Season, DetCount) %>%
  filter(Season == "Spring 2018") %>%
  mutate(SpringDeployCt = 27) %>%
  mutate(SpringLabel = paste(Species,"|", DetCount, "/", SpringDeployCt) ) %>%
  rename(SprDetCount = DetCount) 

#Merge proportion dataframe with rest of data
CR3<-merge(CR1, NumPres, by = "Species")

CR3["tot"]<-length(unique(CR1$Deployment))

#Adds column label which outputs the desired text labels for the boxplot. 
#NEED TO MAKE LABELS FOR OTHER SEASONS HERE. I THINK I HAVE ALL IN NEED in CR3 to DO IT DIRECTLY
CR3["labelAll"]<-paste(CR3$Species,"|", CR3$NumPres, "/", CR3$tot)
Labels_all <- CR3 %>%
  select(Species, labelAll) %>%
  distinct

camdata_summaryL <- camdata_summary %>%
  left_join(Labels_all, by = "Species") %>%
  left_join(labelsSpring,by = "Species") %>%
  left_join(labelsFall,by = "Species") %>%
  left_join(labelsSum,by = "Species") %>%
  left_join(labelsWinter,by = "Species") %>%
  select(-c(Season.y, Season.x.x, Season.y.y, Season))          
  
  


#Boxplot showing median capture rate per species, all seasons. Remember I can use "subset" here to just do a boxplot for summer, but it much more complicated to get the proportion of stations with presence in this way. Might need to redo all code above by season?? SHould consider doing this with mean in a bar lot with whiskers, instead of median.
par(mar=c(9,17,4,2))
plot<-boxplot(camdata_summaryL$CR~camdata_summaryL$labelAll,
              cex.main = 2.5,
              main = "",
              cex.lab = 1.0,
              cex.axis = 1.0,
              horizontal = T, las = 2, cex.axis = 1,
              names.arg = camdata_summaryL$Species, at=rank(tapply(camdata_summaryL$CR,camdata_summaryL$Species,median), ties.method = "random"))
mtext(expression(bold("Species")), side = 2, line = 15, cex = 1.7)
mtext(expression(bold("Total Capture Rate (events per 100 camera-nights)")), side = 1 , line = 4, cex = 1.3)

#Plots for the 4 seasons capture rate by species#### 
#Good to go, just need to figure out how to order the species either the same every time, or in order by season.
par(mfrow = c(2,2))
par(mar=c(9,17,4,2))
plot<-boxplot(camdata_summaryL$CR~camdata_summaryL$SumLabel,
              cex.main = 2.5,
              main = "",
              cex.lab = 1.0,
              cex.axis = 1.0,
              horizontal = T, las = 2, cex.axis = 1,
              ylim = c(0,450),
              subset = camdata_summaryL$Season.x == "Summer 2017",
              names.arg = camdata_summaryL$Species, at=rank(tapply(camdata_summaryL$CR,camdata_summaryL$Species,median), ties.method = "random"))
mtext(expression(bold("Species")), side = 2, line = 15, cex = 1.7)
mtext(expression(bold("Summer CR")), side = 1 , line = 4, cex = 1.3)

plot<-boxplot(data = camdata_summaryL, CR~FallLabel,
              cex.main = 2.5,
              main = "",
              cex.lab = 1.0,
              cex.axis = 1.0,
              horizontal = T, las = 2, cex.axis = 1,
              ylim = c(0,450),
              subset = Season.x == "Fall 2017",
              names.arg = camdata_summaryL$Species, at=rank(tapply(camdata_summaryL$CR,camdata_summaryL$Species,median), ties.method = "random"))
mtext(expression(bold("Fall CR")), side = 1 , line = 4, cex = 1.3)

plot<-boxplot(data = camdata_summaryL, CR~WinterLabel,
              cex.main = 2.5,
              main = "",
              cex.lab = 1.0,
              cex.axis = 1.0,
              horizontal = T, las = 2, cex.axis = 1,
              ylim = c(0,450),
              subset = Season.x == "Winter 2017",
              names.arg = camdata_summaryL$Species, at=rank(tapply(camdata_summaryL$CR,camdata_summaryL$Species,median), ties.method = "random"))
mtext(expression(bold("Species")), side = 2, line = 15, cex = 1.7)
mtext(expression(bold("Winter CR")), side = 1 , line = 4, cex = 1.3)

plot<-boxplot(data = camdata_summaryL, CR~SpringLabel,
              cex.main = 2.5,
              main = "",
              cex.lab = 1.0,
              cex.axis = 1.0,
              horizontal = T, las = 2, cex.axis = 1,
              subset = Season.x == "Spring 2018",
              ylim = c(0,450),
              names.arg = camdata_summaryL$Species, at=rank(tapply(camdata_summaryL$CR,camdata_summaryL$Species,median), ties.method = "random"))
mtext(expression(bold("Spring CR")), side = 1 , line = 4, cex = 1.3)

#______________________________________________
#Plot CapRates by Point Size on Grid####
#Remember to somehow label cameras that were not operating. This goes for Winter and Fall figures.
#_____________________________________________

#Bring in XY data for grid
gridXY <- read.csv("data/Grid_Coordinates.csv")
camdata_summary <- merge(camdata_summary, gridXY, by.x = "Deployment_Name2", by.y = "Deployment")

#_____________________________________________
##CapRates Point Size Figures- DEER----
#_____________________________________________
#Would be nice to add proportion of cameras with hits, or to label cameras that had no detections. Would also be good to have a small key to circle size, maybe showing lower quartile, median and upper quartile values? Final figures will probably just have Season names, with the species in the caption?

par(mar = c(1,1,1,1))
par(mfrow = c(2,2))

#DEER SUMMER
s=12
camdata_summary %>%
  filter(Species == "Odocoileus virginianus") %>%
  filter(Season == "Summer 2017") %>%
  with(plot(NAD83_X, NAD83_Y,main = "",
            axes = F,
            cex.main = 2.2,
            xlab = "", ylab = "",
            cex = CR/s, 
            pch = 20,  
            xlim = c(747430, 747550),
            ylim = c(4308910,4309030)))
axis(side = 2,tck = 0.02, labels = F)
axis(side = 1, tck = 0.02, labels = F)
box()
points(camdata_summary$NAD83_X, camdata_summary$NAD83_Y,
       xlab = "", ylab = "",
       cex = .8, 
       pch = 3,  
       xlim = c(747430, 747550),
       ylim = c(4308910,4309030))
text(labels = "Deer - Summer", y = 4309025, x = 747430, cex = 1.6, pos = 4)

#DEER Fall
s=12
camdata_summary %>%
  filter(Species == "Odocoileus virginianus") %>%
  filter(Season == "Fall 2017") %>%
  with(plot(NAD83_X, NAD83_Y,main = "",
            axes = F,
            cex.main = 2.2,
            xlab = "", ylab = "",
            cex = CR/s, 
            pch = 20,  
            xlim = c(747430, 747550),
            ylim = c(4308910,4309030)))
axis(side = 2,tck = 0.02, labels = F)
axis(side = 1, tck = 0.02, labels = F)
box()
points(camdata_summary$NAD83_X, camdata_summary$NAD83_Y,
       xlab = "", ylab = "",
       cex = .8, 
       pch = 3,  
       xlim = c(747430, 747550),
       ylim = c(4308910,4309030))
text(labels = "Deer - Fall", y = 4309025, x = 747430, cex = 1.6, pos = 4)

#DEER WINTER
s=12
camdata_summary %>%
  filter(Species == "Odocoileus virginianus") %>%
  filter(Season == "Winter 2017") %>%
  with(plot(NAD83_X, NAD83_Y,main = "",
            axes = F,
            cex.main = 2.2,
            xlab = "", ylab = "",
            cex = CR/s, 
            pch = 20,  
            xlim = c(747430, 747550),
            ylim = c(4308910,4309030)))
axis(side = 2,tck = 0.02, labels = F)
axis(side = 1, tck = 0.02, labels = F)
box()
points(camdata_summary$NAD83_X, camdata_summary$NAD83_Y,
       xlab = "", ylab = "",
       cex = .8, 
       pch = 3,  
       xlim = c(747430, 747550),
       ylim = c(4308910,4309030))
text(labels = "Deer - Winter", y = 4309025, x = 747430, cex = 1.6, pos = 4)

#DEER Spring
s=12
camdata_summary %>%
  filter(Species == "Odocoileus virginianus") %>%
  filter(Season == "Spring 2018") %>%
  with(plot(NAD83_X, NAD83_Y,main = "",
            axes = F,
            cex.main = 2.2,
            xlab = "", ylab = "",
            cex = CR/s, 
            pch = 20,  
            xlim = c(747430, 747550),
            ylim = c(4308910,4309030)))
axis(side = 2,tck = 0.02, labels = F)
axis(side = 1, tck = 0.02, labels = F)
box()
points(camdata_summary$NAD83_X, camdata_summary$NAD83_Y,
       xlab = "", ylab = "",
       cex = .8, 
       pch = 3,  
       xlim = c(747430, 747550),
       ylim = c(4308910,4309030))
text(labels = "Deer - Spring", y = 4309025, x = 747430, cex = 1.6, pos = 4)



#_____________________________________________
##CapRates Point Size Figures- Fox Squirrel----
#_____________________________________________
#Might want to map unknown squirrels to see if "expert review" varied in calling squirrel ID over seasons

#Fox Squirrel SUMMER
s=6
camdata_summary %>%
  filter(Species == "Sciurus niger") %>%
  filter(Season == "Summer 2017") %>%
  with(plot(NAD83_X, NAD83_Y,main = "",
            axes = F,
            cex.main = 2.2,
            xlab = "", ylab = "",
            cex = CR/s, 
            pch = 20,  
            xlim = c(747430, 747550),
            ylim = c(4308910,4309030)))
axis(side = 2,tck = 0.02, labels = F)
axis(side = 1, tck = 0.02, labels = F)
box()
points(camdata_summary$NAD83_X, camdata_summary$NAD83_Y,
       xlab = "", ylab = "",
       cex = .8, 
       pch = 3,  
       xlim = c(747430, 747550),
       ylim = c(4308910,4309030))
text(labels = "Fox Sq. - Summer", y = 4309025, x = 747430, cex = 1.6, pos = 4)

#Fox Squirrel Fall
s=6
camdata_summary %>%
  filter(Species == "Sciurus niger") %>%
  filter(Season == "Fall 2017") %>%
  with(plot(NAD83_X, NAD83_Y,main = "",
            axes = F,
            cex.main = 2.2,
            xlab = "", ylab = "",
            cex = CR/s, 
            pch = 20,  
            xlim = c(747430, 747550),
            ylim = c(4308910,4309030)))
axis(side = 2,tck = 0.02, labels = F)
axis(side = 1, tck = 0.02, labels = F)
box()
points(camdata_summary$NAD83_X, camdata_summary$NAD83_Y,
       xlab = "", ylab = "",
       cex = .8, 
       pch = 3,  
       xlim = c(747430, 747550),
       ylim = c(4308910,4309030))
text(labels = "Fox sq. - Fall", y = 4309025, x = 747430, cex = 1.6, pos = 4)

#Fox Squirrel Winter
s=6
camdata_summary %>%
  filter(Species == "Sciurus niger") %>%
  filter(Season == "Winter 2017") %>%
  with(plot(NAD83_X, NAD83_Y,main = "",
            axes = F,
            cex.main = 2.2,
            xlab = "", ylab = "",
            cex = CR/s, 
            pch = 20,  
            xlim = c(747430, 747550),
            ylim = c(4308910,4309030)))
axis(side = 2,tck = 0.02, labels = F)
axis(side = 1, tck = 0.02, labels = F)
box()
points(camdata_summary$NAD83_X, camdata_summary$NAD83_Y,
       xlab = "", ylab = "",
       cex = .8, 
       pch = 3,  
       xlim = c(747430, 747550),
       ylim = c(4308910,4309030))
text(labels = "Fox Sq. - Winter", y = 4309025, x = 747430, cex = 1.6, pos = 4)

#Fox Squirrel - Spring
s=6
camdata_summary %>%
  filter(Species == "Sciurus niger") %>%
  filter(Season == "Spring 2018") %>%
  with(plot(NAD83_X, NAD83_Y,main = "",
            axes = F,
            cex.main = 2.2,
            xlab = "", ylab = "",
            cex = CR/s, 
            pch = 20,  
            xlim = c(747430, 747550),
            ylim = c(4308910,4309030)))
axis(side = 2,tck = 0.02, labels = F)
axis(side = 1, tck = 0.02, labels = F)
box()
points(camdata_summary$NAD83_X, camdata_summary$NAD83_Y,
       xlab = "", ylab = "",
       cex = .8, 
       pch = 3,  
       xlim = c(747430, 747550),
       ylim = c(4308910,4309030))
text(labels = "Fox sq. - Spring", y = 4309025, x = 747430, cex = 1.6, pos = 4)

#_____________________________________________
##CapRates Point Size Figures- Gray Squirrel----
#_____________________________________________

par(mfrow = c(2,2))

#Gray Squirrel Summer
s=6
camdata_summary %>%
  filter(Species == "Sciurus carolinensis") %>%
  filter(Season == "Summer 2017") %>%
  with(plot(NAD83_X, NAD83_Y,main = "",
            axes = F,
            cex.main = 2.2,
            xlab = "", ylab = "",
            cex = CR/s, 
            pch = 20,  
            xlim = c(747430, 747550),
            ylim = c(4308910,4309030)))
axis(side = 2,tck = 0.02, labels = F)
axis(side = 1, tck = 0.02, labels = F)
box()
points(camdata_summary$NAD83_X, camdata_summary$NAD83_Y,
       xlab = "", ylab = "",
       cex = .8, 
       pch = 3,  
       xlim = c(747430, 747550),
       ylim = c(4308910,4309030))
text(labels = "Gray Sq. - Summer", y = 4309025, x = 747430, cex = 1.6, pos = 4)

#Gray Squirrel FALL
s=6
camdata_summary %>%
  filter(Species == "Sciurus carolinensis") %>%
  filter(Season == "Fall 2017") %>%
  with(plot(NAD83_X, NAD83_Y,main = "",
            axes = F,
            cex.main = 2.2,
            xlab = "", ylab = "",
            cex = CR/s, 
            pch = 20,  
            xlim = c(747430, 747550),
            ylim = c(4308910,4309030)))
axis(side = 2,tck = 0.02, labels = F)
axis(side = 1, tck = 0.02, labels = F)
box()
points(camdata_summary$NAD83_X, camdata_summary$NAD83_Y,
       xlab = "", ylab = "",
       cex = .8, 
       pch = 3,  
       xlim = c(747430, 747550),
       ylim = c(4308910,4309030))
text(labels = "Gray sq. - Fall", y = 4309025, x = 747430, cex = 1.6, pos = 4)

#Gray Squirrel Winter
s=6
camdata_summary %>%
  filter(Species == "Sciurus carolinensis") %>%
  filter(Season == "Winter 2017") %>%
  with(plot(NAD83_X, NAD83_Y,main = "",
            axes = F,
            cex.main = 2.2,
            xlab = "", ylab = "",
            cex = CR/s, 
            pch = 20,  
            xlim = c(747430, 747550),
            ylim = c(4308910,4309030)))
axis(side = 2,tck = 0.02, labels = F)
axis(side = 1, tck = 0.02, labels = F)
box()
points(camdata_summary$NAD83_X, camdata_summary$NAD83_Y,
       xlab = "", ylab = "",
       cex = .8, 
       pch = 3,  
       xlim = c(747430, 747550),
       ylim = c(4308910,4309030))
text(labels = "Gray sq. - Winter", y = 4309025, x = 747430, cex = 1.6, pos = 4)

#Gray Squirrel Spring
s=6
camdata_summary %>%
  filter(Species == "Sciurus carolinensis") %>%
  filter(Season == "Spring 2018") %>%
  with(plot(NAD83_X, NAD83_Y,main = "",
            axes = F,
            cex.main = 2.2,
            xlab = "", ylab = "",
            cex = CR/s, 
            pch = 20,  
            xlim = c(747430, 747550),
            ylim = c(4308910,4309030)))
axis(side = 2,tck = 0.02, labels = F)
axis(side = 1, tck = 0.02, labels = F)
box()
points(camdata_summary$NAD83_X, camdata_summary$NAD83_Y,
       xlab = "", ylab = "",
       cex = .8, 
       pch = 3,  
       xlim = c(747430, 747550),
       ylim = c(4308910,4309030))
text(labels = "Gray sq. - Spring", y = 4309025, x = 747430, cex = 1.6, pos = 4)


#_____________________________________________
##CapRates Point Size Figures- Raccoon----
#_____________________________________________

par(mfrow = c(2,2))

#Raccoon Summer
s=5
camdata_summary %>%
  filter(Species == "Procyon lotor") %>%
  filter(Season == "Summer 2017") %>%
  with(plot(NAD83_X, NAD83_Y,main = "",
            axes = F,
            cex.main = 2.2,
            xlab = "", ylab = "",
            cex = CR/s, 
            pch = 20,  
            xlim = c(747430, 747550),
            ylim = c(4308910,4309030)))
axis(side = 2,tck = 0.02, labels = F)
axis(side = 1, tck = 0.02, labels = F)
box()
points(camdata_summary$NAD83_X, camdata_summary$NAD83_Y,
       xlab = "", ylab = "",
       cex = .8, 
       pch = 3,  
       xlim = c(747430, 747550),
       ylim = c(4308910,4309030))
text(labels = "Raccoon - Summer", y = 4309025, x = 747430, cex = 1.6, pos = 4)

#Raccoon FALL
s=5
camdata_summary %>%
  filter(Species == "Procyon lotor") %>%
  filter(Season == "Fall 2017") %>%
  with(plot(NAD83_X, NAD83_Y,main = "",
            axes = F,
            cex.main = 2.2,
            xlab = "", ylab = "",
            cex = CR/s, 
            pch = 20,  
            xlim = c(747430, 747550),
            ylim = c(4308910,4309030)))
axis(side = 2,tck = 0.02, labels = F)
axis(side = 1, tck = 0.02, labels = F)
box()
points(camdata_summary$NAD83_X, camdata_summary$NAD83_Y,
       xlab = "", ylab = "",
       cex = .8, 
       pch = 3,  
       xlim = c(747430, 747550),
       ylim = c(4308910,4309030))
text(labels = "Raccoon - Fall", y = 4309025, x = 747430, cex = 1.6, pos = 4)

#Raccoon Winter
s=5
camdata_summary %>%
  filter(Species == "Procyon lotor") %>%
  filter(Season == "Winter 2017") %>%
  with(plot(NAD83_X, NAD83_Y,main = "",
            axes = F,
            cex.main = 2.2,
            xlab = "", ylab = "",
            cex = CR/s, 
            pch = 20,  
            xlim = c(747430, 747550),
            ylim = c(4308910,4309030)))
axis(side = 2,tck = 0.02, labels = F)
axis(side = 1, tck = 0.02, labels = F)
box()
points(camdata_summary$NAD83_X, camdata_summary$NAD83_Y,
       xlab = "", ylab = "",
       cex = .8, 
       pch = 3,  
       xlim = c(747430, 747550),
       ylim = c(4308910,4309030))
text(labels = "Raccoon - Winter", y = 4309025, x = 747430, cex = 1.6, pos = 4)

#Raccoon Spring
s=5
camdata_summary %>%
  filter(Species == "Procyon lotor") %>%
  filter(Season == "Spring 2018") %>%
  with(plot(NAD83_X, NAD83_Y,main = "",
            axes = F,
            cex.main = 2.2,
            xlab = "", ylab = "",
            cex = CR/s, 
            pch = 20,  
            xlim = c(747430, 747550),
            ylim = c(4308910,4309030)))
axis(side = 2,tck = 0.02, labels = F)
axis(side = 1, tck = 0.02, labels = F)
box()
points(camdata_summary$NAD83_X, camdata_summary$NAD83_Y,
       xlab = "", ylab = "",
       cex = .8, 
       pch = 3,  
       xlim = c(747430, 747550),
       ylim = c(4308910,4309030))
text(labels = "Raccoon - Spring", y = 4309025, x = 747430, cex = 1.6, pos = 4)

#_____________________________________________
##CapRates Point Size Figures- Bear----
#_____________________________________________

par(mfrow = c(2,2))

#Black Bear Summer
s=1
camdata_summary %>%
  filter(Species == "Ursus americanus") %>%
  filter(Season == "Summer 2017") %>%
  with(plot(NAD83_X, NAD83_Y,main = "",
            axes = F,
            cex.main = 2.2,
            xlab = "", ylab = "",
            cex = CR/s, 
            pch = 20,  
            xlim = c(747430, 747550),
            ylim = c(4308910,4309030)))
axis(side = 2,tck = 0.02, labels = F)
axis(side = 1, tck = 0.02, labels = F)
box()
points(camdata_summary$NAD83_X, camdata_summary$NAD83_Y,
       xlab = "", ylab = "",
       cex = .8, 
       pch = 3,  
       xlim = c(747430, 747550),
       ylim = c(4308910,4309030))
text(labels = "Black Bear - Summer", y = 4309025, x = 747430, cex = 1.6, pos = 4)

#Black Bear FALL
s=1
camdata_summary %>%
  filter(Species == "Ursus americanus") %>%
  filter(Season == "Fall 2017") %>%
  with(plot(NAD83_X, NAD83_Y,main = "",
            axes = F,
            cex.main = 2.2,
            xlab = "", ylab = "",
            cex = CR/s, 
            pch = 20,  
            xlim = c(747430, 747550),
            ylim = c(4308910,4309030)))
axis(side = 2,tck = 0.02, labels = F)
axis(side = 1, tck = 0.02, labels = F)
box()
points(camdata_summary$NAD83_X, camdata_summary$NAD83_Y,
       xlab = "", ylab = "",
       cex = .8, 
       pch = 3,  
       xlim = c(747430, 747550),
       ylim = c(4308910,4309030))
text(labels = "Black Bear - Fall", y = 4309025, x = 747430, cex = 1.6, pos = 4)

#Black Bear Winter
s=1
camdata_summary %>%
  filter(Species == "Ursus americanus") %>%
  filter(Season == "Winter 2017") %>%
  with(plot(NAD83_X, NAD83_Y,main = "",
            axes = F,
            cex.main = 2.2,
            xlab = "", ylab = "",
            cex = CR/s, 
            pch = 20,  
            xlim = c(747430, 747550),
            ylim = c(4308910,4309030)))
axis(side = 2,tck = 0.02, labels = F)
axis(side = 1, tck = 0.02, labels = F)
box()
points(camdata_summary$NAD83_X, camdata_summary$NAD83_Y,
       xlab = "", ylab = "",
       cex = .8, 
       pch = 3,  
       xlim = c(747430, 747550),
       ylim = c(4308910,4309030))
text(labels = "Black Bear - Winter", y = 4309025, x = 747430, cex = 1.6, pos = 4)

#Black Bear Spring
s=1
camdata_summary %>%
  filter(Species == "Ursus americanus") %>%
  filter(Season == "Spring 2018") %>%
  with(plot(NAD83_X, NAD83_Y,main = "",
            axes = F,
            cex.main = 2.2,
            xlab = "", ylab = "",
            cex = CR/s, 
            pch = 20,  
            xlim = c(747430, 747550),
            ylim = c(4308910,4309030)))
axis(side = 2,tck = 0.02, labels = F)
axis(side = 1, tck = 0.02, labels = F)
box()
points(camdata_summary$NAD83_X, camdata_summary$NAD83_Y,
       xlab = "", ylab = "",
       cex = .8, 
       pch = 3,  
       xlim = c(747430, 747550),
       ylim = c(4308910,4309030))
text(labels = "Black Bear - Spring", y = 4309025, x = 747430, cex = 1.6, pos = 4)
