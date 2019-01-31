#This file calculates various summary statistics and figures for the SCBI Camera Grid project. It should call in the clean data object called "camdataAllClean" from the data folder in the git repository (which brings in a df called "camdataMSeq". This object is generated in the script 'CamDataPrep.R". Note that sequences here have already been defined with a 10min independence threshold.

#This script calculated a range of summary statistics for the project, and creates a boxplot of capture rates overall and for each season. It finally creates point maps for each species and each season, with point size of each camera location based on capture rate.

#The script also saves the object "camdata_summary" which is a key object used in other analysis scripts for this project.

library(reshape)
library(dplyr)
library(tidyverse)
load("data/CamdataAllClean.RData")

alldata <- camdataMSeq #shorter name for ease of coding

#Here I remove unwanted species and categories right away.
badSP <- c("Animal Not on List", "Camera Misfire", "Camera Trapper", "Corvus brachyrhynchos", "Corvus corax", "Cyanocitta cristata", "Homo sapiens", "Meleagris gallopavo", "No Animal", "Other Bird species", "Owl species", "Raptor Species", "Time Lapse", "Unknown Animal", "Unknown Bird", "Unknown Canid", "Unknown Skunk_Badger")
alldata <- subset(alldata, !alldata$Species.Name %in% badSP)
levels(as.factor(alldata$Species.Name))

alldata$Species.Name<-sub(pattern = 'Unknown Small Rodent','Unknown small rodent',alldata$Species.Name)
alldata$Species.Name<-sub(pattern = 'Unknown Squirrel','Unknown squirrel',alldata$Species.Name)

#Summarize Sequence Counts and Capture Rates for Whole Grid, All Seasons --------

#create summary table of frequencies of sequences and convert it to a data frame
camdata_summary<-as.data.frame(table(alldata$Deployment_Name, 
                             alldata$Species.Name))

#change column names
colnames(camdata_summary)<-c("Deployment_Name", "Species", "Number_of_Sequences")

#Pull out season IDs for each Deployment Name
seasonID <- alldata %>%
  dplyr::select(Deployment_Name, Season) %>%
  distinct

#get unique camera effort values from sequence data for each deployment. This immediately excludes cameras that had no sequences recorded at all, which is OK.
effort <- alldata %>%
  dplyr::select(Deployment_Name, Deploy.Duration) %>%
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

CRSummary <- merge(CRSeasonTots, CRAll, by = "row.names") #Final Season by Species summary of capture rates including overall pooled sequences capture rate over whole study. Later in script I add in overall means across all seasons as well.

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
# Summary Stats - Mean CapRates Overall and By Season ------------------------
#____________________________________________________

#calculating average capture rate for each species across the grid. This can be seen as a baseline number against which we can compare numbers at individual cameras. Note that in this case, this is across all 4 seasons, so we are averaging over more than 100 data points. Can probably do this using summarize like I did for seasonal calcs. Note that these differ from the files above in that they are averages, not overall CR calculations.

baselinesAll<- as.data.frame(tapply(camdata_summary$CR, camdata_summary$Species, mean))
names(baselinesAll)[1] <- "Mean_CR"

#SD of average number of sequences across the grid for each species
baselinesAll_sd <- tapply(camdata_summary$CR, camdata_summary$Species, sd)
baselinesAll$sd<- baselinesAll_sd
names(baselinesAll)[2] <- "StDev_CR"

#Capture rate by species by camera, here with species as rows
CR_bD_bS <- as.data.frame(tapply(camdata_summary$CR, list(camdata_summary$Species, camdata_summary$Deployment_Name), mean))

#Min and max CR for each species across the grid. I think I could have done this the same as I did with mean and sd above.
baselinesAll$minCR <- apply(CR_bD_bS,1,min)
names(baselinesAll)[3]<- "Min_CR"
baselinesAll$maxCR <- apply(CR_bD_bS,1,max)
names(baselinesAll)[4]<-"Max_CR"

#Add columns to CR_bD_bS for standard deviation, mean, max and min.
#CR_bD_bS$Mean_CR <- baselinesAll$Mean_CR
#CR_bD_bS$StDev_CR<-baselinesAll$StDev_CR
#CR_bD_bS$Min_CR<-baselinesAll$Min_CR
#CR_bD_bS$Max_CR<-baselinesAll$Max_CR

baselinesAll <- merge(baselinesAll, CRSummary, by.x = "row.names", by.y = "Row.names")
baselinesAll$TotCV <- baselinesAll$Mean_CR/baselinesAll$StDev_CR

#Capture Rates and Proportion of Detections, by Season and Species####

#Note that the total count for cameras in Fall was 25, winter was 25, summer and spring were both 27, so need to take this into account when lookign at proportion of cameras with detections. 

SeasonsCR <- camdata_summary %>%
  group_by(.dots = c("Season", "Species")) %>%
  summarize(meanCR = mean(CR),
            sdCR = sd(CR),
            minCR = min(CR),
            maxCR = max(CR),
            DetCount = length(which(Number_of_Sequences >0)))
            
FWCalcs <- camdata_summary %>%
  filter(Season == "Fall 2017" | Season == "Winter 2017") %>%
  group_by(.dots = c("Season", "Species")) %>%
     summarize(seCR = sd(CR)/sqrt(25),
            lowCI = mean(CR)-(1.96*(sd(CR)/sqrt(25))),
            highCI = mean(CR)+(1.96*(sd(CR)/sqrt(25))))            
          
SSCalcs <- camdata_summary %>%
  filter(Season == "Summer 2017" | Season == "Spring 2018") %>%
  group_by(.dots = c("Season", "Species")) %>%
  summarize(seCR = sd(CR)/sqrt(27),
            lowCI = mean(CR)-(1.96*(sd(CR)/sqrt(27))),
            highCI = mean(CR)+(1.96*(sd(CR)/sqrt(27))))      

NewCalcs <- rbind(FWCalcs, SSCalcs)
SeasonsCR <- SeasonsCR %>%
  left_join(NewCalcs, by = c("Season", "Species")) 
  
SeasonsCR$lowCI <- ifelse(SeasonsCR$lowCI >= 0, SeasonsCR$lowCI, 0) #make negative CI limits 0.
SeasonsCR$CV <- SeasonsCR$sdCR/SeasonsCR$meanCR



#_______________________________
#Calculate Proportion of Camera Deployments Confirming Each Species####
#_______________________________

#Create replicate data frame and remove mean,sd,min,max columns
#CR_bD_bS2<- dplyr::select(CR_bD_bS,-c("Mean_CR", "StDev_CR", "Min_CR", "Max_CR"))

#convert all values >0 to 1
CR_bD_bS2<-as.data.frame((ifelse(CR_bD_bS==0,0,1)))
CR_bD_bS2<- CR_bD_bS2[,-c(105:107)]

#create column of camera capture proportion per species for whole grid to original data frame. SHOULD THINK ABOUT ADDING THIS INFO TO CR SUMMARY NOT CR_bD_bS
CR_bD_bS$Prop_Grid<-rowSums(CR_bD_bS2)/length(unique(alldata$Deployment_Name))
CR_bD_bS$DetCount <- rowSums(CR_bD_bS2)
CR_bD_bS$TotDepl <- length(unique(camdata_summary$Deployment_Name))

#At this point I just have baseline summary values for capture rates by species, across all seasons, in a few difffert formats for various uses later.

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

#Same but for each season. Would be better to have the deployment count generated instead of manually inserted
labelsFall <- SeasonsCR %>%
  dplyr::select(Species, Season, DetCount) %>%
  filter(Season == "Fall 2017") %>%
  mutate(FallDeployCt = 25) %>%
  mutate(FallLabel = paste(Species,"(", DetCount, "/", FallDeployCt, ")") ) %>%
  dplyr::rename(FallDetCount = DetCount)  

labelsSum <- SeasonsCR %>%
  dplyr::select(Species, Season, DetCount) %>%
  filter(Season == "Summer 2017") %>%
  mutate(SumDeployCt = 27) %>%
  mutate(SumLabel = paste(Species,"(", DetCount, "/", SumDeployCt, ")") ) %>%
  dplyr::rename(SumDetCount = DetCount) 

labelsWinter <- SeasonsCR %>%
  dplyr::select(Species, Season, DetCount) %>%
  filter(Season == "Winter 2017") %>%
  mutate(WinterDeployCt = 25) %>%
  mutate(WinterLabel = paste(Species,"(", DetCount, "/", WinterDeployCt, ")") ) %>%
  dplyr::rename(WinDetCount = DetCount) 

labelsSpring <- SeasonsCR %>%
  dplyr::select(Species, Season, DetCount) %>%
  filter(Season == "Spring 2018") %>%
  mutate(SpringDeployCt = 27) %>%
  mutate(SpringLabel = paste(Species,"(", DetCount, "/", SpringDeployCt, ")") ) %>%
  dplyr::rename(SprDetCount = DetCount) 

#Merge proportion dataframe with rest of data
CR3<-merge(CR1, NumPres, by = "Species")

CR3["tot"]<-length(unique(CR1$Deployment))

#Adds column label which outputs the desired text labels for the boxplot. 

CR3["labelAll"]<-paste(CR3$Species,"(", CR3$NumPres, "/", CR3$tot, ")")
Labels_all <- CR3 %>%
  dplyr::select(Species, labelAll) %>%
  distinct

camdata_summaryL <- camdata_summary %>%
  left_join(Labels_all, by = "Species") %>%
  left_join(labelsSpring,by = "Species") %>%
  left_join(labelsFall,by = "Species") %>%
  left_join(labelsSum,by = "Species") %>%
  left_join(labelsWinter,by = "Species") %>%
  dplyr::select(-c(Season.y, Season.x.x, Season.y.y, Season))          
  
  


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
              names.arg = camdata_summaryL$Species, at=rank(tapply(camdata_summaryL$CR,camdata_summaryL$Species, median), ties.method = "random"))
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
              names.arg = camdata_summaryL$Species, at=rank(tapply(camdata_summaryL$CR,camdata_summaryL$Species, median), ties.method = "random"))
mtext(expression(bold("Fall CR")), side = 1 , line = 4, cex = 1.3)

plot<-boxplot(data = camdata_summaryL, CR~WinterLabel,
              cex.main = 2.5,
              main = "",
              cex.lab = 1.0,
              cex.axis = 1.0,
              horizontal = T, las = 2, cex.axis = 1,
              ylim = c(0,450),
              subset = Season.x == "Winter 2017",
              names.arg = camdata_summaryL$Species, at=rank(tapply(camdata_summaryL$CR,camdata_summaryL$Species, median), ties.method = "random"))
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
              names.arg = camdata_summaryL$Species, at=rank(tapply(camdata_summaryL$CR,camdata_summaryL$Species, median), ties.method = "random"))
mtext(expression(bold("Spring CR")), side = 1 , line = 4, cex = 1.3)


#______________________________________________
#Plot CapRates for Each Season and Species in One Plot####
#_____________________________________________
#Need to figure out if I want just the one graph or if I want to show both the whole graph, and the zoomed into the rare species graphs as well.

library(ggplot2)
SeasonsCR$Species <- factor(SeasonsCR$Species, 
                           levels = baselinesAll$Row.names[order(baselinesAll$Mean_CR)])
SeasonsCR$Species

a <-
  ggplot(data = SeasonsCR, aes(y = meanCR, x = Species, color = Season)) +
  geom_point(shape = 1.4,
             size = 2,
             position = position_dodge(.2)) +
  geom_errorbar(
    aes(ymin = lowCI, ymax = highCI),
    size = .5,
    width = .1,
    position = position_dodge(.2)) +
  coord_flip()+
  theme_classic() +
  ylab("Mean capture rate") +
  theme(panel.grid.major.x = element_line(color="gray", size = .5, linetype = "dashed"),
        panel.grid.minor.x = element_line(color="gray", size = .5,linetype = "dashed"), legend.position = "none")
  

a
ggsave("results/SeasonCRAll.tiff", width = 6.5, height = 4)

rare <- c("Mustela frenata", "Glaucomys volans", "Vulpes vulpes", "Lynx rufus", "Mephitis mephitis", "Urocyon cinereoargenteus", "Canis latrans", "Didelphis virginianus")
b <-
  ggplot(data = filter(SeasonsCR, Species %in% rare), aes(y = meanCR, x = Species, color = Season)) +
  geom_point(shape = 1.4,
             size = 2,
             position = position_dodge(.2)) +
  geom_errorbar(
    aes(ymin = lowCI, ymax = highCI),
    size = .5,
    width = .1,
    position = position_dodge(.2)) +
  coord_flip()+
  theme_classic() +
  xlab("")+
  ylab("Mean capture rate") +
  theme(panel.grid.major.x = element_line(color="gray", size = .5, linetype = "dashed"),
        panel.grid.minor.x = element_line(color="gray", size = .5,linetype = "dashed"))

b
ggsave("results/SeasonCR_Rare.tiff", width = 6.5, height = 4)

multiplot(a,b, cols = 2) #http://www.cookbook-r.com/Graphs/Multiple_graphs_on_one_page_(ggplot2)/


#______________________________________________
#Plot CapRates by Point Size on Grid####
#Remember to somehow label cameras that were not operating. This goes for Winter and Fall figures.
#_____________________________________________

#Bring in XY data for grid
gridXY <- read.csv("data/Grid_Coordinates.csv")
#Create file with coding in each type for deployments
deployCodes <- alldata %>%
  dplyr::select(Deployment_Name, Deployment_Name2) %>%
  distinct
#Add basic deployment code into camdata_summary
camdata_summary <- merge(camdata_summary, deployCodes, by.x = "Deployment_Name")

camdata_summary <- merge(camdata_summary, gridXY, by.x = "Deployment_Name2", by.y = "Deployment")

save(camdata_summary, file = "data/camdata_summary")

#_____________________________________________
##CapRates Point Size Figures- DEER----
#_____________________________________________
#Added info to labels in deer figure indicating proportion of cameras with detections - need to add to other species. Would also be good to have a small key to circle size, maybe showing lower quartile, median and upper quartile values? Final figures will probably just have Season names, with the species in the caption? 

#call out failed points for Fall deployments
idxf <- which(gridXY$Deployment == "505" | gridXY$Deployment == "206")
#call out failed points for Winter deployments
idxw <- which(gridXY$Deployment == "101" | gridXY$Deployment == "401")

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
text(labels = paste("Deer - Summer", "(", labelsSum$SumDetCount[which(labelsSum$Species == "Odocoileus virginianus")],"/", labelsSum$SumDeployCt[1], ")"), y = 4309030, x = 747425, cex = 1.6, pos = 4)

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
text(labels = paste("Deer - Fall", "(", labelsFall$FallDetCount[which(labelsFall$Species == "Odocoileus virginianus")],"/", labelsFall$FallDeployCt[1], ")"), y = 4309030, x = 747425, cex = 1.6, pos = 4)
points(gridXY$NAD83_X[idxf], gridXY$NAD83_Y[idxf], cex = 1, pch = 8, col = "red")

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
text(labels = paste("Deer - Winter", "(", labelsWinter$WinDetCount[which(labelsWinter$Species == "Odocoileus virginianus")],"/", labelsWinter$WinterDeployCt[1], ")"), y = 4309030, x = 747425, cex = 1.6, pos = 4)
points(gridXY$NAD83_X[idxw], gridXY$NAD83_Y[idxw], cex = 1, pch = 8, col = "red")

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
text(labels = paste("Deer - Spring", "(", labelsSpring$SprDetCount[which(labelsSpring$Species == "Odocoileus virginianus")],"/", labelsSpring$SpringDeployCt[1], ")"), y = 4309030, x = 747425, cex = 1.6, pos = 4)



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
text(labels = paste("Fox Squirrel - Summer", "(", labelsSum$SumDetCount[which(labelsSum$Species == "Sciurus niger")],"/", labelsSum$SumDeployCt[1], ")"), y = 4309030, x = 747425, cex = 1.6, pos = 4)

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
text(labels = paste("Fox squirrel - Fall", "(", labelsFall$FallDetCount[which(labelsFall$Species == "Sciurus niger")],"/", labelsFall$FallDeployCt[1], ")"), y = 4309030, x = 747425, cex = 1.6, pos = 4)
points(gridXY$NAD83_X[idxf], gridXY$NAD83_Y[idxf], cex = 1, pch = 8, col = "red")

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
text(labels = paste("Fox squirrel - Winter", "(", labelsWinter$WinDetCount[which(labelsWinter$Species == "Sciurus niger")],"/", labelsWinter$WinterDeployCt[1], ")"), y = 4309030, x = 747425, cex = 1.6, pos = 4)
points(gridXY$NAD83_X[idxw], gridXY$NAD83_Y[idxw], cex = 1, pch = 8, col = "red")

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
text(labels = paste("Fox Squirrel - Spring", "(", labelsSpring$SprDetCount[which(labelsSpring$Species == "Sciurus niger")],"/", labelsSpring$SpringDeployCt[1], ")"), y = 4309030, x = 747425, cex = 1.6, pos = 4)

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
text(labels = paste("Gray Squirrel - Summer", "(", labelsSum$SumDetCount[which(labelsSum$Species == "Sciurus carolinensis")],"/", labelsSum$SumDeployCt[1], ")"), y = 4309030, x = 747425, cex = 1.6, pos = 4)

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
text(labels = paste("Gray Squirrel - Fall", "(", labelsFall$FallDetCount[which(labelsFall$Species == "Sciurus carolinensis")],"/", labelsFall$FallDeployCt[1], ")"), y = 4309030, x = 747425, cex = 1.6, pos = 4)
points(gridXY$NAD83_X[idxf], gridXY$NAD83_Y[idxf], cex = 1, pch = 8, col = "red")

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
text(labels = paste("Gray Squirrel - Winter", "(", labelsWinter$WinDetCount[which(labelsWinter$Species == "Sciurus carolinensis")],"/", labelsWinter$WinterDeployCt[1], ")"), y = 4309030, x = 747425, cex = 1.6, pos = 4)
points(gridXY$NAD83_X[idxw], gridXY$NAD83_Y[idxw], cex = 1, pch = 8, col = "red")

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
text(labels = paste("Gray Squirrel - Spring", "(", labelsSpring$SprDetCount[which(labelsSpring$Species == "Sciurus carolinensis")],"/", labelsSpring$SpringDeployCt[1], ")"), y = 4309030, x = 747425, cex = 1.6, pos = 4)


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
text(labels = paste("Raccoon - Summer", "(", labelsSum$SumDetCount[which(labelsSum$Species == "Procyon lotor")],"/", labelsSum$SumDeployCt[1], ")"), y = 4309030, x = 747425, cex = 1.6, pos = 4)

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
text(labels = paste("Raccoon - Fall", "(", labelsFall$FallDetCount[which(labelsFall$Species == "Procyon lotor")],"/", labelsFall$FallDeployCt[1], ")"), y = 4309030, x = 747425, cex = 1.6, pos = 4)
points(gridXY$NAD83_X[idxf], gridXY$NAD83_Y[idxf], cex = 1, pch = 8, col = "red")

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
text(labels = paste("Raccoon - Winter", "(", labelsWinter$WinDetCount[which(labelsWinter$Species == "Procyon lotor")],"/", labelsWinter$WinterDeployCt[1], ")"), y = 4309030, x = 747425, cex = 1.6, pos = 4)
points(gridXY$NAD83_X[idxw], gridXY$NAD83_Y[idxw], cex = 1, pch = 8, col = "red")

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
text(labels = paste("Raccoon - Spring", "(", labelsSpring$SprDetCount[which(labelsSpring$Species == "Procyon lotor")],"/", labelsSpring$SpringDeployCt[1], ")"), y = 4309030, x = 747425, cex = 1.6, pos = 4)

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
text(labels = paste("Black Bear - Summer", "(", labelsSum$SumDetCount[which(labelsSum$Species == "Ursus americanus")],"/", labelsSum$SumDeployCt[1], ")"), y = 4309030, x = 747425, cex = 1.6, pos = 4)

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
text(labels = paste("Black Bear - Fall", "(", labelsFall$FallDetCount[which(labelsFall$Species == "Ursus americanus")],"/", labelsFall$FallDeployCt[1], ")"), y = 4309030, x = 747425, cex = 1.6, pos = 4)
points(gridXY$NAD83_X[idxf], gridXY$NAD83_Y[idxf], cex = 1, pch = 8, col = "red")

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
text(labels = paste("Black Bear - Winter", "(", labelsWinter$WinDetCount[which(labelsWinter$Species == "Ursus americanus")],"/", labelsWinter$WinterDeployCt[1], ")"), y = 4309030, x = 747425, cex = 1.6, pos = 4)
points(gridXY$NAD83_X[idxw], gridXY$NAD83_Y[idxw], cex = 1, pch = 8, col = "red")

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
text(labels = paste("Black Bear - Spring", "(", labelsSpring$SprDetCount[which(labelsSpring$Species == "Ursus americanus")],"/", labelsSpring$SpringDeployCt[1], ")"), y = 4309030, x = 747425, cex = 1.6, pos = 4)
