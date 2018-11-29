#This file calculates various summary statistics for the SCBI Camera Grid project. It should call in the clean data object called "camdataAllClean" from the data folder in the git repository

library(reshape)
library(dplyr)
library(tidyverse)
load("data/CamdataAllClean.RData")

alldata <- camdataMSeq #shorter name for ease of coding

#Here I should remove unwanted species and categories right away.
badSP <- c("Animal Not on List", "Camera Misfire", "Camera Trapper", "Corvus brachyrhynchos", "Corvus corax", "Cyanocitta cristata", "Homo sapiens", "Meleagris gallopavo", "No Animal", "Other Bird species", "Owl species", "Raptor Species", "Time Lapse", "Unknown Animal", "Unknown Bird", "Unknown Canid", "Unknown Skunk_Badger")
alldata <- subset(alldata, !alldata$Species.Name %in% badSP)
levels(as.factor(alldata$Species.Name))

#Summarize Sequence Counts and Capture Rates for Whole Grid, All Seasons --------

#create summary table of frequencies of sequences and convert it to a data frame
camdata_summary<-table(alldata$Deployment_Name, 
                             alldata$Species.Name)
camdata_summary<-as.data.frame(camdata_summary)

#change column names
colnames(camdata_summary)<-c("Deployment_Name", "Species", "Number_of_Sequences")

#get unique camera effort values from sequence data for each deployment. This immediately excludes cameras that had no sequences recorded, which is OK.
effort <- camdataMSeq %>%
  select(Deployment_Name, Deploy.Duration) %>%
  distinct
totEffort <- sum(as.numeric(effort$Deploy.Duration))

#Overall Capture Rate per species, across everything (NOT THE MEAN)
CRAll <- as.data.frame(tapply(camdata_summary$Number_of_Sequences, camdata_summary$Species, sum))
names(CRAll)[1] <- "TotalSeqs"
CRAll$TotCR <- (CRAll$TotalSeqs/totEffort)*100

#Merge in the camera effort and calculate capture rate for each species and each deployment, as rate of sequences per 100 camera nights
camdata_summary <- merge(camdata_summary, effort, by = "Deployment_Name" )
camdata_summary$Deploy.Duration <- as.numeric(camdata_summary$Deploy.Duration)
camdata_summary$CR <- (camdata_summary$Number_of_Sequences/camdata_summary$Deploy.Duration) *100

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



# Summary Stats - CapRates All Seasons All Species ------------------------
#calculating average capture rate for each species across the grid. This can be seen as a baseline number against which we can compare numbers at individual cameras. Note that in this case, this is across all 4 seasons, so we are averaging over more than 100 data points.

baselines<- tapply(camdata_summary$CR, camdata_summary$Species, mean)
baselines<- as.data.frame(baselines)
names(baselines)[1] <- "Mean_CR"

#SD of average number of sequences across the grid for each species
baselines_sd <- tapply(camdata_summary$CR, camdata_summary$Species, sd)
baselines$sd<- baselines_sd
names(baselines)[2] <- "StDev_CR"

#Capture rate by species by camera, here with species as rows
CR_bD_bS <- tapply(camdata_summary$CR, list(camdata_summary$Species, camdata_summary$Deployment_Name), mean)
CR_bD_bS <- as.data.frame(CR_bD_bS)

#Min and max CR for each species across the grid
baselines$minCR <- apply(CR_bD_bS,1,min)
names(baselines)[3]<- "Min_CR"
baselines$maxCR <- apply(CR_bD_bS,1,max)
names(baselines)[4]<-"Max_CR"

#Add columns to CR_bD_bS for standard deviation, mean, max and min.
CR_bD_bS$Mean_CR <- baselines$Mean_CR
CR_bD_bS$StDev_CR<-baselines$StDev_CR
CR_bD_bS$Min_CR<-baselines$Min_CR
CR_bD_bS$Max_CR<-baselines$Max_CR

baselines <- merge(baselines, CRAll, by = "row.names")

#At this point I just have baseline summary values for capture rates by species, across all seasons, in a few difffert formats for various uses later.



#Calculate Proportion of Species Captures for Each Camera####


#Create replicate data frame and remove mean,sd,min,max columns
CR_bD_bS2<- dplyr::select(CR_bD_bS,-c("Mean_CR", "StDev_CR", "Min_CR", "Max_CR"))

#convert all values >0 to 1
CR_bD_bS2<-as.data.frame((ifelse(CR_bD_bS2==0,0,1)))

#create column of camera capture proportion per species for whole grid to original data frame
CR_bD_bS$Prop_Grid<-rowSums(CR_bD_bS2)/length(unique(alldata$Deployment_Name))


#Boxplot of CapRates and Proportions of Deployments----


#Transpose headers
CR<-t(CR_bD_bS2)

#reshape data so that each row is per deployment, per species, per capture rate
#FIGURE OUT IF I ALREADY HAVE ALL THIS. MIGHT BE REDUNDANT
CR1<-melt(CR)
colnames(CR1)<-c ("Deployment", "Species", "Detected")
#Calculates proportion of cameras that captured each species
CR2<-subset(CR1, Detected !=0)
NumPres<-as.data.frame(table(CR2$Species))
colnames(NumPres)<-c("Species","NumPres")

#Merge proportion dataframe with rest of data
CR3<-merge(CR1, NumPres, by = "Species")

CR3["tot"]<-length(unique(CR1$Deployment))

#Adds column label which outputs the desired text labels for the boxplot. 
#NEED TO MAKE LABELS FOR OTHER SEASONS HERE. I THINK I HAVE ALL IN NEED in CR3 to DO IT DIRECTLY
CR3["labelAll"]<-paste(CR3$Species,"|", CR3$NumPres, "/", CR3$tot)
Labels_all <- CR3 %>%
  select(Species, labelAll) %>%
  distinct

camdata_summary <- merge(camdata_summary, Labels_all, by = "Species")
Season <- alldata %>%
  select(Deployment_Name, Deployment_Name2, Season) %>%
  distinct
camdata_summary <- merge(camdata_summary, Season, by = "Deployment_Name")

save(camdata_summary, file = "data/camdata_summary")

#Boxplot showing median capture rate per species, all seasons. Remember I can use "subset" here to just do a boxplot for summer, but it much more complicated to get the proportion of stations with presence in this way. Might need to redo all code above by season??
par(mar=c(9,17,4,2))
plot<-boxplot(camdata_summary$CR~camdata_summary$labelAll,
              cex.main = 2.5,
              main = "",
              cex.lab = 1.0,
              cex.axis = 1.0,
              horizontal = T, las = 2, cex.axis = 1,
              names.arg = camdata_summary$Species, at=rank(tapply(camdata_summary$CR,camdata_summary$Species,median), ties.method = "random"))
mtext(expression(bold("Species")), side = 2, line = 15, cex = 1.7)
mtext(expression(bold("Total Capture Rate (events per 100 camera-nights)")), side = 1 , line = 4, cex = 1.3)

#Need to bring in Circle Graphs here



