#This script is designed to take my cleaned data from the SCBI 1-year testing grid and convert it into usable data for occupancy models. 

#This data file has already been cleaned, in that sequences have already been identified to be those separated by at least 10 minutes. The duration of the deployment has already been calculated based on any potential problem days, and any photos taken during problem events have already been removed. The cleaning happened in the script "CamDataPrep.R".

#Note that cameras 0101_W17, 0401_W17, 0505_F17 and 0206_F17 all failed entirely or were lost. I need to remove these from the camera operation files because no photos should exist in the sequence data from these after the cleaning.

#Here I will use camtrapR to prepare the data into a detection history matrix. I'll need to sequence data, which exists as the object "CamdataAllClean.Rdata", and I'll need the camera operation information, which sits for each season in an excel file.

library(camtrapR)
library(tidyverse)

# Bring in Camera Sequence Information ------------------------------------


load("data/CamdataAllClean.RData")

alldata <- camdataMSeq #shorter name for ease of coding

#This should be done ideally in CamDataPrep.R but oh well.
#Here I remove unwanted species and categories right away.
badSP <- c("Animal Not on List", "Camera Misfire", "Camera Trapper", "Corvus brachyrhynchos", "Corvus corax", "Cyanocitta cristata", "Homo sapiens", "Meleagris gallopavo", "No Animal", "Other Bird species", "Owl species", "Raptor Species", "Time Lapse", "Unknown Animal", "Unknown Bird", "Unknown Canid", "Unknown Skunk_Badger")
alldata <- subset(alldata, !alldata$Species.Name %in% badSP)
levels(as.factor(alldata$Species.Name))

alldata$Species.Name<-sub(pattern = 'Unknown Small Rodent','Unknown small rodent',alldata$Species.Name)
alldata$Species.Name<-sub(pattern = 'Unknown Squirrel','Unknown squirrel',alldata$Species.Name)

#Species ID needs to be called "Species". 
names(alldata)
names(alldata)[10] <- "Species"

#Remove unnecessary columns so camtrapR doesn't get confused
alldata_cut <- alldata %>%
  select(c(Deployment_Name, Begin.Time, Species, Deployment_Name2, Season))
  

# Bring in Camera Operation Information -----------------------------------


#bringing in deployment length information
camnightdata_SI17<-read.csv("data/CameraOperation_SI17.csv", stringsAsFactors = FALSE)
names(camnightdata_SI17)[4:7] <- c("Problem1_from", "Problem1_to", "Problem2_from", "Problem2_to")


camnightdata_F17<-read.csv("data/CameraOperation_F17.csv", stringsAsFactors = FALSE)
camnightdata_F17 <- subset(camnightdata_F17, camnightdata_F17$Deployment_Name != "0505_F17" & camnightdata_F17$Deployment_Name != "0206_F17")
names(camnightdata_F17)[4:7] <- c("Problem1_from", "Problem1_to", "Problem2_from", "Problem2_to")


camnightdata_W17<-read.csv("data/CameraOperation_W17.csv", stringsAsFactors = FALSE)
camnightdata_W17 <- subset(camnightdata_W17, camnightdata_W17$Deployment_Name != "0101_W17" & camnightdata_W17$Deployment_Name != "0401_W17")
names(camnightdata_W17)[4:7] <- c("Problem1_from", "Problem1_to", "Problem2_from", "Problem2_to")


camnightdata_Sp18<-read.csv("data/CameraOperation_Sp18.csv", stringsAsFactors = FALSE)
names(camnightdata_Sp18)[4:7] <- c("Problem1_from", "Problem1_to", "Problem2_from", "Problem2_to")

#Final steps for Summer
#Convert times to proper format in camera operation csv files
#Convert time in and out date columns to read as dates
#Replace "T" in the time columns with a space first

camnightdata_SI17$Date.out<-sub(pattern = 'T',' ', camnightdata_SI17$Date.out)
camnightdata_SI17$Date.in<-sub(pattern = 'T', ' ', camnightdata_SI17$Date.in)
camnightdata_SI17$Problem1_from<-sub(pattern = 'T', ' ', camnightdata_SI17$Problem1_from)
camnightdata_SI17$Problem1_to<-sub(pattern = 'T', ' ', camnightdata_SI17$Problem1_to)
camnightdata_SI17$Problem2_from<-sub(pattern = 'T', ' ', camnightdata_SI17$Problem2_from)
camnightdata_SI17$Problem2_to<-sub(pattern = 'T', ' ', camnightdata_SI17$Problem2_to)
camnightdata_SI17$Date.out<-as.POSIXct(as.character(camnightdata_SI17$Date.out,"%Y-%b-%d %H:%M:%S", tz = "EST5EDT"))
camnightdata_SI17$Date.in<-as.POSIXct(as.character(camnightdata_SI17$Date.in,"%Y-%b-%d %H:%M:%S", tz = "EST5EDT"))
camnightdata_SI17$Problem1_from<-as.POSIXct(as.character(camnightdata_SI17$Problem1_from, "%Y-%b-%d %H:%M:%S", tz = "EST5EDT"))
camnightdata_SI17$Problem1_to<-as.POSIXct(as.character(camnightdata_SI17$Problem1_to, "%Y-%b-%d %H:%M:%S", tz = "EST5EDT"))
camnightdata_SI17$Problem2_from<-as.POSIXct(as.character(camnightdata_SI17$Problem2_from, "%Y-%b-%d %H:%M:%S", tz = "EST5EDT"))
camnightdata_SI17$Problem2_to<-as.POSIXct(as.character(camnightdata_SI17$Problem2_to, "%Y-%b-%d %H:%M:%S", tz = "EST5EDT"))

camopS <- cameraOperation(CTtable = camnightdata_SI17,
                          stationCol = "Deployment_Name",
                          setupCol = "Date.out",
                          retrievalCol = "Date.in",
                          writecsv = FALSE,
                          hasProblems = TRUE,
                          dateFormat = "%Y-%m-%d %H:%M:%S"
)
head(camopS)

#Final steps for Fall
#Convert times to proper format in camera operation csv files
#Convert time in and out date columns to read as dates
#Replace "T" in the time columns with a space first

camnightdata_F17$Date.out<-sub(pattern = 'T',' ', camnightdata_F17$Date.out)
camnightdata_F17$Date.in<-sub(pattern = 'T', ' ', camnightdata_F17$Date.in)
camnightdata_F17$Problem1_from<-sub(pattern = 'T', ' ', camnightdata_F17$Problem1_from)
camnightdata_F17$Problem1_to<-sub(pattern = 'T', ' ', camnightdata_F17$Problem1_to)
camnightdata_F17$Problem2_from<-sub(pattern = 'T', ' ', camnightdata_F17$Problem2_from)
camnightdata_F17$Problem2_to<-sub(pattern = 'T', ' ', camnightdata_F17$Problem2_to)
camnightdata_F17$Date.out<-as.POSIXct(as.character(camnightdata_F17$Date.out,"%Y-%b-%d %H:%M:%S", tz = "EST5EDT"))
camnightdata_F17$Date.in<-as.POSIXct(as.character(camnightdata_F17$Date.in,"%Y-%b-%d %H:%M:%S", tz = "EST5EDT"))
camnightdata_F17$Problem1_from<-as.POSIXct(as.character(camnightdata_F17$Problem1_from, "%Y-%b-%d %H:%M:%S", tz = "EST5EDT"))
camnightdata_F17$Problem1_to<-as.POSIXct(as.character(camnightdata_F17$Problem1_to, "%Y-%b-%d %H:%M:%S", tz = "EST5EDT"))
camnightdata_F17$Problem2_from<-as.POSIXct(as.character(camnightdata_F17$Problem2_from, "%Y-%b-%d %H:%M:%S", tz = "EST5EDT"))
camnightdata_F17$Problem2_to<-as.POSIXct(as.character(camnightdata_F17$Problem2_to, "%Y-%b-%d %H:%M:%S", tz = "EST5EDT"))

#Need to change problem start dates for 0401 because they start before the deployment begins. Looks like NA works, because it just ignores the problem dates columns
camnightdata_F17[16,4:7] <- NA

camopF <- cameraOperation(CTtable = camnightdata_F17,
                          stationCol = "Deployment_Name",
                          setupCol = "Date.out",
                          retrievalCol = "Date.in",
                          writecsv = FALSE,
                          hasProblems = TRUE,
                          dateFormat = "%Y-%m-%d %H:%M:%S"
)
head(camopF)

#Final steps for Winter
#Convert times to proper format in camera operation csv files
#Convert time in and out date columns to read as dates
#Replace "T" in the time columns with a space first

camnightdata_W17$Date.out<-sub(pattern = 'T',' ', camnightdata_W17$Date.out)
camnightdata_W17$Date.in<-sub(pattern = 'T', ' ', camnightdata_W17$Date.in)
camnightdata_W17$Problem1_from<-sub(pattern = 'T', ' ', camnightdata_W17$Problem1_from)
camnightdata_W17$Problem1_to<-sub(pattern = 'T', ' ', camnightdata_W17$Problem1_to)
camnightdata_W17$Problem2_from<-sub(pattern = 'T', ' ', camnightdata_W17$Problem2_from)
camnightdata_W17$Problem2_to<-sub(pattern = 'T', ' ', camnightdata_W17$Problem2_to)
camnightdata_W17$Date.out<-as.POSIXct(as.character(camnightdata_W17$Date.out,"%Y-%b-%d %H:%M:%S", tz = "EST5EDT"))
camnightdata_W17$Date.in<-as.POSIXct(as.character(camnightdata_W17$Date.in,"%Y-%b-%d %H:%M:%S", tz = "EST5EDT"))
camnightdata_W17$Problem1_from<-as.POSIXct(as.character(camnightdata_W17$Problem1_from, "%Y-%b-%d %H:%M:%S", tz = "EST5EDT"))
camnightdata_W17$Problem1_to<-as.POSIXct(as.character(camnightdata_W17$Problem1_to, "%Y-%b-%d %H:%M:%S", tz = "EST5EDT"))
camnightdata_W17$Problem2_from<-as.POSIXct(as.character(camnightdata_W17$Problem2_from, "%Y-%b-%d %H:%M:%S", tz = "EST5EDT"))
camnightdata_W17$Problem2_to<-as.POSIXct(as.character(camnightdata_W17$Problem2_to, "%Y-%b-%d %H:%M:%S", tz = "EST5EDT"))

#Need to change problem start dates for 0205 and 0206 because they start before the deployment begins. Looks like NA works, because it just ignores the problem dates columns
camnightdata_W17[9:10,4:7] <- NA

camopW <- cameraOperation(CTtable = camnightdata_W17,
                          stationCol = "Deployment_Name",
                          setupCol = "Date.out",
                          retrievalCol = "Date.in",
                          writecsv = FALSE,
                          hasProblems = TRUE,
                          dateFormat = "%Y-%m-%d %H:%M:%S"
)
head(camopW)

#Final steps for Spring
#Convert times to proper format in camera operation csv files
#Convert time in and out date columns to read as dates
#Replace "T" in the time columns with a space first

camnightdata_Sp18$Date.out<-sub(pattern = 'T',' ', camnightdata_Sp18$Date.out)
camnightdata_Sp18$Date.in<-sub(pattern = 'T', ' ', camnightdata_Sp18$Date.in)
camnightdata_Sp18$Problem1_from<-sub(pattern = 'T', ' ', camnightdata_Sp18$Problem1_from)
camnightdata_Sp18$Problem1_to<-sub(pattern = 'T', ' ', camnightdata_Sp18$Problem1_to)
camnightdata_Sp18$Problem2_from<-sub(pattern = 'T', ' ', camnightdata_Sp18$Problem2_from)
camnightdata_Sp18$Problem2_to<-sub(pattern = 'T', ' ', camnightdata_Sp18$Problem2_to)
camnightdata_Sp18$Date.out<-as.POSIXct(as.character(camnightdata_Sp18$Date.out,"%Y-%b-%d %H:%M:%S", tz = "EST5EDT"))
camnightdata_Sp18$Date.in<-as.POSIXct(as.character(camnightdata_Sp18$Date.in,"%Y-%b-%d %H:%M:%S", tz = "EST5EDT"))
camnightdata_Sp18$Problem1_from<-as.POSIXct(as.character(camnightdata_Sp18$Problem1_from, "%Y-%b-%d %H:%M:%S", tz = "EST5EDT"))
camnightdata_Sp18$Problem1_to<-as.POSIXct(as.character(camnightdata_Sp18$Problem1_to, "%Y-%b-%d %H:%M:%S", tz = "EST5EDT"))
camnightdata_Sp18$Problem2_from<-as.POSIXct(as.character(camnightdata_Sp18$Problem2_from, "%Y-%b-%d %H:%M:%S", tz = "EST5EDT"))
camnightdata_Sp18$Problem2_to<-as.POSIXct(as.character(camnightdata_Sp18$Problem2_to, "%Y-%b-%d %H:%M:%S", tz = "EST5EDT"))

#Need to change problem start dates for 0205 and 0206 because they start before the deployment begins. Looks like NA works, because it just ignores the problem dates columns
#camnightdata_Sp18[9:10,4:7] <- NA

camopSp <- cameraOperation(CTtable = camnightdata_Sp18,
                          stationCol = "Deployment_Name",
                          setupCol = "Date.out",
                          retrievalCol = "Date.in",
                          writecsv = FALSE,
                          hasProblems = TRUE,
                          dateFormat = "%Y-%m-%d %H:%M:%S"
)
head(camopSp)


species <-
  c(
    "Odocoileus virginianus",
    "Sciurus carolinensis",
    "Ursus americanus" ,
    "Procyon lotor" ,
    "Sciurus niger"
  )

# Detection History Creation for Summer ------------------------------------------

#Loop to make detection history for each of 5 species in Summer
alldata_cutS <- alldata_cut %>%
  filter(Season == "Summer 2017")


for (i in species) {

#Code for generation of Detection History
  #Instead of the temp and save I could write to csv and load it back in later I guess
DetHist_temp <- detectionHistory(recordTable = alldata_cutS,
                                 camOp = camopS,
                                 stationCol = "Deployment_Name",
                                 speciesCol = "Species",
                                 recordDateTimeCol = "Begin.Time",
                                 species = i,
                                 occasionLength = 1,
                                 day1 = "station",
                                 datesAsOccasionNames = FALSE,
                                 includeEffort = FALSE,
                                 occasionStartTime = 12,
                                 timeZone = "US/Eastern",
                                 writecsv = FALSE
)
temp <- DetHist_temp$detection_history
save(temp, file = paste("results/DetHistSum",i, sep = "")) 
}

#Loading them back in to check them
load("results/DetHistSumSciurus niger")
DHSumSN <- temp
rowSums(DHSumSN, na.rm = T) #checks out with other data summaries

load("results/DetHistSumSciurus carolinensis")
DHSumSC <- temp
rowSums(DHSumSC, na.rm = T) #checks out with other data summaries

load("results/DetHistSumOdocoileus virginianus")
DHSumOV <- temp
rowSums(DHSumOV, na.rm = T) #checks out with other data summaries

load("results/DetHistSumUrsus americanus")
DHSumUA <- temp
rowSums(DHSumUA, na.rm = T) #checks out with other data summaries


load("results/DetHistSumProcyon lotor")
DHSumPL <- temp
rowSums(DHSumPL, na.rm = T) #checks out with other data summaries


# Detection Histories for Fall ------------------------------------------

alldata_cutF <- alldata_cut %>%
  filter(Season == "Fall 2017")

#Code for generation of Detection History
for (i in species) {
  
  #Code for generation of Detection History
  #Instead of the temp and save I could write to csv and load it back in later I guess
  DetHist_temp <- detectionHistory(recordTable = alldata_cutF,
                                   camOp = camopF,
                                   stationCol = "Deployment_Name",
                                   speciesCol = "Species",
                                   recordDateTimeCol = "Begin.Time",
                                   species = i,
                                   occasionLength = 1,
                                   day1 = "station",
                                   datesAsOccasionNames = FALSE,
                                   includeEffort = FALSE,
                                   occasionStartTime = 12,
                                   timeZone = "US/Eastern",
                                   writecsv = FALSE
  )
  temp <- DetHist_temp$detection_history
  save(temp, file = paste("results/DetHistFall",i, sep = "")) 
}

#Loading them back in to check them
load("results/DetHistFallSciurus niger")
DHFallSN <- temp
rowSums(DHFallSN, na.rm = T) #checks out with other data summaries

load("results/DetHistFallSciurus carolinensis")
DHFallSC <- temp
rowSums(DHFallSC, na.rm = T) #checks out with other data summaries

load("results/DetHistFallOdocoileus virginianus")
DHFallOV <- temp
rowSums(DHFallOV, na.rm = T) #checks out with other data summaries

load("results/DetHistFallUrsus americanus")
DHFallUA <- temp
rowSums(DHFallUA, na.rm = T) #checks out with other data summaries

load("results/DetHistFallProcyon lotor")
DHFallPL <- temp
rowSums(DHFallPL, na.rm = T) #checks out with other data summaries


# Detection Histories for Winter ------------------------------------------

alldata_cutW <- alldata_cut %>%
  filter(Season == "Winter 2017")

#Code for generation of Detection History
for (i in species) {
  
  #Code for generation of Detection History
  #Instead of the temp and save I could write to csv and load it back in later I guess
  DetHist_temp <- detectionHistory(recordTable = alldata_cutW,
                                   camOp = camopW,
                                   stationCol = "Deployment_Name",
                                   speciesCol = "Species",
                                   recordDateTimeCol = "Begin.Time",
                                   species = i,
                                   occasionLength = 1,
                                   day1 = "station",
                                   datesAsOccasionNames = FALSE,
                                   includeEffort = FALSE,
                                   occasionStartTime = 12,
                                   timeZone = "US/Eastern",
                                   writecsv = FALSE
  )
  temp <- DetHist_temp$detection_history
  save(temp, file = paste("results/DetHistWin",i, sep = "")) 
}

#Loading them back in to check them
load("results/DetHistWinSciurus niger")
DHWinSN <- temp
rowSums(DHWinSN, na.rm = T) #checks out with other data summaries

load("results/DetHistWinSciurus carolinensis")
DHWinSC <- temp
rowSums(DHWinSC, na.rm = T) #checks out with other data summaries

load("results/DetHistWinOdocoileus virginianus")
DHWinOV <- temp
rowSums(DHWinOV, na.rm = T) #checks out with other data summaries

load("results/DetHistWinUrsus americanus")
DHWinUA <- temp
rowSums(DHWinUA, na.rm = T) #checks out with other data summaries

load("results/DetHistWinProcyon lotor")
DHWinPL <- temp
rowSums(DHWinPL, na.rm = T) #checks out with other data summaries

# Detection Histories for Spring ------------------------------------------

alldata_cutSp <- alldata_cut %>%
  filter(Season == "Spring 2018")

#Code for generation of Detection History
for (i in species) {
  
  #Code for generation of Detection History
  #Instead of the temp and save I could write to csv and load it back in later I guess
  DetHist_temp <- detectionHistory(recordTable = alldata_cutSp,
                                   camOp = camopSp,
                                   stationCol = "Deployment_Name",
                                   speciesCol = "Species",
                                   recordDateTimeCol = "Begin.Time",
                                   species = i,
                                   occasionLength = 1,
                                   day1 = "station",
                                   datesAsOccasionNames = FALSE,
                                   includeEffort = FALSE,
                                   occasionStartTime = 12,
                                   timeZone = "US/Eastern",
                                   writecsv = FALSE
  )
  temp <- DetHist_temp$detection_history
  save(temp, file = paste("results/DetHistSp",i, sep = "")) 
}

#Loading them back in to check them
load("results/DetHistSpSciurus niger")
DHSpSN <- temp
rowSums(DHSpSN, na.rm = T) #checks out with other data summaries

load("results/DetHistSpSciurus carolinensis")
DHSpSC <- temp
rowSums(DHSpSC, na.rm = T) #checks out with other data summaries

load("results/DetHistSpOdocoileus virginianus")
DHSpOV <- temp
rowSums(DHSpOV, na.rm = T) #checks out with other data summaries

load("results/DetHistSpUrsus americanus")
DHSpUA <- temp
rowSums(DHSpUA, na.rm = T) #checks out with other data summaries

load("results/DetHistSpProcyon lotor")
DHSpPL <- temp
rowSums(DHSpPL, na.rm = T) #checks out with other data summaries


# Detection Histories for Pooled Squirrel Data Summer ----------------------------

AS_alldata_cutS <- alldata_cutS %>%
  filter(Species == "Sciurus carolinensis" |
           Species == "Sciurus niger" |
           Species == "Unknown squirrel") %>%
  mutate(Species = "squirrel")

DHSumSQ <- detectionHistory(recordTable = AS_alldata_cutS,
                                 camOp = camopS,
                                 stationCol = "Deployment_Name",
                                 speciesCol = "Species",
                                 recordDateTimeCol = "Begin.Time",
                                 species = "squirrel",
                                 occasionLength = 1,
                                 day1 = "station",
                                 datesAsOccasionNames = FALSE,
                                 includeEffort = FALSE,
                                 occasionStartTime = 12,
                                 timeZone = "US/Eastern",
                                 writecsv = FALSE)
DHSumSQ <- DHSumSQ$detection_history
save(DHSumSQ, file = "results/DetHistSum_AllSquirrel") 

# Detection Histories for Pooled Squirrel Data Fall ----------------------------

AS_alldata_cutF <- alldata_cutF %>%
  filter(Species == "Sciurus carolinensis" |
           Species == "Sciurus niger" |
           Species == "Unknown squirrel") %>%
  mutate(Species = "squirrel")

DHFallSQ <- detectionHistory(recordTable = AS_alldata_cutF,
                            camOp = camopF,
                            stationCol = "Deployment_Name",
                            speciesCol = "Species",
                            recordDateTimeCol = "Begin.Time",
                            species = "squirrel",
                            occasionLength = 1,
                            day1 = "station",
                            datesAsOccasionNames = FALSE,
                            includeEffort = FALSE,
                            occasionStartTime = 12,
                            timeZone = "US/Eastern",
                            writecsv = FALSE)
DHFallSQ <- DHFallSQ$detection_history
save(DHFallSQ, file = "results/DetHistFall_AllSquirrel") 

# Detection Histories for Pooled Squirrel Data Winter ----------------------------

AS_alldata_cutW <- alldata_cutW %>%
  filter(Species == "Sciurus carolinensis" |
           Species == "Sciurus niger" |
           Species == "Unknown squirrel") %>%
  mutate(Species = "squirrel")

DHWinSQ <- detectionHistory(recordTable = AS_alldata_cutW,
                             camOp = camopW,
                             stationCol = "Deployment_Name",
                             speciesCol = "Species",
                             recordDateTimeCol = "Begin.Time",
                             species = "squirrel",
                             occasionLength = 1,
                             day1 = "station",
                             datesAsOccasionNames = FALSE,
                             includeEffort = FALSE,
                             occasionStartTime = 12,
                             timeZone = "US/Eastern",
                             writecsv = FALSE)
DHWinSQ <- DHWinSQ$detection_history
save(DHWinSQ, file = "results/DetHistWin_AllSquirrel") 

# Detection Histories for Pooled Squirrel Data Spring ----------------------------

AS_alldata_cutSp <- alldata_cutSp %>%
  filter(Species == "Sciurus carolinensis" |
           Species == "Sciurus niger" |
           Species == "Unknown squirrel") %>%
  mutate(Species = "squirrel")

DHSpSQ <- detectionHistory(recordTable = AS_alldata_cutSp,
                            camOp = camopSp,
                            stationCol = "Deployment_Name",
                            speciesCol = "Species",
                            recordDateTimeCol = "Begin.Time",
                            species = "squirrel",
                            occasionLength = 1,
                            day1 = "station",
                            datesAsOccasionNames = FALSE,
                            includeEffort = FALSE,
                            occasionStartTime = 12,
                            timeZone = "US/Eastern",
                            writecsv = FALSE)
DHSpSQ <- DHSpSQ$detection_history
save(DHSpSQ, file = "results/DetHistSpring_AllSquirrel") 
