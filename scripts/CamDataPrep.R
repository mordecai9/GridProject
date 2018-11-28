#Script for formatting and cleaning raw sequence camera trap data from high resolution camera grid established in Posey Hollow SCBI end of May 2017.  Twenty-seven Reconyx cameras approximately 20 meters apart.Cameras run in 2 month seasonal periods and distance markers at 2.5 meters apart up to 10 meters.Photos uploaded to eMammal and data downloaded from eMammal as csv file. The raw data was downloaded and saved as "rawGridData_4S_final.csv". We also kept track of when cameras were running, and when problems occurred. These occur in 4 separate csv files starting with "CameraOperation_". 

#This script removes a few errors in deployment names, converts the date/times into proper format, merges in the camera operation information and uses it to calculate deployment length (i.e. camera nights/effort), removes sequences that occurred when the camera was not functioning properly, removes sequences that started less than 10 minutes from the end of the previous sequence for the same species at the same camera (independence threshold) and then removes unnecessary columns. In the end, we get an R object dataframe called "CamdataAllClean.RData" which should be the starting point for all other scripts.


library(dplyr)


#Bringing in emammal raw data from all deployments. it looks like for some reason we are missing all sequences from 0503_W17. Not sure why. Sent email to Jen Zhao
camdataRaw<-read.csv("data/rawGridData_4S_final.csv", stringsAsFactors = FALSE)
levels(as.factor(camdataRaw$Deployment.Name))

levels(as.factor(camdataRaw$Deployment.Name))
#Substitute Deployment named '0104-SI17-2' to '0104-SI17'
#This was the result of uploading this deployment on two seperate days to eMammal
camdataRaw$Deployment.Name<-sub(pattern = '0104-SI17-2','0104-SI17',camdataRaw$Deployment.Name)

#Substitute Deployment named '0403_F17-1' to '0403_F17'
#This was the result of uploading this deployment on two seperate days to eMammal
camdataRaw$Deployment.Name<-sub(pattern = '0403_F17-1','0403_F17',camdataRaw$Deployment.Name)

#There should now be 27*4 = 108 levels of the "Deployment Name", but no sequences were recorded for 0101_W17 and 0401_W17. Photos were also lost for 0505_F17. So we should end up with 105 levels.
levels(as.factor(camdataRaw$Deployment.Name))

#Convert Begin.Time and End.Time columns to date types
#replace "T" in time columns with a space first
camdataRaw$Begin.Time<-sub(pattern = 'T',' ', camdataRaw$Begin.Time)
camdataRaw$End.Time<-sub(pattern = 'T', ' ', camdataRaw$End.Time)
camdataRaw$Begin.Time<-as.POSIXct(as.character(camdataRaw$Begin.Time, "%Y-%b-%d %H:%M:%S", tz = "EST5EDT"))
camdataRaw$End.Time<-as.POSIXct(as.character(camdataRaw$End.Time, "%Y-%b-%d %H:%M:%S", tz = "EST5EDT"))

#bringing in deployment length information
camnightdata_SI17<-read.csv("data/CameraOperation_SI17.csv", stringsAsFactors = FALSE)
camnightdata_F17<-read.csv("data/CameraOperation_F17.csv", stringsAsFactors = FALSE)
camnightdata_W17<-read.csv("data/CameraOperation_W17.csv", stringsAsFactors = FALSE)
camnightdata_Sp18<-read.csv("data/CameraOperation_Sp18.csv", stringsAsFactors = FALSE)

#Combine all 4 seasons camera operation information
camopsAll <- rbind(camnightdata_SI17, camnightdata_F17, camnightdata_W17, camnightdata_Sp18)

#Convert time in and out date columns to read as dates
#Replace "T" in the time columns with a space first
camopsAll$Date.out<-sub(pattern = 'T',' ', camopsAll$Date.out)
camopsAll$Date.in<-sub(pattern = 'T', ' ', camopsAll$Date.in)
camopsAll$Shift.begin.1<-sub(pattern = 'T', ' ', camopsAll$Shift.begin.1)
camopsAll$Shift.end.1<-sub(pattern = 'T', ' ', camopsAll$Shift.end.1)
camopsAll$Shift.begin.2<-sub(pattern = 'T', ' ', camopsAll$Shift.begin.2)
camopsAll$Shift.end.2<-sub(pattern = 'T', ' ', camopsAll$Shift.end.2)
camopsAll$Date.out<-as.POSIXct(as.character(camopsAll$Date.out,"%Y-%b-%d %H:%M:%S", tz = "EST5EDT"))
camopsAll$Date.in<-as.POSIXct(as.character(camopsAll$Date.in,"%Y-%b-%d %H:%M:%S", tz = "EST5EDT"))
camopsAll$Shift.begin.1<-as.POSIXct(as.character(camopsAll$Shift.begin.1, "%Y-%b-%d %H:%M:%S", tz = "EST5EDT"))
camopsAll$Shift.end.1<-as.POSIXct(as.character(camopsAll$Shift.end.1, "%Y-%b-%d %H:%M:%S", tz = "EST5EDT"))
camopsAll$Shift.begin.2<-as.POSIXct(as.character(camopsAll$Shift.begin.2, "%Y-%b-%d %H:%M:%S", tz = "EST5EDT"))
camopsAll$Shift.end.2<-as.POSIXct(as.character(camopsAll$Shift.end.2, "%Y-%b-%d %H:%M:%S", tz = "EST5EDT"))

#Calculate Deployment Duration by subtracting shift length from deployment length
#Shift means a period of days where the camera was shifted and is unusable data
camopsAll$Deploy.Length<-difftime(camopsAll$Date.in ,camopsAll$Date.out , units = c("days"))
camopsAll$Deploy.Shifts.1<-difftime(camopsAll$Shift.end.1 ,camopsAll$Shift.begin.1 , units = c("days"))
camopsAll$Deploy.Shifts.2<-difftime(camopsAll$Shift.end.2 ,camopsAll$Shift.begin.2 , units = c("days"))
camopsAll$Deploy.Duration<-camopsAll$Deploy.Length-camopsAll$Deploy.Shifts.1-camopsAll$Deploy.Shifts.2

#Substitute periods in column titles to underscores (if needed)
names(camdataRaw) <- sub(pattern = 'Deployment.', 'Deployment_',names(camdataRaw))

#Merge camdata and camera operation information by deployment name. Should not loose any sequences here.
camdataM<-merge(camdataRaw, camopsAll, by = "Deployment_Name", all.x = TRUE)

#remove rows that occur after end dates for deployments
#Need to ensure there are no sequences for 0206 and 0505 for Fall. 0101 for Winter, 0401 for Winter, 
camdataM<-subset(camdataM, Begin.Time>=Date.out&Begin.Time<=Date.in)
camdataM<-subset(camdataM, End.Time<Shift.begin.1|Begin.Time>=Shift.end.1)
camdataM<-subset(camdataM, End.Time<Shift.begin.2|Begin.Time>=Shift.end.2)

#Here is where we would set our threshold of sequence independence. Automatically in emammal it is set to 1 minute. Here I am setting it to 10 minutes. There is a very small chance that this code below, because it does all species together, would remove a sequence in error. I think this would happen only when the first sequence of a species in a deployment occurs within 10 minutes of the last sequence of another species, that comes before it in alphabetical order. This seems pretty unlikely.

camdataM <- camdataM[order(camdataM$Deployment_Name,camdataM$Species.Name, camdataM$Begin.Time), ] #sort rows
camdataM$series <- NA #create new column for time differences


#calculate time difference in minutes for each row, compared to previous row
for (i in 1:length(camdataM$Begin.Time)) 
{camdataM$series[i+1] <- difftime(camdataM$Begin.Time[i+1],camdataM$End.Time[i],units="mins")}

camdataM$event <- NA #create new column, which we will give a 1 or 0 to, based on whether we want to keep it
camdataM$event <- ifelse(camdataM$series <= 10 & camdataM$series >= 0, 0, 1)
camdataM$event[1] <- 1
camdataMSeq <-subset(camdataM,camdataM$event == 1) #keep only rows where the event column is "1"

#Remove unnecessary columns and save final clean R object for use in other scripts
colnames(camdataMSeq)
badcols <- c("Project", "Subproject", "Treatment", "ID.Type", "Age", "Sex", "Individual", "Count", "Actual.Lat", "Actual.Lon", "Fuzzed", "Notes", "Deploy.Length", "Deploy.Shifts.1", "Deploy.Shifts.2", "series", "event" )
camdataMSeq <- dplyr::select(camdataMSeq, -(badcols))

save(camdataMSeq, file = "data/CamdataAllClean.RData")


