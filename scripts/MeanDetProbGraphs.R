#Graphing Average Detection Probabilities for Species across seasons. values were calculated by using the best model is each season for each species, and using mean values for any covariates, with no log present if log was a covariate. Tables with final values were generated in each species occupancy script (e.g. OccMods_Deer.R) and final tables were saved as R objects

# Clear Environment and Load Packages -------------------------------------

rm(list = ls())
library(tidyverse)


# Load individual Species Summary Tables------------------------------------------------

load("results/MeanPSeasons")
load("results/MeanPSeasons_AllSq")
load("results/MeanPSeasons_FoxSq")
load("results/MeanPSeasons_GrSq")
load("results/MeanPSeasons_Rac")
load("results/MeanPSeasonsBear")

DetTable <- rbind(MeanpDeer,MeanpSq, MeanpFoxSq, MeanpGrSq, MeanpRac,MeanpBear)
Species <- rep(c("Deer", "All Squirrel", "Fox Squirrel", "Gray Squirrel", "Raccoon", "Bear"), each = 4)
DetTable$Species <- Species[1:22]
Season <- rep(c("Spring", "Winter", 'Summer', "Fall"), times = 5)
DetTable$Season <- c(Season, "Summer", "Fall")

#Create Graph
#using code here from GridSummaryStats.R where I made cap rate comparison graph
library(ggplot2)
DetTable$Species <- factor(DetTable$Species, 
                            levels = c("Bear", "Raccoon", "Gray Squirrel", "Fox Squirrel", "All Squirrel","Deer"))


a <-
  ggplot(data = DetTable, aes(y = est, x = Species, color = Season)) +
  geom_point(shape = 1.4,
             size = 2.5,
             position = position_dodge(.2)) +
  geom_errorbar(
    aes(ymin = lower, ymax = upper),
    size = 1,
    width = 0,
    position = position_dodge(.2))+
  coord_flip() +
  theme_classic() +
  ylab("Mean detection probability") +
  theme(panel.grid.major.y = element_line(color="gray", size = .5, linetype = "dashed"),
        panel.grid.minor.y = element_line(color="gray", size = .5,linetype = "dashed"),
        legend.position = "right")


a
#ggsave("results/SeasonDetProb.tiff", width = 5, height = 4, units = "in")
