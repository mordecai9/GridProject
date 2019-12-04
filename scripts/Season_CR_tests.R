# Seasonal Comparisons of Capture Rates -----------------------------------

load("data/camdata_summary")

library(MASS)
library(AICcmodavg)



# Correlations of CR at cameras, across seasons ---------------------------
#Here I'm highlighting pairs with correlation values above 0.50
#Bears
camdata_summary %>%
  dplyr::select(Deployment_Name2, CR, Season, Species) %>%
  filter(Species == "Ursus americanus") %>%
  spread(key = Season, value = CR) %>%
  dplyr::select(3:6) %>%
  cor(use = "complete.obs")

#Nothing above 0.50 for bears

#Deer
camdata_summary %>%
  dplyr::select(Deployment_Name2, CR, Season, Species) %>%
  filter(Species == "Odocoileus virginianus") %>%
  spread(key = Season, value = CR) %>%
  dplyr::select(3:6) %>%
  cor(use = "complete.obs")

#Correlated seasons are Fall and Winter (0.79) and Summer and Winter (0.53) but not Fall and Summer (0.48) but it is close.

#Plotting fall vs. Winter, highest correlation out of all species and speason
camdata_summary %>%
  dplyr::select(Deployment_Name2, CR, Season, Species) %>%
  filter(Species == "Odocoileus virginianus") %>%
  spread(key = Season, value = CR) %>%
  dplyr::select(4,5) %>%
  plot()
  

#Raccoon
camdata_summary %>%
  dplyr::select(Deployment_Name2, CR, Season, Species) %>%
  filter(Species == "Procyon lotor") %>%
  spread(key = Season, value = CR) %>%
  dplyr::select(3:6) %>%
  cor(use = "complete.obs")
#Nothing above 0.50 for Raccoons

#Fox squirrel
camdata_summary %>%
  dplyr::select(Deployment_Name2, CR, Season, Species) %>%
  filter(Species == "Sciurus niger") %>%
  spread(key = Season, value = CR) %>%
  dplyr::select(3:6) %>%
  cor(use = "complete.obs")

#Summer and fall are correlated (0.57) and Fall and Winter (0.64) but not summer and winter (0.30) so you can't pool them all. 

#Gray squirrel
camdata_summary %>%
  dplyr::select(Deployment_Name2, CR, Season, Species) %>%
  filter(Species == "Sciurus carolinensis") %>%
  spread(key = Season, value = CR) %>%
  dplyr::select(3:6) %>%
  cor(use = "complete.obs") 

#Fall and Spring (0.54) and Fall and Winter (0.58) and Spring and Winter (0.66) and Summer and Winter(0.52). Only bad one is Summer and Fall (0.36). 

camdata_summary %>%
  dplyr::select(Deployment_Name2, CR, Season, Species) %>%
  filter(Species == "Sciurus carolinensis") %>%
  spread(key = Season, value = CR) %>%
  dplyr::select(4,6) %>%
  plot()

  
# GLM NB test of Season effects on CR -------------------------------------

#Note that this analysis does not test for seasonal change at specific camreas, it does not take camera ID into account, and so is only really testing overall grid mean values

#Deer Season Test
deer_season_test <- camdata_summary %>%
  filter(Species == "Odocoileus virginianus") %>%
  glm.nb(data = ., Number_of_Sequences ~ Season + offset(log(Deploy.Duration)) )

deer_null <- camdata_summary %>%
  filter(Species == "Odocoileus virginianus") %>%
  glm.nb(data = ., Number_of_Sequences ~ 1 + offset(log(Deploy.Duration)) )

SeasonResults <- aictab(cand.set = list(deer_season_test, deer_null), modnames = c("White-tailed deer - Season", "White-tailed deer - Intercept"))

#Fox Squirrel Season Test
foxSq_season_test <- camdata_summary %>%
  filter(Species == "Sciurus niger") %>%
  glm.nb(data = ., Number_of_Sequences ~ Season + offset(log(Deploy.Duration)) )

foxSq_null <- camdata_summary %>%
  filter(Species == "Sciurus niger") %>%
  glm.nb(data = ., Number_of_Sequences ~ 1 + offset(log(Deploy.Duration)) )

SeasonResults <- rbind(SeasonResults, aictab(cand.set = list(foxSq_season_test, foxSq_null), modnames = c("Fox squirrel - Season", "Fox squirrel - Intercept")))

#Gray Squirrel Season Test
graySq_season_test <- camdata_summary %>%
  filter(Species == "Sciurus carolinensis") %>%
  glm.nb(data = ., Number_of_Sequences ~ Season + offset(log(Deploy.Duration)) )

graySq_null <- camdata_summary %>%
  filter(Species == "Sciurus carolinensis") %>%
  glm.nb(data = ., Number_of_Sequences ~ 1 + offset(log(Deploy.Duration)) )

SeasonResults <- rbind(SeasonResults, aictab(cand.set = list(graySq_season_test, graySq_null), modnames = c("Gray squirrel - Season", "Gray squirrel - Intercept")))

#Combined Squirrels. I also create this combined squirrel category in the RegFile Creation Script, so this is redundant, but I didn't feel like cleaning it all up.

foxsqrlData<-subset(camdata_summary, Species == "Sciurus niger")
grsqrlData<-subset(camdata_summary, Species == "Sciurus carolinensis")
unknownsqrlData<-subset(camdata_summary, Species == "Unknown squirrel")

sqrlData<- foxsqrlData %>%
  left_join(grsqrlData, by = "Deployment_Name") %>%
  left_join(unknownsqrlData, by = "Deployment_Name" )

sqrlData <- sqrlData %>%
  mutate(allSqSeqs = Number_of_Sequences.x + Number_of_Sequences.y + Number_of_Sequences) %>%
  dplyr::select(c(allSqSeqs, Season, Deploy.Duration, Species)) %>%
  mutate(Species = "AllSquirrels")

AllSq_season_test <- sqrlData %>%
  glm.nb(data = ., allSqSeqs ~ Season + offset(log(Deploy.Duration)) )

AllSq_null <- sqrlData %>%
  glm.nb(data = ., allSqSeqs ~ 1 + offset(log(Deploy.Duration)) )

SeasonResults <- rbind(SeasonResults, aictab(cand.set = list(AllSq_season_test, AllSq_null), modnames = c("All squirrels - Season", "All squirrels - Intercept")))

#Raccoon Season Test
raccoon_season_test <- camdata_summary %>%
  filter(Species == "Procyon lotor") %>%
  glm.nb(data = ., Number_of_Sequences ~ Season + offset(log(Deploy.Duration)) )

raccoon_null <- camdata_summary %>%
  filter(Species == "Procyon lotor") %>%
  glm.nb(data = ., Number_of_Sequences ~ 1 + offset(log(Deploy.Duration)) )

SeasonResults <- rbind(SeasonResults, aictab(cand.set = list(raccoon_season_test, raccoon_null), modnames = c("Northern raccoon - Season", "Northern raccoon - Intercept")))


#Black Bear Season Test
bear_season_test <- camdata_summary %>%
  filter(Species == "Ursus americanus") %>%
  glm.nb(data = ., Number_of_Sequences ~ Season + offset(log(Deploy.Duration)) )

bear_null <- camdata_summary %>%
  filter(Species == "Ursus americanus") %>%
  glm.nb(data = ., Number_of_Sequences ~ 1 + offset(log(Deploy.Duration)) )

SeasonResults <- rbind(SeasonResults, aictab(cand.set = list(bear_season_test, bear_null), modnames = c("Black bear - Season", "Black bear - Intercept")))

SeasonResultsF <- SeasonResults %>%
  mutate(AICc = round(AICc, 2),
         Delta_AICc = round(Delta_AICc, 2),
         AICcWt = round(AICcWt, 2),
         LL = round(LL, 2)) %>%
  rename(Model = Modnames) %>%
  dplyr::select(-c(Cum.Wt,ModelLik))

#Used for Table S1 in manuscript
#write.csv(SeasonResultsF, file = "results/SeasonResultsTable.csv")

#Unknown Squirrel Season Test
unknownSq_season_test <- camdata_summary %>%
  filter(Species == "Unknown squirrel") %>%
  glm.nb(data = ., Number_of_Sequences ~ Season + offset(log(Deploy.Duration)) )

unknownSq_null <- camdata_summary %>%
  filter(Species == "Unknown squirrel") %>%
  glm.nb(data = ., Number_of_Sequences ~ 1 + offset(log(Deploy.Duration)) )

#Small Rodent Season Test
rodent_season_test <- camdata_summary %>%
  filter(Species == "Unknown small rodent") %>%
  glm.nb(data = ., Number_of_Sequences ~ Season + offset(log(Deploy.Duration)) )

rodent_null <- camdata_summary %>%
  filter(Species == "Unknown small rodent") %>%
  glm.nb(data = ., Number_of_Sequences ~ 1 + offset(log(Deploy.Duration)) )

aictab(cand.set = list(rodent_season_test, rodent_null), modnames = c("rodent Season", "rodent Intercept"))
summary(rodent_season_test)



#Chipmunk Season Test
chip_season_test <- camdata_summary %>%
  filter(Species == "Tamias striatus") %>%
  glm.nb(data = ., Number_of_Sequences ~ Season + offset(log(Deploy.Duration)) )

chip_null <- camdata_summary %>%
  filter(Species == "Tamias striatus") %>%
  glm.nb(data = ., Number_of_Sequences ~ 1 + offset(log(Deploy.Duration)) )

aictab(cand.set = list(chip_season_test, chip_null), modnames = c("Chipmunk Season", "Chipmunk Intercept"))





#Repeated measures approach to the seasonal tests of capture rates####
library(lme4)


#Deer
deer_season_test2 <- camdata_summary_temp %>%
  filter(Species == "Odocoileus virginianus") %>%
  glmer.nb(data = ., Number_of_Sequences ~ Season + offset(log(Deploy.Duration)) + (1|Deployment_Name2) )

deer_null2 <- camdata_summary_temp %>%
  filter(Species == "Odocoileus virginianus") %>%
  glmer.nb(data = ., Number_of_Sequences ~ 1 + offset(log(Deploy.Duration)) + (1|Deployment_Name2) )

SeasonResults <- aictab(cand.set = list(deer_season_test, deer_null), modnames = c("White-tailed deer - Season", "White-tailed deer - Intercept"))
