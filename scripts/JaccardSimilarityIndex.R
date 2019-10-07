#Jaccard Similarity Index for Species Detection History

library(jaccard)

##Load species files for Spring
load("results/DetHistSpSciurus niger")
DHSpSN <- temp

load("results/DetHistSpSciurus carolinensis")
DHSpSC <- temp

load("results/DetHistSpOdocoileus virginianus")
DHSpOV <- temp

load("results/DetHistSpUrsus americanus")
DHSpUA <- temp

load("results/DetHistSpProcyon lotor")
DHSpPL <- temp

##Load species files for Summer
load("results/DetHistSumSciurus niger")
DHSumSN <- temp

load("results/DetHistSumSciurus carolinensis")
DHSumSC <- temp

load("results/DetHistSumOdocoileus virginianus")
DHSumOV <- temp

load("results/DetHistSumUrsus americanus")
DHSumUA <- temp

load("results/DetHistSumProcyon lotor")
DHSumPL <- temp

##Load species for Fall
load("results/DetHistFallSciurus niger")
DHFallSN <- temp

load("results/DetHistFallSciurus carolinensis")
DHFallSC <- temp

load("results/DetHistFallOdocoileus virginianus")
DHFallOV <- temp

load("results/DetHistFallUrsus americanus")
DHFallUA <- temp

load("results/DetHistFallProcyon lotor")
DHFallPL <- temp

##Load species for Winter
load("results/DetHistWinSciurus niger")
DHWinSN <- temp

load("results/DetHistWinSciurus carolinensis")
DHWinSC <- temp

load("results/DetHistWinOdocoileus virginianus")
DHWinOV <- temp

load("results/DetHistWinUrsus americanus")
DHWinUA <- temp

load("results/DetHistWinProcyon lotor")
DHWinPL <- temp

#Changing Temp Files to Data Frames
DFDHFALLOV<-as.data.frame(t(DHFallOV))

##Jaccard Index between Cameras (same season and species)

library('clusteval')
cluster_similarity(DFDHFALLOV$`0101_F17` , DFDHFALLOV$`0102_F17`, similarity="jaccard", method="independence")


library(dplyr)
library(magrittr)
jaccard <- function(DFDHFALLOV, margin=2) {
  if ( margin == 2) {
    M_00 <- apply(DFDHFALLOV, margin, sum) == 0
    M_11 <- apply(DFDHFALLOV, margin, sum) == 2
    if (margin == 2) {
      DFDHFALLOV <- DFDHFALLOV[!M_00, ]
      JSim <- sum(M_11) / ncol(DFDHFALLOV)
    } else {
      
    }
    JDist <- 2 - JSim
    return(c(JSim = JSim, JDist = JDist))
  } else break
}

jaccard(DFDHFALLOV[1:2,], margin=2)


jaccard_per_column <- function(DFDHFALLOV, margin=2){
  require(magrittr)
  require(dplyr)
  key_pairs <- expand.grid(colnames(DFDHFALLOV), colnames(DFDHFALLOV))
  results <- t(apply(key_pairs, 2, function(column) jaccard(DFDHFALLOV[c(column[1], column[2]),], margin=margin)))
  key_pair <- key_pairs %>% mutate(pair = paste(Var1,"_",Var2,sep=""))
  results <- data.frame(results)
  colnames(results) <- key_pair$pair
  results
}

jaccard_per_column(DFDHFALLOV, margin=2)

