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

#Testing one way of calculating Jaccard
#seems this can only be used when the length of the data is the same
#Most of the results from different seasons have different lengths
library('clusteval')
cluster_similarity(DHWinSN[c("0102_W17"),] ,DHWinSN[c("0103_W17"),], similarity="jaccard", method="independence")
