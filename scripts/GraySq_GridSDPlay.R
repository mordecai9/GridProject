#Trying to test how many cameras it takes, on average, randomly selected, to get a mean capture rate estimate within 1 SD of the overal grid average, here basically assumed to be "correct".

load("data/grsqDataSum.RData")
load("data/grsqDataFall.RData")
load("data/grsqDataWin.RData")
load("data/grsqDataSpr.RData")

#Starting with Summer First
gridMean <- mean(grsqDataSum$CR)
gridSD <- sd(grsqDataSum$CR)
gridCRSDHigh <- gridMean + gridSD
gridCRSDLow <- gridMean - gridSD

df <- data.frame("Number" = 1:27, "MeanEstCR" = NA, "SD_EstCR" = NA)
temp <- data.frame("Number" = 1:1000, "means" = NA)

for (i in 1:27) {
  
  for (j in 1:1000) {
    
  r <- round(runif(n= i, min = 1, max = 27));
  temp$means[j] <- mean(grsqDataSum$CR[c(r)]);
  
  
  }
  df$MeanEstCR[i] <- mean(temp$means)
  df$SD_EstCR[i] <- sd(temp$means)
}