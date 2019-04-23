#Trying to test how many cameras it takes, on average, randomly selected, to get a mean capture rate estimate within 1 SD of the overal grid average (or something other metric), here basically assumed to be "correct".

load("data/grsqDataSum.RData")
load("data/grsqDataFall.RData")
load("data/grsqDataWin.RData")
load("data/grsqDataSpr.RData")

#Starting with Summer First
gridMean <- mean(grsqDataSum$CR)
gridSD <- sd(grsqDataSum$CR)
gridCRSDHigh <- gridMean + gridSD
gridCRSDLow <- gridMean - gridSD

#The loop below technically works, but doing this 1000 times just ends up with means that are close to the overall grid mean. And in the end the SD of the 27 camera grid is quite high, so getting within 1 SD of the overall grid mean is not a high bar. It even overlaps zero.

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

# A second option proposed by Valentine was to take the 27 cameras, in order, and accumulate mean CR values. First with the first 2 cameras, then the first 3, then the first 4 etc, until you have a vector of means, likely approaching the grid mean. Then rerarrange the 27 and repeat. You could maybe do this something like 50 times, and plot all these lines on a graph, with the overal grid mean line displayed. I don't have any idea what metric to use to know when the line is "close enough" but perhaps it will be visually helpful?