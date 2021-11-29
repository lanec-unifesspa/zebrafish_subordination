#Load necessary libraries to download data (RCurl), run POSET tests (coin), and produce plots (dabestr)
library(dabestr)
library(RCurl)
library(coin)

#Set seed for randomization analyses
set.seed(42)

#Load data
x1 <- getURL("https://raw.githubusercontent.com/lanec-unifesspa/zebrafish_subordination/main/NTT.csv")
NTT <- read.csv(text = x1)
View(NTT)

#Plot data
NTT_TB <- NTT %>% dabest(Status, TimeBottomNTT, idx = list(c("FemaleControl", "FemaleDominant", "FemaleSubordinate"), c("MaleControl", "MaleDominant", "MaleSubordinate")), paired = FALSE) #Organize data in tidy format
multi.group.mean_TB <- NTT_TB %>% mean_diff() #Compute mean differences between groups, with controls as references
plot(multi.group.mean_TB, color.column = Sex, rawplot.ylim = c(0, 360), rawplot.ylabel = "Time spent on bottom third (s)") #Create a Gardner-Altman plot for the data

NTT_TT <- NTT %>% dabest(Status, TimeTopNTT, idx = list(c("FemaleControl", "FemaleDominant", "FemaleSubordinate"), c("MaleControl", "MaleDominant", "MaleSubordinate")), paired = FALSE)
multi.group.mean_TT <- NTT_TT %>% mean_diff()
plot(multi.group.mean_TT, color.column = Sex, rawplot.ylim = c(0, 360), rawplot.ylabel = "Time spent on top third (s)")

NTT_ES <- NTT %>% dabest(Status, AbsoluteTurnAngle, idx = list(c("FemaleControl", "FemaleDominant", "FemaleSubordinate"), c("MaleControl", "MaleDominant", "MaleSubordinate")), paired = FALSE)
multi.group.mean_ES <- NTT_ES %>% mean_diff()
plot(multi.group.mean_ES, color.column = Sex, rawplot.ylim = c(0, 90), rawplot.ylabel = "Absolute turn angle (º)")

NTT_FR <- NTT %>% dabest(Status, Freezing, idx = list(c("FemaleControl", "FemaleDominant", "FemaleSubordinate"), c("MaleControl", "MaleDominant", "MaleSubordinate")), paired = FALSE)
multi.group.mean_FR <- NTT_FR %>% mean_diff()
plot(multi.group.mean_FR, color.column = Sex, rawplot.ylim = c(0, 360), rawplot.ylabel = "Freezing (s)")

NTT_SP <- NTT %>% dabest(Status, Speed, idx = list(c("FemaleControl", "FemaleDominant", "FemaleSubordinate"), c("MaleControl", "MaleDominant", "MaleSubordinate")), paired = FALSE)
multi.group.mean_SP <- NTT_SP %>% mean_diff()
plot(multi.group.mean_SP, color.column = Sex, rawplot.ylim = c(0, 20), rawplot.ylabel = "Swimming speed (cm/s)")

#Define coherence criterion for POSET test
coherence <- function(data) {
  x <- as.matrix(data)
  matrix(apply(x, 1, function(y)
    sum(colSums(t(x) < y) == ncol(x)) -
      sum(colSums(t(x) > y) == ncol(x))), ncol = 1)
}

#Asymptotic POSET test
poset <- independence_test(TimeBottomNTT + AbsoluteTurnAngle + Freezing + Speed ~ Sex:Status, data = NTT, ytraf = coherence)
poset

# Step-down adjusted p-values
pvalue(poset, method = "step-down")