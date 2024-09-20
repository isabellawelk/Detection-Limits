
# Code developed by Bryn Chaffee 
# 
# Intent is to process data from our Metrohm Ion Chromatographs in order to return detection limit and error
#
# Detection limit equation: 
# DL = 3.3 * sigma / m
#  sigma = standard deviation of the residuals
#  m = slope of calibration curve
#  residuals = difference in observed results from linear model
#
# Error is calculated as the mean difference between residuals (duplicates)

##clear workspace
rm(list = ls()) 

##set working directory
setwd("C:/Users/bryn_/OneDrive/Documents/5108 LAB/Data_analysis_detect_limit/IC_Detection-Limits")
list.files()

#load libraries
library(dplyr)
library(reshape)

#import needed .csv files
anions <- read.csv("Run data/DCEW001-589_Anions.csv")   #anion run file
cations <- read.csv("Run data/DCEW001-589_Cations.csv")   #cation run file

an_calibration <- read.csv("Run data/7Jun2024_Anion_Calibration.csv")     #import calibration used for current anion samples
cat_calibration <- read.csv("Run data/29May2024_Cation_Calibration.csv")  #import calibration used for current cation samples

#specify names
an.names <- c("Fluoride","Chloride","Nitrite","Nitrate","Phosphate","Sulfate")
cat.names <- c("Lithium","Sodium","Ammonium","Potassium","Magnesium","Calcium")

##clean the data to remove standards, blanks, and invalid values
an.clean <- anions[which(anions$Sample.type == "Sample"),]    #checks the sample.type column and pulls only "Sample" values
cat.clean <- cations[which(cations$Sample.type == "Sample"),]

an.clean <- replace(an.clean, an.clean == 'invalid', NA)   #replaces any 'invalid' values with NA
cat.clean <- replace(cat.clean, cat.clean == 'invalid', NA)

#remove '-duplicate' from any sample names
an.clean$Info.1 <- gsub("-duplicate", "", as.character(an.clean$Info.1))
an.clean$Info.1 <- gsub("-Duplicate", "", as.character(an.clean$Info.1))
an.clean$Info.1 <- gsub("-DUPLICATE", "", as.character(an.clean$Info.1))
colnames(an.clean) <- c("Sample.type","Info.1","Fluoride Concentration","Chloride Concentration","Nitrite Concentration",
                        "Nitrate Concentration","Phosphate Concentration","Sulfate Concentration")

cat.clean$Info.1 <- gsub("-duplicate", "", as.character(an.clean$Info.1))
cat.clean$Info.1 <- gsub("-Duplicate", "", as.character(an.clean$Info.1))
cat.clean$Info.1 <- gsub("-DUPLICATE", "", as.character(an.clean$Info.1))
colnames(cat.clean) <- c("Sample.type","Info.1","Lithium Concentration","Sodium Concentration","Ammonium Concentration",
                         "Potassium Concentration","Magnesium Concentration","Calcium Concentration")

#create cation and anion concentration variables
an.conc <- subset(an.clean, select = -Sample.type)
cat.conc <- subset(cat.clean, select = -Sample.type)

#standard deviation workflow:

#find duplicates
an.dups <- an.clean[duplicated(an.clean$Info.1) | duplicated(an.clean$Info.1, fromLast = TRUE),]  #locates any duplicated samples within an.clean
cat.dups <- cat.clean[duplicated(cat.clean$Info.1) | duplicated(cat.clean$Info.1, fromLast = TRUE),]  #locates any duplicated samples within cat.clean

#calculate standard deviation (error) using duplicates
#for anions
an.std <- data.frame(NA, nrow = 1, ncol =6) #create empty data frame
for (i in 1:6){
  an.makeodd <- seq_len(nrow(an.dups)) %% 2                           #creating a sequence for the number of rows of an.dups
  an.diff <- abs(diff(as.numeric(an.dups[,i+2])))[an.makeodd == 1]    #taking the difference between duplicates
  an.diff <- na.omit(an.diff)                                         #omit NA's
  an.std[1,i] <- matrix(mean(an.diff))                                #insert error (mean of the dup diffs) into the matrix
}
an.std <- replace(an.std, an.std == "NaN", NA)
colnames(an.std) <- an.names

#for cations
cat.std <- data.frame(NA, nrow = 1, ncol =6)
for (i in 1:6){
  cat.makeodd <- seq_len(nrow(cat.dups)) %% 2
  cat.diff <- abs(diff(as.numeric(cat.dups[,i+2])))[cat.makeodd == 1]
  cat.diff <- na.omit(cat.diff)
  cat.std[1,i] <- matrix(mean(cat.diff))
}
cat.std <- replace(cat.std, cat.std == "NaN", NA)
colnames(cat.std) <- cat.names

#create std columns
an.std.rep <- matrix(rep(an.std, each = nrow(an.conc)), ncol = 6)        #repeats the an.std values vertically for the length of an.clean
cat.std.rep <- matrix(rep(cat.std, each = nrow(cat.conc)), ncol = 6)        #repeats the cat.std values vertically for the length of cat.clean

#name them
colnames(an.std.rep) <- c("Fluoride STD","Chloride STD","Nitrite STD","Nitrate STD","Phosphate STD","Sulfate STD")
colnames(cat.std.rep) <- c("Lithium STD","Sodium STD","Ammonium STD","Potassium STD","Magnesium STD","Calcium STD")

#detection limit workflow: 

#split the calibration data into concentration and area
#for anions
an.cal.conc <- select(an_calibration, Anions.Fluoride.Concentration, Anions.Chloride.Concentration, Anions.Nitrite.as.N.Concentration,
                      Anions.Nitrate.as.N.Concentration, Anions.Phosphate.as.P.Concentration, Anions.Sulfate.Concentration)
an.cal.area <- select(an_calibration, Anions.Fluoride.Area, Anions.Chloride.Area, Anions.Nitrite.as.N.Area,
                      Anions.Nitrate.as.N.Area, Anions.Phosphate.as.P.Area, Anions.Sulfate.Area)

#for cations
cat.cal.conc <- select(cat_calibration, Cations.Lithium.Concentration, Cations.Sodium.Concentration, Cations.Ammonium.as.N.Concentration,
                       Cations.Potassium.Concentration, Cations.Magnesium.Concentration, Cations.Calcium.Area)
cat.cal.area <- select(cat_calibration, Cations.Lithium.Area, Cations.Sodium.Area, Cations.Ammonium.as.N.Area,
                       Cations.Potassium.Area, Cations.Magnesium.Area, Cations.Calcium.Area)

#create a linear model
#for anions
lm.an <- list()
for (i in 1:6){
  lm.an[[i]] <- lm(an.cal.area[,i] ~ an.cal.conc[,i])     #creates a linear model where area is the response/dependent variable,
}                                                   #and concentration is independent

#for cations
lm.cat <- list()
for (i in 1:6){
  lm.cat[[i]] <- lm(cat.cal.area[,i] ~ cat.cal.conc[,i])
}

#pull residuals and coefficients from the linear model
#for anions
sd.resid.an <- matrix(nrow = 1, ncol = length(lm.an)) #empty matrix
slope.an <- matrix(nrow = 1, ncol = length(lm.an))    #empty matrix
colnames(sd.resid.an) <- an.names
colnames(slope.an) <- an.names

for (t in 1:length(lm.an)){
  sd.resid.an[,t] <- sd(lm.an[[t]][["residuals"]])  #pulls residuals from the linear model
  slope.an[,t] <- lm.an[[t]][[1]][[2]]              #pulls slope coefficients from the linear model
}

#for cations
sd.resid.cat <- matrix(nrow = 1, ncol = length(lm.cat))
slope.cat <- matrix(nrow = 1, ncol = length(lm.cat))
colnames(sd.resid.cat) <- cat.names
colnames(slope.cat) <- cat.names

for (j in 1:6){
  sd.resid.cat[,j] <- sd(lm.cat[[j]][["residuals"]])
  slope.cat[,j] <- lm.cat[[j]][[1]][[2]]
}

#calculate detection limit
#for anions
DL.an <- matrix(nrow = 1, ncol = length(lm.an))  #empty matrix
colnames(DL.an) <- an.names
for (k in 1:length(lm.an)){
  DL.an[,k] <- 3.3 * sd.resid.an[k]/slope.an[k]  #calculating DL using slope and sd of residuals
}

#for cations
DL.cat <- matrix(nrow = 1, ncol = length(lm.cat))
colnames(DL.cat) <- cat.names
for (k in 1:length(lm.cat)){
  DL.cat[,k] <- 3.3 * sd.resid.cat[k]/slope.cat[k]
}

#reprocess data
#for anions
for (y in 1:length(DL.an)){
  an.clean[y+2][an.clean[y+2] < DL.an[y]] <- NA    #inserts NA if any samples are below DL
}

#for cations
for (x in 1:length(DL.cat))
  cat.clean[x+2][cat.clean[x+2] < DL.cat[x]] <- NA    #inserts NA if any samples are below DL

#repeating DL columns for output
#for anions
DL.an.rep <- matrix(rep(DL.an, each = nrow(an.conc)), ncol = 6)
colnames(DL.an.rep) <- c("Fluoride DL","Chloride DL","Nitrite DL","Nitrate DL","Phosphate DL","Sulfate DL")

#for cations
DL.cat.rep <- matrix(rep(DL.cat, each = nrow(cat.conc)), ncol = 6)
colnames(DL.cat.rep) <- c("Lithium DL","Sodium DL","Ammonium DL","Potassium DL","Magnesium DL","Calcium DL")

#full outputs
#anions
anion_processed_output <- cbind(an.conc, an.std.rep, DL.an.rep)
an.sorted.colnames <- sort(names(anion_processed_output))
anion_processed_output <- anion_processed_output[,an.sorted.colnames]
columns.to.move <- c("Info.1")
anion_processed_output <- anion_processed_output[, c(columns.to.move, setdiff(names(anion_processed_output), columns.to.move))]

#cations
cation_processed_output <- cbind(cat.conc, cat.std.rep, DL.cat.rep)
cat.sorted.colnames <- sort(names(cation_processed_output))
cation_processed_output <- cation_processed_output[,cat.sorted.colnames]
columns.to.move <- c("Info.1")
cation_processed_output <- cation_processed_output[, c(columns.to.move, setdiff(names(cation_processed_output), columns.to.move))]


#example export lines:

#write.csv(anion_processed_output,"C:/Users/bryn_/OneDrive/Documents/5108 LAB/Data_analysis_detect_limit/DL_New/Outputs/FILE NAME.csv")

#write.csv(cation_processed_output,"C:/Users/bryn_/OneDrive/Documents/5108 LAB/Data_analysis_detect_limit/DL_New/Outputs/FILE NAME.csv")
