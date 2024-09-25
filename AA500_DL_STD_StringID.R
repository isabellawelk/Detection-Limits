
# Detection limit calculation for output from AA500 
# Modified from original code written by Hannah Richardson

# This code is a modified version of 'AA500_DL_STD.R' made to work with string variable smaple ID numbers. 
# Use AA500_DL_STD.R if your sample IDs are numeric

#Detection Limit Equation:
# DL = 3.3 * sigma / m
#  sigma = standard deviation of the residuals
#  m = slope of calibration curve
#  residuals = difference in observed results from linear model

rm(list=ls()) #clear previous variables

#### INPUT information to set up script and load data: ####

setwd("/Users/isabellawelk/Documents/GitHub/Detection-Limits/") #copy file path
list.files()
data <- read.csv("NARunData/IW_SiO2_29Aug24.csv", row.names=NULL) #file name for input data (must be .csv format)

numchan<- 1                #Still using for column labeling
NOx<-c()          #INPUT sample number range for NOx, PO4, & NH3
PO4<-c()
NH3<-c()
SiO2<-c("SC-GrF15_C","SC-DF15_A")

samp_range <- c("SC-GrF15_C","SC-DF15_A")  #INPUT FIRST AND LAST SAMPLE NUMBERS

#####load necessary packages

library(dplyr)
library(reshape)

#####find full sample list

#all_samples <- data[which(data$ANAL == samp_range[1]) : which(data$ANAL == samp_range[2]),'ANAL']
#all_samples <- data.frame(as.numeric(all_samples))
#colnames(all_samples) <- ("ANAL")
#all_samples <- na.omit(all_samples)
#all_samples <- all_samples %>% group_by(ANAL) %>% summarise_all(funs(mean)) #Summarise duplicate sample #'s

all_samples <- data[which(data$ANAL == samp_range[1]) : which(data$ANAL == samp_range[2]),'ANAL']
all_samples <- data.frame(all_samples)
colnames(all_samples) <- ("ANAL")
all_samples <- na.omit(all_samples)
all_samples <- all_samples %>% group_by(ANAL) %>% summarise_all(funs(mean)) #Summarise duplicate sample #'s

#####extract standards and column names

standards <- data[data$ANAL %in% c("Sample ID","100", "50", "25", "10", "5", "2", "blank"),]

#####set column names that make sense

names(standards) <- lapply(standards[1, ], as.character)
standards <- standards[-1,] 

#####remove NA columns

colSums(is.na(standards) | standards=="")
emptycols<- colSums(is.na(standards) | standards =="")==nrow(standards)
standards<- standards[,!emptycols]

#####Extract concentrations and absorbance

chan1.conc<-as.numeric(standards$`Results 1`)
chan1.absorb<-as.numeric(standards$Absorbance)
chan2.conc<-as.numeric(standards$`Results 2`)
chan2.absorb<-as.numeric(standards$Absorbance.1)
chan3.conc<-as.numeric(standards$`Results 3`)
chan3.absorb<-as.numeric(standards$Absorbance.2)

#####bind all together in a matrix

conc<-t(rbind(chan1.conc, chan2.conc, chan3.conc))
absorb<- t(rbind(chan1.absorb, chan2.absorb, chan3.absorb))

#####fit a regression line to the calibration data, generating a list for each ion of regression statistics

lm <-list() #allocate list
for (i in 1:numchan){
  lm[[i]] <- lm(absorb[,i] ~ conc[,i])
}

#####calculate detection limit

DL<-matrix(nrow=1, ncol=numchan)
for (i in 1:numchan){
  DL[i]=as.numeric(3.3*sd(lm[[i]][["residuals"]])/lm[[i]][["coefficients"]][2])
}

#####read the channel from data for proper labeling 

if (numchan == 1){
  colnames(DL)<-c(data[8,5])
} else {
  if (numchan == 2){
    colnames(DL)<-c(data[8,5],data[8,9])
  } else {
    colnames(DL)<-c(data[8,5],data[8,9],data[8,13])
  } 
}

#####calculate error

if (length(NOx) > 0)  {
  selected_data_NOx <- data[which(data$ANAL == NOx[1]) : which(data$ANAL == NOx[2]),c(1,which(data[8,] == "NOx"))]
}
if (length(PO4) > 0) {
  selected_data_PO4 <- data[which(data$ANAL == PO4[1]) : which(data$ANAL == PO4[2]),c(1,which(data[8,] == "PO4"))]
}
if (length(NH3) > 0) {
  selected_data_NH3 <- data[which(data$ANAL == NH3[1]) : which(data$ANAL == NH3[2]),c(1,which(data[8,] == "NH3"))]
}
if (length(SiO2) > 0) {
  selected_data_SiO2 <- data[which(data$ANAL == SiO2[1]) : which(data$ANAL == SiO2[2]),c(1,which(data[8,] == "SiO2"))]
}


#####function to clean data, find duplicates, calculate standard deviation, and create output

dups_std <- function(selected_data,col_label,col_select){
  
  #####clean data
  selected_data <- selected_data[!(selected_data$ANAL == "Drift" | selected_data$ANAL == "Baseline" | selected_data$ANAL == "Blank"),]
  selected_data[,1] <- (selected_data[,1]) 
  selected_data[,2] <- as.numeric(selected_data[,2]) 
  
  #####find duplicates
  dups <- selected_data[duplicated(selected_data$ANAL) |duplicated(selected_data$ANAL, fromLast=TRUE),]
  
  #####use duplicates to calculate standard deviation
  makeodd <- seq_len(nrow(dups)) %% 2
  std <- matrix(nrow=1, ncol=length(dups))
  colnames(std) <- c("Sample",toString(col_label))
  
  #####calculate standard deviation   
  diff<-abs(diff(as.numeric(dups[,2])))[makeodd == 1] #difference between rows, just saving the pairs of the same samples. 
  diff<-na.omit(diff) #omit na's 
  
  #####calculate standard deviation
  std <- mean(diff)
  
  #####create data output  
  data_output <- selected_data %>% group_by(ANAL) %>% summarise_all(funs(mean))
  rep_std <- rep(std,nrow(data_output))
  selected_DL <- subset(DL, select = col_select)
  rep_DL <- rep(selected_DL, length(data_output))
  matrix1 <- cbind(data_output,rep_std,selected_DL)
  
}

if (length(NOx) > 0) {
  std_NOx <- dups_std(selected_data_NOx,'Nitrate',"NOx")
} 
if (length(PO4) > 0) {
  std_PO4 <- dups_std(selected_data_PO4,'Phosphate',"PO4")
}
if (length(NH3) > 0) {
  std_NH3 <- dups_std(selected_data_NH3,'Ammonia',"NH3")
}
if (length (SiO2) > 0) {
  std_SiO2 <- dups_std(selected_data_SiO2, 'Silica',"SiO2")
}

full_output <- data.frame(all_samples$ANAL)
colnames(full_output) <- ('ANAL')

#####conditional for number of times to merge data 

if (length(NOx) > 0) {
  full_output <- merge(full_output, std_NOx, by = 'ANAL', all.x = TRUE) 
} 
if (length(PO4) > 0) {
  full_output <- merge(full_output, std_PO4, by = 'ANAL', all.x = TRUE)
}
if (length(NH3) > 0) {
  full_output <- merge(full_output, std_NH3, by = 'ANAL', all.x = TRUE)
}
if (length(SiO2) > 0) {
  #full_output <- merge(full_output, std_SiO2, by = 'ANAL', all.x = TRUE)
  #full_output <- cbind(selected_data_SiO2, std_SiO2) # doesnt work 
}

#colnames(full_output) <- c("Sample","NOx Conc.","NOx STD","NOx DL","PO4 Conc.","PO4 STD","PO4 DL",
#                           "NH3 Conc.","NH3 STD","NH3 DL")
######For de-bugging  
#selected_data <- selected_data_NOx 
#col_label <- "Nitrate"

#colnames(full_output) <- c("Sample", "SiO2 Conc", "SiO2 STD", "SiO2 DL")
#write.csv(full_output,"/Users/isabellawelk/Documents/GitHub/Detection-Limits/Output/FILE NAME.csv")
