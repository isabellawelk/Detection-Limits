
# Detection limit calculation for output from AA500 
# Modified from original code written by Hannah Richardson

#Detection Limit Equation:
# DL = 3.3 * sigma / m
#  sigma = std deviation of the residuals
#  m = slope of calibration curve
#  residuals = difference in observed results from linear model

###### Input information to set up script and load data: 
 setwd("C:/Users/bryn_/OneDrive/Documents/5108 LAB/Data_analysis_detect_limit") #change this to indicate location of files on your machine
 list.files()
 data <- read.csv("Run data/MDV_NOx_NH3_PO4_9Feb24_Corr.Vals.csv", row.names=NULL) #Change this for the file name you want to calcualte a detection limit for
 
 numchan<- 3 #Will this work? 
 
 ###### Clean Data: 
 library(dplyr)
 
 #extract standards and column names
 standards <- data[data$ANAL %in% c("Sample ID","100", "50", "25", "10", "5", "2", "blank"),]
 
 #set column names that make sense
 names(standards) <- lapply(standards[1, ], as.character)
 standards <- standards[-1,] 
 
 #remove NA columns
 colSums(is.na(standards) | standards=="")
 emptycols<- colSums(is.na(standards) | standards =="")==nrow(standards)
 standards<-standards[,!emptycols]
 
 #Extract concentrations and absorbance
 chan1.conc<-as.numeric(standards$`Results 1`)
 chan1.absorb<-as.numeric(standards$Absorbance)
 chan2.conc<-as.numeric(standards$`Results 2`)
 chan2.absorb<-as.numeric(standards$Absorbance.1)
 chan3.conc<-as.numeric(standards$`Results 3`)
 chan3.absorb<-as.numeric(standards$Absorbance.2)
 
 #bind all together in a matrix
 conc<-t(rbind(chan1.conc, chan2.conc, chan3.conc))
 absorb<- t(rbind(chan1.absorb, chan2.absorb, chan3.absorb))
 
 #fit a regression line to the calibration data, generating a list for each ion of regression statistics
 lm <-list() #allocate list
 for (i in 1:numchan){
   lm[[i]] <- lm(absorb[,i] ~ conc[,i])
 }
 
 #calculate detection limit
 DL<-matrix(nrow=1, ncol=numchan)
 for (i in 1:numchan){
   DL[i]=as.numeric(3.3*sd(lm[[i]][["residuals"]])/lm[[i]][["coefficients"]][2])
 }
 #read the channel from data for proper labeling 
 if (numchan == 2) {
   colnames(DL)<-c(data[8,5],data[8,9])
 } else {
   colnames(DL)<-c(data[8,5],data[8,9],data[8,13])
 } 

 #calculate error
 #find duplicates
 dupsinit<- data[duplicated(data$ANAL) |duplicated(data$ANAL, fromLast=TRUE),]
 #remove non-sample duplicates
 dups<-dupsinit[!(dupsinit$ANAL == "Drift" | dupsinit$ANAL == "Low" | dupsinit$ANAL == "Baseline"),]
 
 #rownames(dups)<-(dups$ANAL)
 
 #read the channel from data for proper labelling 
 if (numchan == 2) {
   dupsclean<-dups[c("ANAL","X.2", "X.6")]
 } else {
   dupsclean<-dups[c("ANAL","X.2", "X.6", "X.10")]
 }
 
 makeodd<-seq_len(nrow(dupsclean)) %% 2
 std<-matrix(nrow=1, ncol=length(DL))
 colnames(std)<-c(data[8,5],data[8,9],data[8,13]) #needs to be fixed for flexible number of channels
 for (i in 1:length(DL)){
   diff<-abs(diff(as.numeric(dupsclean[,i+1])))[makeodd == 1] #difference between rows, just saving the pairs of the same samples. 
   diff<-na.omit(diff) #omit na's 
   std[i] <- sum(diff^2)/(2*(length(dupsclean)/2))
 }
 std 