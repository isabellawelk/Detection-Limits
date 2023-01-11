AA500_detectionlimit <- function (data,numchan) {
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

#create linear model
lm <-list() #allocate list
for (i in 1:numchan){
  lm[[i]] <- lm(absorb[,i] ~ conc[,i])
}

#Calculate detection limit
DL<-matrix(nrow=1, ncol=numchan)
for (i in 1:numchan){
  DL[i]=as.numeric(3.3*sd(lm[[i]][["residuals"]])/lm[[i]][["coefficients"]][2])
}
#read names from spreadsheet to determine channels used, add to output
colnames(DL)<-c(data[8,5],data[8,9],data[8,13])
return(DL) 
}

