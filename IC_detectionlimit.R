IC_detectionlimit <- function (cat,an) {
  cat.clean <- cat[4:10,]
  an.clean <- an[7:13,]
  
  cat.names<-c("Na", "Ammonium", "K", "Mg", "Ca")
  an.names <-c("F", "Cl", "Nitrite", "Nitrate", "Phosphate", "Sulfate")
  
  cat.conc <- cat.clean[,3:7] #separate the x and y data
  cat.area<-cat.clean[,8:12]
  colnames(cat.conc)<- cat.names
  colnames(cat.area)<- cat.names
  
  an.conc <-an.clean[,3:8]
  an.area<-an.clean[,9:14]
  colnames(an.conc)<- an.names
  colnames(an.area)<- an.names
  #cations
  lm.cat <-list() #allocate list
  for (i in 1:5){
    lm.cat[[i]] <- lm(cat.area[ ,i] ~ cat.conc[ ,i])
  }
  
  #anions
  lm.an <-list() #allocate list
  for (i in 1:6){
    lm.an[[i]] <- lm(an.area[ ,i] ~ an.conc[ ,i])
  }
  #cations
  sd.resid.cat <- matrix(nrow = 1, ncol=length(lm.cat))
  slope.cat<-matrix(nrow = 1, ncol=length(lm.cat))
  colnames(sd.resid.cat)<-cat.names
  colnames(slope.cat)<-cat.names
  
  for (n in 1:length(lm.cat)){
    sd.resid.cat[,n]<-sd(lm.cat[[n]][["residuals"]])
    slope.cat[,n]<-lm.cat[[n]][["coefficients"]][["cat.conc[, i]"]]
  }
  
  #anions
  sd.resid.an <- matrix(nrow = 1, ncol=length(lm.an))
  slope.an<-matrix(nrow = 1, ncol=length(lm.an))
  colnames(sd.resid.an)<-an.names
  colnames(slope.an)<-an.names
  
  for (n in 1:length(lm.an)){
    sd.resid.an[,n]<-sd(lm.an[[n]][["residuals"]])
    slope.an[,n]<-lm.an[[n]][["coefficients"]][["an.conc[, i]"]]
  }
  #Cations
  DL.cat<-matrix(nrow = 1, ncol=length(lm.cat))
  colnames(DL.cat)<-cat.names
  for (n in 1:length(lm.cat)){
    DL.cat[,n]=3.3*sd.resid.cat[n]/slope.cat[n]
  }
  
  #Anions
  DL.an<-matrix(nrow = 1, ncol=length(lm.an))
  colnames(DL.an)<-an.names
  for (n in 1:length(lm.an)){
    DL.an[,n]=3.3*sd.resid.an[n]/slope.an[n]
  }
  DL <-list("Cation Detection Limit"=DL.cat, "Anion Detection Limit"=DL.an)
  return(DL)
}


#reprocess data
#cations
for (i in 1:length(DL.cat)){
  cat.samples[i+3][cat.samples[i+3]< DL.cat[i]] <- "<LOD"
}
#anions
for (i in 1:length(DL.an)){
  an.samples[i+3][an.samples[i+3]< DL.an[i]] <- "<LOD"
}

fileName <- paste(wd, 'Processed Anion Data')
write.csv(cat.samples, fileName)
write.csv(an.samples, file = wd) 