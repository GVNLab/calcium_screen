##Load Packages
library(dplyr)
library(ggplot2)
library(pracma)
library(tidyr)
library(stats)

normalize <- function(smoothed.data){
  
  relative_difference <- function(x) 
    (x - min(x))/min(x)
  
  norm.data <- lapply(smoothed.data, relative_difference)
  norm.data$Time <- smoothed.data$Time
  
  norm.data<-as.data.frame(norm.data)
  
  return(norm.data)
  
}

processdata <- function(x){
  normalized.data <- normalize(smoothed.data = x)
  
  
  
  return(normalized.data)
  
}

## Get csv file inputs into workspace.
files <- list.files("./Input/", pattern="*.csv")
sampName <- gsub("*.csv", "", files)
indFiles <- list()
for(k in 1:length(files)){
  indFiles[[sampName[k]]] <- read.table(paste("./Input/",files[k],sep=""), header = TRUE, sep = ",", fileEncoding = "UTF-8-BOM")
}
list2env(indFiles, env = .GlobalEnv) #load into global environment so preprocessing can occur

processedsampNames <- gsub("_t","_processed",sampName) #name of procssessed variables
indProcessed <- list() #make an empty list
for (k in 1:length(processedsampNames)){
  indProcessed[[processedsampNames[k]]] <- processdata(get(sampName[k])) #load each variable named in sampName and run processdata
  
} 

dir.create("./output/")
dir.create("./output/processed/", showWarnings = FALSE)
processedoutputdr <- "./output/processed/"
sapply(names(indProcessed), 
       function (x) write.csv(indProcessed[[x]], file=paste(processedoutputdr,x,".csv"), row.names = TRUE )   )
rm(processedoutputdr)

timetopeak.fun <- function(findpeakinput){ #findpeakinput is a dataframe containing reads from 1 plate values with Time in the first column (output of processdata.function)
  findpeakoutput <- lapply(findpeakinput, pracma::findpeaks, nups = 1, ndowns = 3, minpeakheight = 0.05) #findpeakoutput is a list including Time, of 384 well peaks
  findpeakoutput <- findpeakoutput[-1] #get rid of Time to not analyze.
  n <- length(findpeakoutput)  #number of wells in the input datafile
  time <- data.frame(findpeakinput[1])[,1] #store the time points of the well as a vector - 
  ttpmean <- vector()
  ttpsd <- vector()
  for (k in 1:n){ #for loop through all wells k, 
    
    wellpeaks <- data.frame(findpeakoutput[k]) #peakdata for the well
    if (isempty(wellpeaks)==1){ #if cannot find any peaks, set all values to NA. Else, calculate.
      ttpmean[k] <- NA
      ttpsd[k] <-NA
    }
    else{
      npeaks <- dim(wellpeaks) # number of peaks in well k - returns n x 4 (4 column output of findpeaks)
      npeaks <- npeaks[1] #number of peaks
      ttps = vector() #empty vector of a list of ttps for each well
      
      for (i in 1:npeaks){ #for each peak calculate time to peak
        ttps[i] <- time[wellpeaks[i,2]] - time[wellpeaks[i,3]]
      }
      
      ttpmean[k] <- mean(ttps) #mean ttps for the well
      ttpsd[k] <- sd(ttps)/sqrt(npeaks) #se of the ttps for the well
    }
  }
  
  res_ttp <- data.frame(matrix(nrow=n)) #data frame with 384 wells
  res_ttp$TTP <- ttpmean
  res_ttp$TTP_SE <- ttpsd
  
  name <- colnames(data.frame(findpeakinput))
  name <- gsub(".*\\.","",name) 
  name <- name[-1] #grab row names, remove the time column, and add to dataframe as the name.
  rownames(res_ttp) <- name
  res_ttp <- res_ttp[-1]
  
  
  return(res_ttp)
  
}

amplitude.fun <- function(findpeakinput){
  findpeakoutput <- lapply(findpeakinput, pracma::findpeaks, nups = 1, ndowns = 3, minpeakheight = 0.05) #findpeakoutput is a list including Time, of 384 well peaks
  findpeakoutput <- findpeakoutput[-1] #get rid of Time to not analyze.
  n <- length(findpeakoutput)  #number of wells in the input datafile
  time <- data.frame(findpeakinput[1])[,1] #store the time points of the well as a vector - 
  ampmean <- vector()
  ampsd <- vector()
  bpmmean <- vector()
  bpmsd <- vector()
  number <- vector()
  b2bmean <- vector()
  b2bsd <- vector()
  for (k in 1:n){ #for loop through all wells k, 
    wellpeaks <- data.frame(findpeakoutput[k]) #peakdata for the well
    if (isempty(wellpeaks)==1){ #if cannot find any peaks, set all values to NA. Else, calculate.
      ampmean[k] <- NA
      ampsd[k] <- NA
      bpmmean[k] <- NA
      bpmsd[k] <- NA
      number[k] <- NA
      b2bmean[k] <- NA
      b2bsd[k] <- NA
    }
    else{
      npeaks <- dim(wellpeaks) # number of peaks in well k - returns n x 4 (4 column output of findpeaks)
      npeaks <- npeaks[1] #number of peaks
      
      
      
      
      b2b <- vector() #empty vector for beat-beat time.
      
      for (i in 2:npeaks){ #loop through each peak to find beat to beat time.
        if(isempty(wellpeaks)==1){
          bpm[i-1] <- NA
        }
        else{
          b2b[i-1] <- time[wellpeaks[i,2]]-  time[wellpeaks[i-1,2]]
        }
        
      }
      bpm<-60/b2b
      ampmean[k] <- mean(wellpeaks[,1])
      ampsd[k] <- sd(wellpeaks[,1])/sqrt(npeaks)
      bpmmean[k] <- 60/mean(b2b)
      #bpmsd[k] <- 60/sd(b2b)
      b2bmean[k] <- mean(b2b)
      b2bsd[k] <- sd(b2b)/sqrt(npeaks)
      number[k] <- npeaks[1]
    }
  }
  res_amp <- data.frame(matrix(nrow=n)) #data frame with 384 wells
  #res_amp$BPM_SD <- bpmsd
  res_amp$Amp <- ampmean
  res_amp$Amp_SE <- ampsd
  res_amp$Beat2Beat <- b2bmean
  res_amp$Beat2Beat_SE <- b2bsd
  res_amp$BPM <- bpmmean
  res_amp$PeaksID <- number
  
  name <- colnames(data.frame(findpeakinput))
  name <- gsub(".*\\.","",name) 
  name <- name[-1] #grab row names, remove the time column, and add to dataframe as the name.
  rownames(res_amp) <- name
  res_amp <- res_amp[-1]
  
  
  return(res_amp)
}

tau.fun <- function(normalizeddf) {
  tau <- matrix(nrow =2, ncol = length(normalizeddf)-1)
  time <- normalizeddf[,1]
  for (i in 2:length(normalizeddf)){ #iterates over all wells
    y <- normalizeddf[,i]
    peaks <- pracma::findpeaks(normalizeddf[,i],nups = 1, ndowns = 3, minpeakheight = 0.05) # find peaks of the well
    if (isempty(peaks)==1){ #if cannot find any peaks, set all values to NA. Else, calculate.
      tau[1,i-1] <- NA
      tau[2,i-2] <- NA
    }
    else{
      dim <- size(peaks)
      npeaks <- dim[1]
      tautemp = vector()
      for (idx in 1:dim[1]){
        #store variables for easy access
        
        peakstart <- peaks[idx,3]
        peakloc <- peaks[idx,2]
        peakend <- peaks[idx,4] - 1 #omit last value so that no 0's are in the curve for regression.
        ylocal <- y[peakloc:peakend] #local curve
        timemax <- time[peakloc]
        timelocal <- time[peakloc:peakend]
        results<-lm(formula = log(na.omit(ylocal)) ~ na.omit(timelocal)) #use exponetial regression to find a model
        
        #plot to see how the curve matches
        # xpredict <- data.frame(seq(from = timelocal[1], to = timelocal[length(timelocal)], length.out = 100))
        # ypredict <- exp(predict.lm(results, newdata = xpredict))
        #plot(ylocal~timelocal)
        # plot(ypredict~xpredict)
        #break
        
        tautemp[idx] <- results$coefficients[2]*-1  #extract lambda
        
        
        
      }
      tautemp <- 1/tautemp #convert lambda into tau
      tau[1,i-1] <- mean(tautemp) #mean tau of all peaks in well
      tau[2,i-1] <- sd(tautemp)/sqrt(npeaks)
    }
  }
  c <- colnames(normalizeddf[-1])
  colnames(tau)<-c
  rownames(tau) <- c("Tau", "Tau_SE")
  tau <- t(tau)
  name <- colnames(data.frame(normalizeddf))
  name <- gsub(".*\\.","",name) 
  name <- name[-1] #grab row names, remove the time column, and add to dataframe as the name.
  rownames(tau) <- name #labeled wells
  
  return(tau)
}

master_tau <- list()
for (k in 1:length(indProcessed)){ 
  result<- as.data.frame(tau.fun(data.frame(indProcessed[k])))
  master_tau[k] <- list(result)
} #master list of tau results

master_amp <- list()
for (k in 1:length(indProcessed)){ 
  result<- as.data.frame(amplitude.fun(data.frame(indProcessed[k])))
  master_amp[k] <- list(result)
} #master list of amp results

master_ttp <- list()
for (k in 1:length(indProcessed)){ 
  result<- as.data.frame(timetopeak.fun(data.frame(indProcessed[k])))
  master_ttp[k] <- list(result)
} #master list of ttp results

master_results <- list()
for (k in 1:length(indProcessed)){
  result <- as.data.frame(bind_cols(master_tau[[k]],master_amp[[k]],master_ttp[[k]]))
  master_results[k] <- list(result)
}

resultNames <- gsub("_processed","_results",processedsampNames) #name of procssessed variables
names(master_results) <- resultNames
dir.create("./output/results/", showWarnings = FALSE)
resultoutputdir <- "./output/results/"
sapply(names(master_results), 
       function (x) write.csv(master_results[[x]], file=paste(resultoutputdir,x,".csv"), row.names = TRUE )   ) #write to csv in the working directory
