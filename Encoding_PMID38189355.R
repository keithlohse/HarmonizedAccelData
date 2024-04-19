#This code script was used to process upper limb accelerometry data for the study, Associations Between Coordination and
#Wearable Sensor Variables (PMID38189355) for the project:
#Harmonized Upper and Lower Limb Accelerometry Data
#Participants were instructed to wear bilateral wrist sensors for two days. In-lab time was removed.

#Clear the workspace and console
rm(list = ls(all = TRUE))
cat("\f")

#Load libraries
library(tidyverse) # suite of ggplot2,tibble,tidyr,purrr,dplyr,stringr,forcats
library(magrittr)  # for extra pipes %$% %<>%
library(labelled) # for working with variable labels (e.g. var_label() & val_label())
library(ggthemes)
library(readxl)
library(graphics)
library(MESS)
library(phonTools)
library(svDialogs) #use guis easily
library(pracma) #for sample_entropy function
library(signal) #for butter and filtfilt functions
library(stats) #for density function
library(lubridate)


#Ask user to the participant ID 
INPUT <-dlgInput(message = "Enter the Participant ID\n(Note: Loading the raw files may be time consuming)", default = NULL,gui = .GUI)  #Answer goes into INPUT$res
Participant_ID <- as.character(INPUT$res)
SUBJECTID <- Participant_ID
print(Participant_ID)

INPUT2 <-  dlgInput(message = "What day? Enter: (Day1 or Day2)", gui = .GUI, default = "Day#") 
DAY <-as.character(INPUT2$res)


#####CALCULATE 1 HZ VARIABLES#####
#read in the 1 Hz files

RAccData <- read.csv(Name.R1, skip = 0, header = T)
RAccData <- as.data.frame(RAccData)
RAccData$TIME <- parse_date_time(RAccData$TIME,orders = c("Ymd HMS","mdy HM"), tz = Sys.timezone())
Master_RAccData <- RAccData #creates a copy

LAccData <- read.csv(Name.L1, skip = 0, header = T)
LAccData <- as.data.frame(LAccData)
LAccData$TIME <- parse_date_time(LAccData$TIME,orders = c("Ymd HMS","mdy HM"), tz = Sys.timezone())
Master_LAccData <- LAccData #creates a copy

  LXData <- LAccData[,5]
  LYData <- LAccData[,6]
  LZData <- LAccData[,7]
  LVMData <- sqrt(LXData^2 + LYData^2 + LZData^2)
  
  RXData <- RAccData[,5]
  RYData <- RAccData[,6]
  RZData <- RAccData[,7]
  RVMData <- sqrt(RXData^2 + RYData^2 + RZData^2)


########## Find variables from the 1Hz full day data  ###############


PLOT.1 <- ggplot() + 
  geom_line( aes(y=LVMData, x= seq(length(LVMData)))) + 
  ylab("Activity Counts") +
  xlab("Time(Sec)") + ggtitle("Full Day") +
  geom_line( aes(y=RVMData, x= seq(length(RVMData))),color="red", alpha = 0.7)

#Find values > Threshold.  
Threshold <- 2 
LCount <- which(LVMData >= Threshold)
RCount <- which(RVMData >= Threshold) 

#Duration variables
recording_time <- round((length(LVMData)/60), 0) #Wearing times in minutes, as integer
l_time <- round(((length(LCount)/60/recording_time)), 2)
r_time <- round(((length(RCount)/60/recording_time)), 2)
NoMvt <- which(LVMData == 0 & RVMData == 0)
total_no_mvt_time <- round(((length(NoMvt)/60/recording_time)), 2)
total_mvt_time <- 1 - total_no_mvt_time
TempData <- which(LVMData > Threshold & RVMData > Threshold)
simultaneous_time <- round((length(TempData)/60/recording_time), 2)
TempData <- which(LVMData > Threshold & RVMData < Threshold) 
l_only_time <- round((length(TempData)/60/recording_time), 2)
TempData <- which(LVMData < Threshold & RVMData > Threshold) 
r_only_time <- round((length(TempData)/60/recording_time), 2)

#Plot Duration variables
ChartData <- data.frame(Category = c("TotalNoMvtTime","LOnlyTime", "ROnlyTime", "SimultaneousTime"),
                        value = c(TotalNoMvtTime, LOnlyTime, ROnlyTime, SimultaneousTime))
b1 <- ggplot(ChartData, mapping = aes(x="", y=value, fill = Category))
b1 <- b1 + geom_bar(width = 1, stat = "identity") + scale_fill_manual(values = c("grey90", "grey70", "grey50", "black"))
b1 <- b1 + labs(x = "Time spent moving", y = "Proportion of Time") + ggtitle("Full Day")
b1


#Magnitude variables, activity counts
l_magnitude <- round(median(LVMData[LCount]), 3) #means when moving
r_magnitude <- round(median(RVMData[RCount]), 3)
l_magnitude_sd <- round(sd(LVMData[LCount]), 3) #SDs when moving
r_magnitude_sd <- round(sd(RVMData[RCount]), 3)
bilateral_magnitude <- l_magnitude + r_magnitude
l_peak_magnitude <- round(max(LVMData[LCount]), 3)
r_peak_magnitude <- round(max(RVMData[RCount]), 3)


#Ratios 
#Since <5, use all L/R
use_ratio <- round(l_time/r_time, 3) 
variation_ratio <- round(l_magnitude_sd/r_magnitude_sd, 3)
simple_magnitude_ratio <- round(l_magnitude/r_magnitude, 3)


#Entropy
#Prior to calculating entropy, need to find the most active hour.
#Fnd length of dataset
LENG <- length(LVMData)
#Subtract the last hour because we will only compare whole hours
LENG.LAST <- LENG - (60 * 60) 
#also set the length of an hour in terms of frames - here for 1 sec sampling rate
HOUR <- (60 * 60)

##Loop
#Set some staring values and vector for the loop below
X=0
MAX = 0
MIN = 120
SAMPLE.AVE <- vector()

#For loop to run through each full hour of the data
#LENG.LAST is the length of the data set up to the last whole hour.
for (X in 1:LENG.LAST){
  #HOUR = the frames in one hour based on the capture hertz
  END <- X + (HOUR-1)
  SAMPLE <- LVMData[X:END]
  MEAN <- mean(SAMPLE)
  SAMPLE.AVE[X] <- MEAN
  
  #If loop to store the maximum average (most active hour)
  if (SAMPLE.AVE[X] > MAX){
    MAX <- SAMPLE.AVE[X]
    START.MAX <- X
  }
}

#Find the END row numbers based on the starting row
END.MAX <- START.MAX + HOUR
RVM.MAX <- RVMData[START.MAX:END.MAX]
LVM.MAX <- LVMData[START.MAX:END.MAX]

#Calculate entropy on just max hour
r_entropy <- round(sample_entropy(RVM.MAX, edim = 2, 0.2*sd(RVM.MAX), tau = 1) ,3)
l_entropy <- round(sample_entropy(LVM.MAX, edim = 2, 0.2*sd(LVM.MAX), tau = 1) ,3)


#####CALCULATE 30 HZ VARIABLES#####


#read in the 30 Hz files

RAccData.30 <- read.csv(Name.R30, skip = 0, header = T)
RAccData.30$VM <- sqrt(RAccData.30$X^2 + RAccData.30$Y^2 + RAccData.30$Z^2)
Master_RAccData.30 <- RAccData.30 #creates a copy 

LAccData.30 <- read.csv(Name.L30, skip = 0, header = T)
LAccData.30$VM <- sqrt(LAccData.30$X^2 + LAccData.30$Y^2 + LAccData.30$Z^2)
Master_LAccData.30 <- LAccData.30 #creates a copy 


################### Full Day 30 Hz Calculations ############

#Standardized Vector Magnitudes

LAccData.30$VM.STD <- as.vector( scale(LAccData.30$VM))
RAccData.30$VM.STD <- as.vector( scale(RAccData.30$VM))


#Filter data
bf <- butter(2, c(I(.2/15),I(12/15)), type="pass")
LAccData.30$VM<- filtfilt(x=LAccData.30$VM, filt = bf)
RAccData.30$VM<- filtfilt(x=RAccData.30$VM, filt = bf)

#Jerk
TP = 1/30 #the sampling interval

#Right jerk
JERK.R.30 <- vector()
JERK.R.30[1]=0

for (i in 1:(nrow(RAccData.30)-1)){
  A=as.numeric(RAccData.30[i+1,8] )#the 8th column is the VM column
  B=as.numeric(RAccData.30[(i),8] )
  JERK.R.30[i] <- (A-B) / TP
  
}

JERK.R.30[nrow(RAccData.30)] = 0
RAccData.30$JERK.R.30 <- JERK.R.30

#Left jerk
JERK.L.30 <- vector()
JERK.L.30[1]=0
i=1
for (i in 1:(nrow(LAccData.30)-1)){
  A=LAccData.30[i+1,8] #the 9th column is the VM column
  B=LAccData.30[i,8]
  JERK.L.30[i] = (A-B) / TP
  
}
JERK.L.30[nrow(LAccData.30)] = 0
LAccData.30$JERK.L.30 <- JERK.L.30

#Jerk output variables
NoJerk.R.30 <- which(RAccData.30$JERK.R.30 == 0)
JERK.R.30 <- RAccData.30$JERK.R.30[-NoJerk.R.30]

NoJerk.L.30 <- which(LAccData.30$JERK.L.30 == 0)
JERK.L.30 <- LAccData.30$JERK.L.30[-NoJerk.L.30]

l_jerk_ave <- round(mean(abs(JERK.L.30)), 3)
r_jerk_ave <- round(mean(abs(JERK.R.30)), 3)
jerk_aym <- (R_Jerk_Avg.30 - L_Jerk_Avg.30) / (L_Jerk_Avg.30 + R_Jerk_Avg.30) #ignoring handedness 

#Frequency variables
#Weighted variance functions from https://stat.ethz.ch/pipermail/r-help/2008-July/168762.html
weighted.var <- function(x, w, na.rm = FALSE) {
  if (na.rm) {
    w <- w[i <- !is.na(x)]
    x <- x[i]
  }
  sum.w <- sum(w)
  sum.w2 <- sum(w^2)
  mean.w <- sum(x * w) / sum(w)
  (sum.w / (sum.w^2 - sum.w2)) * sum(w * (x - mean.w)^2, na.rm =
                                       na.rm)
}

#Find the spectrum, kernal filter = 5
SPEC.R.VM30 <- spectrum(RAccData.30$VM, log='no', span=5, plot = F)
SPEC.L.VM30 <- spectrum(LAccData.30$VM, log='no', span=5, plot = F)

#divide spectrum by sampling interval to create per time instead of per sampling interval
freq.R30 <- SPEC.R.VM30$freq/0.033333
freq.L30 <- SPEC.L.VM30$freq/0.033333

#We should also multiply the spectral density by 2 so that
#the area under the periodogram actually equals the variance of the time series
Density.R30 <- 2*SPEC.R.VM30$spec
Density.L30 <- 2*SPEC.L.VM30$spec
 

#weighted variance of the frequencies present
VAR.WT.R30 <- weighted.var(freq.R30,Density.R30) #
VAR.WT.L30 <- weighted.var(freq.L30,Density.L30) #

r_mean_freq <- round(weighted.mean(freq.R30, Density.R30), 3)
l_mean_freq <- round(weighted.mean(freq.L30, Density.L30), 3)
r_sd_freq <- round(sqrt(VAR.WT.R30), 3) 
l_sd_freq <- round(sqrt(VAR.WT.L30), 3)



######################### SAVING FILES #########################################################


#####SAVE TO FILES for full day #####
WriteCSVFilename <- paste0(SUBJECTID, "_", DAY, "_Full_Output.csv", sep = "") 
WriteRDataFileName <- paste0(SUBJECTID, "_", DAY, "_Full_Output.RData", sep = "") 


#For output to RData file
#Generate dataframe to house variables
Hour.Type <- "fullday"
TaskData <- data.frame(SUBJECTID,Hour.Type,DAY, recording_time, total_no_mvt_time, total_mvt_time, 
                       l_time, l_only_time, r_time, r_only_time, simultaneous_time,
                       l_magnitude, r_magnitude, bilateral_magnitude,
                       l_magnitude_sd, r_magnitude_sd,
                       l_peak_magnitude, r_peak_magnitude,
                       use_ratio, variation_ratio, simple_magnitude_ratio,
                       l_entropy, r_entropy,
                       l_jerk_ave, r_jerk_ave, jerk_aym,
                       l_mean_freq, r_mean_freq, l_sd_freq, r_sd_freq)

##Write CSV file of output

            