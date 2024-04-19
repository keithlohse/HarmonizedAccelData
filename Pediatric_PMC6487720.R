#This code script was used to process upper limb accelerometry data for the study, Detection of Pediatric Upper Extremity 
#Motor Activity and Deficits with Accelerometry (PMC6487720) for the project:
#Harmonized Upper and Lower Limb Accelerometry Data
#Accelerometry data from this cohort was processed one day at a time

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

#Ask user to set directory, then execute it so the directory actually changes 
TempPath <- dlg_dir(title = "Set directory that contains the subject's files")  
PATH <-TempPath$res
setwd(PATH)

#set theme to minimal
theme_set(theme_minimal()) #all plots that run below use this theme

#####CALCULATE 1 HZ VARIABLES#####
#Read in 1 Hz csv files generated from the Data Extraction script
dlg_message("Please load a Left, 1 Hz csv file.", gui = .GUI)
L.FILE <- dlg_open(PATH, title = "Select the Left, 1 Hz csv file", filters = dlg_filters[c("All"),])$res
LAccData <- read.csv(L.FILE, skip = 10, header = T)
LAccData <- as.data.frame(LAccData)

#Pull off file name, then find and parse string to create parts of output file name
TempFileName <- strsplit(L.FILE, "/")
TempFileName <- matrix(unlist(TempFileName), ncol=1, byrow=TRUE)
TempLength <- length(TempFileName)
stringtemp <- TempFileName[TempLength] #leaves you with filename
stringtemp <- strsplit(stringtemp, "_")
stringtemp <- matrix(unlist(stringtemp), ncol=1, byrow=TRUE)
SUBJECTID <- as.character(stringtemp[1]) #subject ID
DAY <- as.character(stringtemp[2]) # day
#Change v1 - V4 to day1 - day4
if (DAY == "V1"){
  DAY <- "Day1" }
if (DAY == "V2"){
  DAY <- "Day2" }
if (DAY == "V3"){
  DAY <- "Day3" }
if (DAY == "V4"){
  DAY <- "Day4" }

#Load right side automatically
R.FILE <- paste0(stringtemp[1], "_", stringtemp[2], "_", stringtemp[3], "_R_", stringtemp[5], "_1sec.csv") 
RAccData <- read.csv(R.FILE, skip = 10, header = T)
RAccData <- as.data.frame(RAccData)


#Assign assessment time point (event name) for data harmonization
redcap_event_name <- "baseline"

#Trim off 30 min at beginning and end of files
LAccData <- LAccData[1800:88199, ]
RAccData <- RAccData[1800:88199, ]

#Computing the vector magnitude of the left side data
LXData <- LAccData[,1]
LYData <- LAccData[,2]
LZData <- LAccData[,3]
LVMData <- sqrt(LXData^2 + LYData^2 + LZData^2)

#Computing the vector magnitude of the right side data
RXData <- RAccData[,1]
RYData <- RAccData[,2]
RZData <- RAccData[,3]
RVMData <- sqrt(RXData^2 + RYData^2 + RZData^2)

#Plot to check vector magnitudes, making sure two sides ~match
plot(LVMData, type="l", col="black", 
     ylab="Activity Counts", xlab = "Time(sec)")
matplot(RVMData, type="l", col="red", add = 1)

#Find values > Threshold.  
Threshold <- 2 
LCount <- which(LVMData >= Threshold)
RCount <- which(RVMData >= Threshold) 

#Duration variables
#Most time variables as %, to account for differing durations of task. 
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
b1 <- b1 + labs(x = "Time spent moving", y = "Proportion of Time")
b1


#Calculating magnitude variables from the 1Hz data, in activity counts
l_magnitude <- round(median(LVMData[LCount]), 3) #means when moving
r_magnitude <- round(median(RVMData[RCount]), 3)
l_magnitude_sd <- round(sd(LVMData[LCount]), 3) #SDs when moving
r_magnitude_sd <- round(sd(RVMData[RCount]), 3)
bilateral_magnitude <- l_magnitude + r_magnitude
l_peak_magnitude <- round(max(LVMData[LCount]), 3)
r_peak_magnitude <- round(max(RVMData[RCount]), 3)

#Ratios where affected or non-dominant side comes from file name
if (stringtemp[5] == "LH") {
  AffSide <- "R"
} else {AffSide <- "L"}


if (AffSide == "L") {
  use_ratio <- round(l_time/r_time, 3)
  variation_ratio <- round(l_magnitude_sd/r_magnitude_sd, 3)
  simple_magnitude_ratio <- round(l_magnitude/r_magnitude, 3)
} else {
  use_ratio <- round(r_time/l_time, 3)
  variation_ratio <- round(r_magnitude_sd/l_magnitude_sd, 3)
  simple_magnitude_ratio <- round(r_magnitude/l_magnitude, 3)
}


#Calculating entropy
#Prior to calculating entropy, need to find the most active hour on most active side.
if (AffSide == "R") {
  #Use left side to find length of data
  LENG <- length(LVMData)
} else {
  #Use right side to find length of data
  LENG <- length(RVMData)
}

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
  if (AffSide == "R") {
    SAMPLE <- LVMData[X:END]
  } else {
    SAMPLE <- RVMData[X:END]
  }

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

#Calculate entropy on just the hour of maximum activity
r_entropy <- round(sample_entropy(RVM.MAX, edim = 2, 0.2*sd(RVM.MAX), tau = 1) ,3)
l_entropy <- round(sample_entropy(LVM.MAX, edim = 2, 0.2*sd(LVM.MAX), tau = 1) ,3)


#####CALCULATE 30 HZ VARIABLES#####
##Load in associated 30Hz files###

#Change to 30Hz file folder
#Rebuild path name and switch to HarmonizedOutput directory
RawPath <- file.path(TempFileName[1], TempFileName[2], TempFileName[3],
                        TempFileName[4], TempFileName[5], TempFileName[6],
                        TempFileName[7], TempFileName[8], "BPEAC_30HzData")
setwd(RawPath)

#Read in 30 Hz csv files generated from the Data Extraction script
L.FILE <- paste0(stringtemp[1], "_", stringtemp[2], "_", stringtemp[3], "_L_", stringtemp[5], "_RAW.csv")
LAccData <- read.csv(L.FILE, skip = 11, header = T)

R.FILE <- paste0(stringtemp[1], "_", stringtemp[2], "_", stringtemp[3], "_R_", stringtemp[5], "_RAW.csv")
RAccData <- read.csv(R.FILE, skip = 11, header = T)

#Trim 30 Hz data to 24 hrs
#START_ROW <- START_ROW *30
LAccData <- LAccData[54000:2645970, ]
RAccData <- RAccData[54000:2645970, ]
DATA.1.L <- LAccData
DATA.1.R <- RAccData

#rename columns
colnames(DATA.1.L) <- c("X", "Y", "Z")
colnames(DATA.1.R) <- c("X", "Y", "Z")

#Calculate vector magnitudes
DATA.1.L$VM <-sqrt(DATA.1.L$X^2 + DATA.1.L$Y^2 + DATA.1.L$Z^2 )
DATA.1.L$VM.STD <- as.vector( scale(DATA.1.L$VM))

DATA.1.R$VM <-sqrt(DATA.1.R$X^2 + DATA.1.R$Y^2 + DATA.1.R$Z^2 )
DATA.1.R$VM.STD <- as.vector( scale(DATA.1.R$VM))

#Filter data
bf <- butter(2, c(I(.2/15),I(12/15)), type="pass")
DATA.1.L$VM<- filtfilt(x=DATA.1.L$VM, filt = bf)
DATA.1.R$VM<- filtfilt(x=DATA.1.R$VM, filt = bf)
ycol <- which(colnames(DATA.1.L) == "VM") #finds column number for vector magnitude

#Calculating jerk- a measure of smoothness of movement
TP = 1/30 #the sampling interval

#Right jerk
JERK.R <- vector()
JERK.R[1]=0

for (i in 1:(nrow(DATA.1.R)-1)){
  A=as.numeric(DATA.1.R[i+1,ycol] )
  B=as.numeric(DATA.1.R[(i),ycol] )
  JERK.R[i] <- (A-B) / TP
  
}

JERK.R[nrow(DATA.1.R)] = 0
DATA.1.R$JERK.R <- JERK.R

#Left jerk
JERK.L <- vector()
JERK.L[1]=0
i=1
for (i in 1:(nrow(DATA.1.L)-1)){
  A=DATA.1.L[i+1,ycol] 
  B=DATA.1.L[i,ycol]
  JERK.L[i] = (A-B) / TP
  
}
JERK.L[nrow(DATA.1.L)] = 0
DATA.1.L$JERK.L <- JERK.L

#Jerk output variables
NoJerk.R <- which(DATA.1.R$JERK.R == 0)
JERK.R <- DATA.1.R$JERK.R[-NoJerk.R]

NoJerk.L <- which(DATA.1.L$JERK.L == 0)
JERK.L <- DATA.1.L$JERK.L[-NoJerk.L]

l_jerk_ave <- round(mean(abs(JERK.L)), 3)
r_jerk_ave <- round(mean(abs(JERK.R)), 3)

if (AffSide == "L") {
  jerk_aym <- (l_jerk_ave - r_jerk_ave) / (l_jerk_ave + r_jerk_ave)
} else {
  jerk_aym <- (r_jerk_ave - l_jerk_ave) / (l_jerk_ave + r_jerk_ave)
}


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
SPEC.R.VM <- spectrum(DATA.1.R$VM, log='no', span=5, plot = F)
SPEC.L.VM <- spectrum(DATA.1.L$VM, log='no', span=5, plot = F)

#divide spectrum by sampling interval to create per time instead of per sampling interval
freq.R <- SPEC.R.VM$freq/0.033333
freq.L <- SPEC.L.VM$freq/0.033333

#Multiplying the spectral density by 2 so that the area under the periodogram equals the variance of the time series
Density.R <- 2*SPEC.R.VM$spec
Density.L <- 2*SPEC.L.VM$spec

#weighted variance of the frequencies present
VAR.WT.R <- weighted.var(freq.R,Density.R) #
VAR.WT.L <- weighted.var(freq.L,Density.L) #

r_mean_freq <- round(weighted.mean(freq.R, Density.R), 3)
l_mean_freq <- round(weighted.mean(freq.L, Density.L), 3)
r_sd_freq <- round(sqrt(VAR.WT.R), 3) 
l_sd_freq <- round(sqrt(VAR.WT.L), 3)

#####SAVE TO FILES#####
WriteCSVFilename <- paste0(SUBJECTID, "_", DAY, "_HarmOutput.csv", sep = "") 

#For output to RData file
#Generate dataframe to house variables
TaskData <- data.frame(SUBJECTID, redcap_event_name, recording_time, total_no_mvt_time, total_mvt_time, 
                          l_time, l_only_time, r_time, r_only_time, simultaneous_time,
                          l_magnitude, r_magnitude, bilateral_magnitude,
                          l_magnitude_sd, r_magnitude_sd,
                          l_peak_magnitude, r_peak_magnitude,
                          use_ratio, variation_ratio, simple_magnitude_ratio,
                          l_entropy, r_entropy,
                          l_jerk_ave, r_jerk_ave, jerk_aym,
                          l_mean_freq, r_mean_freq, l_sd_freq, r_sd_freq)

#Switch to Output directory
#Rebuild path name and switch to HarmonizedData/Output directory
OutputPath <- file.path(TempFileName[1], TempFileName[2], TempFileName[3],
                        TempFileName[4], TempFileName[5], TempFileName[6],
                        TempFileName[7], "HarmonizedData", "Output")
setwd(OutputPath)

#For output to csv file
if (DAY == "Day1") {
  names(TaskData) <- c("participant_id", "redcap_event_name", "d1_wear_time", "d1_total_no_mvt_time", "d1_total_mvt_time", 
                     "d1_l_time", "d1_l_only_time", "d1_r_time", "d1_r_only_time", "d1_simultaneous_time",
                     "d1_l_magnitude", "d1_r_magnitude", "d1_bilateral_magnitude",
                     "d1_l_magnitude_sd", "d1_r_magnitude_sd",
                     "d1_l_peak_magnitude", "d1_r_peak_magnitude",
                     "d1_use_ratio", "d1_variation_ratio", "d1_simple_magnitude_ratio",
                     "d1_vm_corrcoef", "d1_l_entropy", "d1_r_entropy",
                     "d1_l_jerk_avg", "d1_r_jerk_avg", "d1_jerk_aym",
                     "d1_l_mean_freq", "d1_r_mean_freq", "d1_l_sd_freq", "d1_r_sd_freq",
                     "d1_l_arc_length", "d1_r_arc_length")
} 
if (DAY == "Day2") {
  names(TaskData) <- c("participant_id", "redcap_event_name", "d2_wear_time", "d2_total_no_mvt_time", "d2_total_mvt_time", 
                       "d2_l_time", "d2_l_only_time", "d2_r_time", "d2_r_only_time", "d2_simultaneous_time",
                       "d2_l_magnitude", "d2_r_magnitude", "d2_bilateral_magnitude",
                       "d2_l_magnitude_sd", "d2_r_magnitude_sd",
                       "d2_l_peak_magnitude", "d2_r_peak_magnitude",
                       "d2_use_ratio", "d2_variation_ratio", "d2_simple_magnitude_ratio",
                       "d2_vm_corrcoef", "d2_l_entropy", "d2_r_entropy",
                       "d2_l_jerk_avg", "d2_r_jerk_avg", "d2_jerk_aym",
                       "d2_l_mean_freq", "d2_r_mean_freq", "d2_l_sd_freq", "d2_r_sd_freq",
                       "d2_l_arc_length", "d2_r_arc_length") 
}
if (DAY == "Day3") {
  names(TaskData) <- c("participant_id", "redcap_event_name", "d3_wear_time", "d3_total_no_mvt_time", "d3_total_mvt_time", 
                       "d3_l_time", "d3_l_only_time", "d3_r_time", "d3_r_only_time", "d3_simultaneous_time",
                       "d3_l_magnitude", "d3_r_magnitude", "d3_bilateral_magnitude",
                       "d3_l_magnitude_sd", "d3_r_magnitude_sd",
                       "d3_l_peak_magnitude", "d3_r_peak_magnitude",
                       "d3_use_ratio", "d3_variation_ratio", "d3_simple_magnitude_ratio",
                       "d3_vm_corrcoef", "d3_l_entropy", "d3_r_entropy",
                       "d3_l_jerk_avg", "d3_r_jerk_avg", "d3_jerk_aym",
                       "d3_l_mean_freq", "d3_r_mean_freq", "d3_l_sd_freq", "d3_r_sd_freq",
                       "d3_l_arc_length", "d3_r_arc_length") 
}
if (DAY == "Day4") {
  names(TaskData) <- c("participant_id", "redcap_event_name", "d4_wear_time", "d4_total_no_mvt_time", "d4_total_mvt_time", 
                       "d4_l_time", "d4_l_only_time", "d4_r_time", "d4_r_only_time", "d4_simultaneous_time",
                       "d4_l_magnitude", "d4_r_magnitude", "d4_bilateral_magnitude",
                       "d4_l_magnitude_sd", "d4_r_magnitude_sd",
                       "d4_l_peak_magnitude", "d4_r_peak_magnitude",
                       "d4_use_ratio", "d4_variation_ratio", "d4_simple_magnitude_ratio",
                       "d4_vm_corrcoef", "d4_l_entropy", "d4_r_entropy",
                       "d4_l_jerk_avg", "d4_r_jerk_avg", "d4_jerk_aym",
                       "d4_l_mean_freq", "d4_r_mean_freq", "d4_l_sd_freq", "d4_r_sd_freq",
                       "d4_l_arc_length", "d4_r_arc_length") 
  }

#Save data to csv file
write.csv(TaskData, file = WriteCSVFilename, row.names = F)

