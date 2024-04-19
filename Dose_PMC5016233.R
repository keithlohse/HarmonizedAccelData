#This code script was used to process upper limb accelerometry data from the study, Does Task-Specific Training Improve Upper
#Limb Performance in Daily Life Poststroke? (PMC5016233) for the project:
#Harmonized Upper and Lower Limb Accelerometry Data
#The code reads in a data processing log that contains the participant's ID number (column 1),
#time point of data collection (column 2), and affected side (column 5). Note that this script processes
#one day of accelerometry data, per the study protocol

library(Rmisc)
library(tidyverse)
library(ggplot2)
library(psych) 
library(lsr)
library(ggrepel)
library(car)
library(irr)
library(svDialogs)
library(forecast)
library(tcltk)
library(shiny)
library(GGally)
library(gWidgets2)
library(signal) #for butter and filtfilt functions
library(pracma) #for sample_entropy function
library(stats) #for density function
library(stringi)


##READ IN DATA PROCESSING FILE
dlg_message("Please load the data processing file.", gui = .GUI)
FILE <- dlg_open(title = "Select *.csv data processing file", filters = dlg_filters[c("All"),])$res
DATAPROCESSINGLOG <- read.csv(FILE, strip.white = TRUE)
colnames(DATAPROCESSINGLOG)[1] <- "StudyID"
colnames(DATAPROCESSINGLOG)[2] <- "TimePoint"
colnames(DATAPROCESSINGLOG)[5] <- "AffectedSide"

#Convert AffectedSide to Left and Right:
DATAPROCESSINGLOG$AffectedSide <- case_match(DATAPROCESSINGLOG$AffectedSide, 1 ~ "Left", 
                                                                             2 ~ "Right")
#Ask user to set directory to where the raw data files are located
TempPath <- dlg_dir(Message = "Set directory for opening files")  
PATH <-TempPath$res
main_dir <- setwd(PATH)

#List of subject folders
SUBJECT_LIST <-list.files(main_dir)

#Prune the non-folders out of the list:
dir.exists(SUBJECT_LIST)
SUBJECT_LIST <- SUBJECT_LIST[dir.exists(SUBJECT_LIST)==TRUE]

#Check the sums to make sure the appropriate number of subjects have been selected
sum(dir.exists(SUBJECT_LIST))

#The for-loop below will loop through each row of the data processing log, identify the subject directory
#and time point to process, read in the 1 Hz data files from the left and right sides, process 
#upper limb sensor variables from the 1 Hz data, read in the 30 Hz data files from the left and right sides,
#process upper limb sensor variables from the 30 Hz data, save all computed variables to a data frame, 
#and write a csv file of the computed variables

for(s in 1:nrow(DATAPROCESSINGLOG)) {
    sub_id <- DATAPROCESSINGLOG$StudyID[s]
    sub_dir <- paste(main_dir, "/",sub_id,"/", sep="")
    list.files(sub_dir)
    print(sub_id)
    
    #Match time point on data processing log to time point location on Box  
    ProcessedTimePoint <- DATAPROCESSINGLOG$TimePoint[s]
    time_folder <- list.files(sub_dir, pattern = ProcessedTimePoint)
    time_dir <- paste(sub_dir, time_folder, "/", sep = "")
    print(ProcessedTimePoint)

    #Extract affected side
    AffSide <- DATAPROCESSINGLOG$AffectedSide[s]
    
    # Inside of each time point
    L_dir <- list.files(time_dir, pattern = "LUE")
    L_dir <- paste(time_dir, L_dir, "/", sep = "")

    R_dir <- list.files(time_dir, pattern = "RUE")
    R_dir <- paste(time_dir, R_dir, "/", sep = "")

    # Inside of each LUE/RUE folder
    L_1Hz <- list.files(L_dir, pattern = "1sec.csv")
    L_1Hz <- paste(L_dir, L_1Hz, sep = "")
    
    ##Pull off start date and time from L 1Hz file (assumes sensors were programmed to start recording
    #at the same date and time for both left and right sensors)
    L_1Hz_AccelData <- read.csv(L_1Hz, header = F)
    LStartTime <- L_1Hz_AccelData[3,1]
    LStartTime <- strsplit(LStartTime, " ")
    LStartTime <- LStartTime[[1]][3]
    LStartTime <- str_split(LStartTime, ":")
    LStartHour <- LStartTime[[1]][1]
    LStartMinute <- LStartTime[[1]][2]
    LStartSecond <- LStartTime[[1]][3]

    LStartDate <- L_1Hz_AccelData[4,1]
    LStartDate <- strsplit(LStartDate, " ")
    LStartDate <- LStartDate[[1]][3]
    LStartDate <- as.Date(LStartDate, format = "%m/%d/%Y")
    LStartYear <- substring(LStartDate, 1, 4)
    LStartMonth <- substring(LStartDate, 6, 7)
    LStartDay <- substring(LStartDate, 9, 10)
    
    L_1Hz_AccelData <- read.csv(L_1Hz, skip = 10, header = F)
    L_1Hz_AccelData <- L_1Hz_AccelData[ ,1:3]
    names(L_1Hz_AccelData) <- c("X.L", "Y.L", "Z.L")

    R_1Hz <- list.files(R_dir, pattern = "1sec.csv")
    R_1Hz <- paste(R_dir, R_1Hz, sep = "")
    
    R_1Hz_AccelData <- read.csv(R_1Hz, skip = 10, header = F)
    R_1Hz_AccelData <- R_1Hz_AccelData[ ,1:3]
    names(R_1Hz_AccelData) <- c("X.R", "Y.R", "Z.R")

    ##If left and right sensors have different recording lengths, set LENG_1 to be the shorter of the two:
    if(nrow(L_1Hz_AccelData) == nrow(R_1Hz_AccelData)){
              LENG_1 <- nrow(L_1Hz_AccelData)
    }
    
    if(nrow(L_1Hz_AccelData) < nrow(R_1Hz_AccelData)){
              LENG_1 <- nrow(L_1Hz_AccelData)
              R_1Hz_AccelData <- R_1Hz_AccelData[1:LENG_1, ]
    }
    
    if(nrow(R_1Hz_AccelData) < nrow(L_1Hz_AccelData)){
              LENG_1 <- nrow(R_1Hz_AccelData)
              L_1Hz_AccelData <- L_1Hz_AccelData[1:LENG_1, ]
    }
    
    #Combining left and right 1Hz files into one data frame
    Combined1Hz <- cbind(L_1Hz_AccelData, R_1Hz_AccelData)
    
    #Establishing the sampling frequency
    HERTZ_1 <- 1
    
    #Adding a time stamp (TIME_STAMP) to the data
    TIME_STAMP <- seq(ISOdate(LStartYear,LStartMonth,LStartDay,LStartHour,LStartMinute,LStartSecond), by = (1/HERTZ_1), length.out = LENG_1)
    TIME_STAMP <- as.character(TIME_STAMP)
    ORIG_COUNT <- seq(LENG_1)
    Combined1Hz = data.frame(ORIG_COUNT,TIME = TIME_STAMP, Combined1Hz)
    
    #Remove the first 1.5 hours of the recording due to in-lab time
    InLabTime_1Hz <- 1.5 * 3600
    Combined1Hz <- Combined1Hz[-c(1:InLabTime_1Hz), ]
    
    #Cut file if > 24 hours
    if(nrow(Combined1Hz) > 86400){
      Combined1Hz <- Combined1Hz[1:86400, ]
    }

    #####CALCULATE 1 HZ VARIABLES#####
    #Computing the vector magnitude of the left side data
    LXData <- Combined1Hz[,"X.L"]
    LYData <- Combined1Hz[,"Y.L"]
    LZData <- Combined1Hz[,"Z.L"]
    Combined1Hz$LVMData <- sqrt(LXData^2 + LYData^2 + LZData^2)

    #Computing the vector magnitude of the right side data
    RXData <- Combined1Hz[,"X.R"]
    RYData <- Combined1Hz[,"Y.R"]
    RZData <- Combined1Hz[,"Z.R"]
    Combined1Hz$RVMData <- sqrt(RXData^2 + RYData^2 + RZData^2)
    
    #Index values > Threshold.  
    Threshold <- 2 
    
    Combined1Hz$LCount <- ifelse(Combined1Hz$LVMData >= Threshold, 1, 0)
    Combined1Hz$RCount <- ifelse(Combined1Hz$RVMData >= Threshold, 1, 0)
    
    #Calculating duration variables from the 1Hz data
    Hz1DurationMetrics <- Combined1Hz %>% dplyr::summarise(recording_time = round((length(LVMData)/60), 0),  #Wearing times in minutes, as integer
                                                           NoMvt = sum(LVMData == 0 & RVMData == 0)/60,
                                                           total_no_mvt_time = round(NoMvt/recording_time, 2),
                                                           total_mvt_time = 1 - total_no_mvt_time, 
                                                           LCount = sum((LCount == 1)/60),
                                                           l_time = round(LCount/recording_time, 2),
                                                           RCount = sum((RCount == 1)/60),
                                                           r_time = round(RCount/recording_time, 2),
                                                           SimTime = sum(LVMData > Threshold & RVMData > Threshold),
                                                           simultaneous_time = round((SimTime/60/recording_time), 2),
                                                           LVMAboveThreshold = sum(LVMData > Threshold & RVMData < Threshold),
                                                           l_only_time = round(LVMAboveThreshold/60/recording_time, 2),
                                                           RVMAboveThreshold = sum(LVMData < Threshold & RVMData > Threshold),
                                                           r_only_time = round(RVMAboveThreshold/60/recording_time, 2))
    
    #Calculating magnitude variables from the 1Hz data, in activity counts
    l_magnitude <- Combined1Hz %>% dplyr::filter(., LCount == "1") %>% dplyr::summarise(l_magnitude = round(median(LVMData), 3))
    r_magnitude <- Combined1Hz %>% dplyr::filter(., RCount == "1") %>% dplyr::summarise(r_magnitude = round(median(RVMData), 3))
    l_magnitude_sd <- Combined1Hz %>% dplyr::filter(., LCount == "1") %>% dplyr::summarise(l_magnitude_sd = round(sd(LVMData), 3))
    r_magnitude_sd <- Combined1Hz %>% dplyr::filter(., RCount == "1") %>% dplyr::summarise(r_magnitude_sd = round(sd(RVMData), 3))
    bilateral_magnitude <- sum(l_magnitude + r_magnitude)
    l_peak_magnitude <- Combined1Hz %>% dplyr::filter(., LCount == "1") %>% dplyr::summarise(l_peak_magnitude = round(max(LVMData), 3))
    r_peak_magnitude <- Combined1Hz %>% dplyr::filter(., RCount == "1") %>% dplyr::summarise(r_peak_magnitude = round(max(RVMData), 3))
    
    
    #Calculating ratios variables from the 1Hz data
    if (AffSide == "Left") {
      use_ratio <- Hz1DurationMetrics %>% dplyr::summarise(use_ratio = round(l_time/r_time, 3))
      variation_ratio <- round(l_magnitude_sd$l_magnitude_sd/r_magnitude_sd$r_magnitude_sd, 3)
      simple_magnitude_ratio <- round(l_magnitude$l_magnitude/r_magnitude$r_magnitude, 3)
      
    } else {
      use_ratio <- Hz1DurationMetrics %>% dplyr::summarise(use_ratio = round(r_time/l_time, 3))
      variation_ratio <- round(r_magnitude_sd$r_magnitude_sd/l_magnitude_sd$l_magnitude_sd, 3)
      simple_magnitude_ratio <- round(r_magnitude$r_magnitude/l_magnitude$l_magnitude, 3)
      
    } 

    
    #Calculating ntropy
    #Create new data frame for entropy
    Day1LVMData <- Combined1Hz %>% select(LVMData) %>% dplyr::rename("LVMData_Day1" = "LVMData")
    Day1RVMData <- Combined1Hz %>% select(RVMData) %>% dplyr::rename("RVMData_Day1" = "RVMData")

    #Prior to calculating entropy, need to find the most active hour on most active side.
    if (AffSide == "Right") {
      #Use left side to find length of dataset
      LENG_DAY1 <- length(Day1LVMData$LVMData_Day1)
      EntropyData_Day1 <- Day1LVMData$LVMData_Day1

    } else {
      #Use right side to find length of data
      LENG_DAY1 <- length(Day1RVMData$RVMData_Day1)
      EntropyData_Day1 <- Day1RVMData$RVMData_Day1

    }
    
    #Subtract the last hour because we will only compare whole hours
    LENG.LAST.Day1 <- LENG_DAY1 - (60 * 60) 
    
    #also set the length of an hour in terms of frames - here for 1 sec sampling rate
    HOUR.SAMPLINGRATE <- (60 * 60)
    
    
    ###DAY 1###
    
    ##Loop
    #Set some staring values and vector for the loop below
    X=0
    MAX = 0
    MIN = 120
    SAMPLE.AVE <- vector()
    
    #For loop to run through each full hour of the data
    #LENG.LAST is the length of the data set up to the last whole hour.
      for (X in 1:LENG.LAST.Day1){
        #HOUR = the frames in one hour based on the capture hertz
        END <- X + (HOUR.SAMPLINGRATE-1)
        
        SAMPLE <- EntropyData_Day1[X:END]
        
        MEAN <- mean(SAMPLE, na.rm = TRUE)
        SAMPLE.AVE[X] <- MEAN
        
        #If loop to store the maximum average (most active hour)
        if (SAMPLE.AVE[X] > MAX){
          MAX <- SAMPLE.AVE[X]
          START.MAX <- X
        }
      }
      
      #Find the END row numbers based on the starting row
      END.MAX <- START.MAX + HOUR.SAMPLINGRATE
      RVM.MAX_Day1 <- Day1RVMData$RVMData_Day1[START.MAX:END.MAX]
      LVM.MAX_Day1 <- Day1LVMData$LVMData_Day1[START.MAX:END.MAX]
      
      #Calculate entropy on just the hour of maximum activity
      r_entropy <- round(sample_entropy(RVM.MAX_Day1, edim = 2, 0.2*sd(RVM.MAX_Day1), tau = 1) ,3)
      l_entropy <- round(sample_entropy(LVM.MAX_Day1, edim = 2, 0.2*sd(LVM.MAX_Day1), tau = 1) ,3)


    
    ##Load in associated 30Hz Files###
    #Load left 30 Hz file
    L_30Hz <- list.files(L_dir, pattern = "RAW.csv")
    L_30Hz <- paste(L_dir, L_30Hz, sep = "")
    L_30Hz_AccelData <- read.csv(L_30Hz, skip = 10)
    names(L_30Hz_AccelData) <- c("X", "Y", "Z")
    
    #Load right 30 Hz file and add time stamp
    R_30Hz <- list.files(R_dir, pattern = "RAW.csv")
    R_30Hz <- paste(R_dir, R_30Hz, sep = "") 
    R_30Hz_AccelData <- read.csv(R_30Hz, skip = 10)
    names(R_30Hz_AccelData) <- c("X", "Y", "Z")
    
    ##If left and right sensors have different recording lengths, set LENG_30 to be the shorter of the two:
    if(nrow(L_30Hz_AccelData) == nrow(R_30Hz_AccelData)){
      LENG_30 <- nrow(L_30Hz_AccelData)
    }
    
    if(nrow(L_30Hz_AccelData) < nrow(R_30Hz_AccelData)){
      LENG_30 <- nrow(L_30Hz_AccelData)
      R_30Hz_AccelData <- R_30Hz_AccelData[1:LENG_30, ]
    }
    
    if(nrow(R_30Hz_AccelData) < nrow(L_30Hz_AccelData)){
      LENG_30 <- nrow(R_30Hz_AccelData)
      L_30Hz_AccelData <- L_30Hz_AccelData[1:LENG_30, ]
    }

    #Create time stamp and add to left and right 30 Hz files
    HERTZ_30 <- 30
    TIME_STAMP <- seq(ISOdate(LStartYear,LStartMonth,LStartDay,LStartHour,LStartMinute,LStartSecond), by = (1/HERTZ_30), length.out = LENG_30)
    TIME_STAMP <- as.character(TIME_STAMP)
    ORIG_COUNT <- seq(LENG_30)
    L_30_DONE_FILE = data.frame(ORIG_COUNT,TIME = TIME_STAMP, L_30Hz_AccelData)
    R_30_DONE_FILE = data.frame(ORIG_COUNT,TIME = TIME_STAMP, R_30Hz_AccelData)

    #Remove the first 1.5 hours of the recording due to in-lab time
    InLabTime_30Hz <- 1.5 * 108000
    L_30_DONE_FILE <- L_30_DONE_FILE[-c(1:InLabTime_30Hz), ]
    R_30_DONE_FILE <- R_30_DONE_FILE[-c(1:InLabTime_30Hz), ]
    
    #Cut file if > 24 hours
    if(nrow(L_30_DONE_FILE) > 2592000){
      L_30_DONE_FILE <- L_30_DONE_FILE[1:2592000, ]
      R_30_DONE_FILE <- R_30_DONE_FILE[1:2592000, ]
    }
    

    #####CALCULATE 30 HZ VARIABLES#####
    #Compute vector magnitudes for the left and right sides
    L_30_XData <- L_30_DONE_FILE[,"X"]
    L_30_YData <- L_30_DONE_FILE[,"Y"]
    L_30_ZData <- L_30_DONE_FILE[,"Z"]
    L_30_DONE_FILE$VM <- sqrt(L_30_XData^2 + L_30_YData^2 + L_30_ZData^2)
    L_30_DONE_FILE$VM.STD <- as.vector( scale(L_30_DONE_FILE$VM))
    
    R_30_XData <- R_30_DONE_FILE[,"X"]
    R_30_YData <- R_30_DONE_FILE[,"Y"]
    R_30_ZData <- R_30_DONE_FILE[,"Z"]
    R_30_DONE_FILE$VM <-sqrt(R_30_XData^2 + R_30_YData^2 + R_30_ZData^2 )
    R_30_DONE_FILE$VM.STD <- as.vector( scale(R_30_DONE_FILE$VM))
    
    #Separate into days (only one recording day in this study)
    Day1_L_30Hz <- L_30_DONE_FILE
    Day1_R_30Hz <- R_30_DONE_FILE

    #Filter data
    bf <- butter(2, c(I(.2/15),I(12/15)), type="pass")
    Day1_L_30Hz$VM<- filtfilt(x=Day1_L_30Hz$VM, filt = bf)
    Day1_R_30Hz$VM<- filtfilt(x=Day1_R_30Hz$VM, filt = bf)
  
    #Use filtered data for Jerk calculations
    L_30Hz_Filtered <- Day1_L_30Hz
    R_30Hz_Filtered <- Day1_R_30Hz
    
    
    #Calculating jerk- a measure of smoothness of movement
    TP = 1/30 #the sampling interval
    
    ##Right jerk
    JERK.R <- vector()
    JERK.R[1]=0
    i <-1
    
    for (i in 1:(nrow(R_30Hz_Filtered)-1)){
      A=as.numeric(R_30Hz_Filtered[i+1,"VM"] )
      B=as.numeric(R_30Hz_Filtered[(i),"VM"] )
      JERK.R[i] <- (A-B) / TP
      
    }
    
    JERK.R[nrow(R_30Hz_Filtered)] = 0
    R_30Hz_Filtered$JERK.R <- JERK.R
    
    #Left jerk
    JERK.L <- vector()
    JERK.L[1]=0
    i <- 1
    
    for (i in 1:(nrow(L_30Hz_Filtered)-1)){
      A=L_30Hz_Filtered[i+1,"VM"] 
      B=L_30Hz_Filtered[(i),"VM"]
      JERK.L[i] = (A-B) / TP
      
    }
    
    JERK.L[nrow(L_30Hz_Filtered)] = 0
    L_30Hz_Filtered$JERK.L <- JERK.L
    
    #Jerk output variables
    r_jerk_ave <- R_30Hz_Filtered %>% dplyr::filter(JERK.R != 0) %>% dplyr::summarise(r_jerk_ave = round(mean(abs(JERK.R)), 3))
    l_jerk_ave <- L_30Hz_Filtered %>% dplyr::filter(JERK.L != 0) %>% dplyr::summarise(l_jerk_ave = round(mean(abs(JERK.L)), 3))
    
    #Combining outputs from the left and right sides into one data frame
    Jerk_Avg_Combined <- merge(r_jerk_ave, l_jerk_ave)
    Jerk_Avg_Combined <- as.data.frame(Jerk_Avg_Combined)
    
    #Calculating jerk asymmetry
    if (AffSide == "Left") {
      jerk_aym <- Jerk_Avg_Combined %>% dplyr::summarise(jerk_aym = (l_jerk_ave - r_jerk_ave) / (l_jerk_ave + r_jerk_ave))
    } else {
      jerk_aym <- Jerk_Avg_Combined %>% dplyr::summarise(jerk_aym = (r_jerk_ave - l_jerk_ave) / (l_jerk_ave + r_jerk_ave))
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
    SPEC.L.VM.D1 <- spectrum(Day1_L_30Hz$VM, log='no', span=5, plot = F)
    SPEC.R.VM.D1 <- spectrum(Day1_R_30Hz$VM, log='no', span=5, plot = F)
    
    #divide spectrum by sampling interval to create per time instead of per sampling interval
    freq.L.1 <- SPEC.L.VM.D1$freq/0.033333
    freq.R.1 <- SPEC.R.VM.D1$freq/0.033333
    
    #Multiplying the spectral density by 2 so that the area under the periodogram equals the variance of the time series
    Density.L.1 <- 2*SPEC.L.VM.D1$spec
    Density.R.1 <- 2*SPEC.R.VM.D1$spec
    
    #weighted variance of the frequencies present
    VAR.WT.L.1 <- weighted.var(freq.L.1,Density.L.1) 
    VAR.WT.R.1 <- weighted.var(freq.R.1,Density.R.1) 
    
    l_mean_freq <- round(weighted.mean(freq.L.1, Density.L.1), 3)
    r_mean_freq <- round(weighted.mean(freq.R.1, Density.R.1), 3)
    
    l_sd_freq <- round(sqrt(VAR.WT.L.1), 3)
    r_sd_freq <- round(sqrt(VAR.WT.R.1), 3)
    
    #Combining outputs from the left and right sides
    WtMeanFreqOutput <- data.frame(l_mean_freq = l_mean_freq, r_mean_freq = r_mean_freq)
    WtSDFreqOutput <- data.frame(l_sd_freq = l_sd_freq, r_sd_freq = r_sd_freq)


    #####SAVE TO FILES#####
    STUDYPROTOCOL <- "Dose"
    #Save to source data participant's folder
    WriteCSVFilename <- paste0(time_dir, sub_id, "_", ProcessedTimePoint, "_", STUDYPROTOCOL, "_SleepIncOutputReprocessed.csv", sep = "") 
   
    TaskData <- cbind(Hz1DurationMetrics, l_magnitude, r_magnitude, bilateral_magnitude, l_magnitude_sd,
                      r_Magnitude_sd, l_peak_magnitude, r_peak_magnitude, use_ratio, variation_ratio,
                      simple_magnitude_ratio, r_entropy, l_entropy, Jerk_Avg_Combined, jerk_aym,
                      WtMeanFreqOutput, WtSDFreqOutput)
    
    TaskData <- data.frame(StudyID = sub_id, TimePoint = ProcessedTimePoint, TaskData)
    
    #Save data to csv file in source data participant's folder
    write.csv(TaskData, file = WriteCSVFilename, row.names = F)
    
}


