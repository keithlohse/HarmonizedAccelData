#This code script was used to process upper limb accelerometry data for the study, A Feasibility Study of Bilateral Wrist Sensors for
#Measuring Motor Traits in Children with Autism (PMC9974780) for the project:
#Harmonized Upper and Lower Limb Accelerometry Data
#Accelerometry data for this cohort were processed individually and one day at a time.
#Participants wore bilateral wrist sensors for 12 hours on any two days.

##### This chunk was used to "manually" combine recorded movement into a 12 hour block (when possible) #######
RIGHT_DATA <- read.csv(R.FILE, skip = 0, header = T)

RIGHT_DATA <- RIGHT_DATA[(1:32000),]
RIGHT_DATA <- RIGHT_DATA[-(15500:18500),]

V.MAG.R <- sqrt(RIGHT_DATA$X^2 + RIGHT_DATA$Y^2 + RIGHT_DATA$Z^2)
plot(V.MAG.R,type="l", col="black", 
     ylab="R Activity Counts", xlab = "Time(sec)")

LEFT_DATA <- read.csv(L.FILE, skip = 0, header = T)
LEFT_DATA <- LEFT_DATA[(1:32000),]
LEFT_DATA <- LEFT_DATA[-(15500:18500),]


V.MAG.L <- sqrt(LEFT_DATA$X^2 + LEFT_DATA$Y^2 + LEFT_DATA$Z^2)
plot(V.MAG.L,type="l", col="black", 
     ylab="L Activity Counts", xlab = "Time(sec)")
######################################################################
library(miscTools)
library(Rmisc)
library(cowplot)
library(varhandle)
library(ggplot2)
library(psych)
library(lsr)
library(ggrepel)
library (car)
library(irr)
library(svDialogs)
library(forecast)
library(tcltk)
library(shiny)
library(GGally)
library(gridExtra)
library(gWidgets2)
library(zoo)
library(stats)
library(tidyverse) # suite of ggplot2,tibble,tidyr,purrr,dplyr,stringr,forcats
library(magrittr)  # for extra pipes %$% %<>%
library(labelled) # for working with variable labels (e.g. var_label() & val_label())
library(ggthemes)
library(readxl)
library(graphics)
library(MESS)
library(phonTools)
library(ggpubr)
library(pracma)

remove(list = ls())

CONTINUE = T
while(CONTINUE == T){
########################################################## GET INFO FROM USER ##############
SUBJECTID.TEMP <-  dlgInput(message = "Enter SUBJECTID ID (C### OR A### or S###) ", default = NULL, gui = .GUI) 
SUBJECTID <-as.character(SUBJECTID.TEMP$res)

HANDEDNESS.TEMP <-  dlgInput(message = "Enter Handedness (R OR L) ", default = "R", gui = .GUI) 
HANDEDNESS <-as.character(HANDEDNESS.TEMP$res)
Handedness <- HANDEDNESS # Maintains compatibility with previous coded, may be removed in the future


DATE.TEMP <-  dlgInput(message = "Enter the date for the file (yearmonthday).", gui = .GUI) 
DATE <-as.character(DATE.TEMP$res)

AGE.TEMP <-  dlgInput(message = paste("Enter the age in months for",SUBJECTID), default="126", gui = .GUI) 
AGE <-as.character(AGE.TEMP$res)
Age = AGE # Maintains compatibility with previous coded, may be removed in the future


VISIT <- dlgInput(message = "Enter visit type. \nEither: wkday, wkend (or day1, day2)")
VISIT <- VISIT$res
Visit = VISIT # Maintains compatibility with previous coded, may be removed in the future

GENDER <- dlgInput(message = "Male or Female: \nEnter M or F",default="M")
GENDER <- GENDER$res

GROUP <- dlgInput(message = "Enter Group: \nControl = C, Autism = A, Sibling = S")
GROUP <- GROUP$res
#############################################




############READING IN DATA##############
dlg_message("Please load the Right csv file.", gui = .GUI)
R.FILE <- dlg_open(title = "Select csv Accelerometry File", filters = dlg_filters[c("All"),])$res
RIGHT_DATA <- read.csv(R.FILE, skip = 0, header = T)

dlg_message("Please load the Left csv file.", gui = .GUI)
L.FILE <- dlg_open(title = "Select csv Accelerometry File", filters = dlg_filters[c("All"),])$res
LEFT_DATA <- read.csv(L.FILE, skip = 0, header = T)


#Pick path for saving
Sys.sleep(2)
setwd(dlg_dir(title = "Set directory for saving Rdata and csv files", default = getwd())$res)
# User prompted to enter file name, age, and, handedness and store those---------------------------------------------------------------------------------------
dlg_message("Remember to enter the 8am to 8pm start/stop rows below.")

break  #stops the code to make sure the correct row numbers are in, can then restart from here

##########################################
#HERE MANUALLY ENTERED ROWS TO KEEP WITHIN 800am - 800pm 
dlg_message("Enter the 8am to 8pm rows below.")
#may have to subset in different ways depending on the period recorded.
#enter the START row here
START <- 	1
#STOP <- 	42030
#STOP <- START + 12*(60*60)
STOP <- nrow(RIGHT_DATA)
RIGHT_DATA <- data.frame(RIGHT_DATA[START:STOP,])
LEFT_DATA <- data.frame(LEFT_DATA[START:STOP,])



###CALC VECTOR MAGNITUE AND ADD IT TO THE DATAFRAME#####
V.MAG.R <- sqrt(RIGHT_DATA$X^2 + RIGHT_DATA$Y^2 + RIGHT_DATA$Z^2)
RIGHT_DATA <- data.frame(RIGHT_DATA, V.MAG.R)

V.MAG.L <- sqrt(LEFT_DATA$X^2 + LEFT_DATA$Y^2 + LEFT_DATA$Z^2)
LEFT_DATA <- data.frame(LEFT_DATA, V.MAG.L)


if (HANDEDNESS == 'R') {V.MAG.INDEX <- V.MAG.R
}else {V.MAG.INDEX <- V.MAG.L}
########################################


############## PLOT TO HELP REMOVE SLEEP TIME ##############

ggplot() + aes(V.MAG.INDEX) + geom_density()


####LENGTH OF DATA SET#####
LENG <- length(V.MAG.INDEX)
if (LENG %% 2 == 0 ){LENG <- LENG - 1} #Creates odd number rows to get an exact median
#############################

######## SET UP TO STOP AFTER THE LAST HOUR #############
#Subtract the last hour because we will only compare whole hours
LENG.LAST <- LENG - (60 * 60) 


#also set the length of an hour in terms of frames - here for 1 sec sampling rate
HOUR <- (60 * 60)

##################################################


######### Loop ###################
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
  SAMPLE <- V.MAG.INDEX[X:END]
  MEAN <- mean(SAMPLE)
  SAMPLE.AVE[X] <- MEAN
  
  #If loop to store the maximum average (most active hour)
  if (SAMPLE.AVE[X] > MAX){
    MAX <- SAMPLE.AVE[X]
    START.MAX <- X
  }
  #If loop to store the minimum average (least active hour)
  if (SAMPLE.AVE[X] < MIN){
    MIN <- SAMPLE.AVE[X]
    START.MIN <- X
  }
  
}
####################################################

MEDIAN <- median(SAMPLE.AVE)
#LOOP TO FIND THE MEDIAN HOURS START ROW, WHICH IS X
for (X in 1:LENG.LAST){
  if(SAMPLE.AVE[X] == MEDIAN){
    START.MEDIAN <- X
  }
}


#Find the END row numbers based on the starting row
END.MEDIAN <- START.MEDIAN + HOUR

END.MAX <- START.MAX + HOUR

END.MIN <-START.MIN + HOUR



###################### NOW START CALCULATIONS FOR ONE SECOND VARIABLES ######



#set theme to minimal
theme_set(theme_minimal()) #all plots that run below use this theme

#PREPARES THE FILES NAMES FOR SAVING:
WRITE_CSV_FILENAME <- paste(SUBJECTID, "_", Visit, "_1sec_BY_ACTIVITY_HOURS.csv", sep = "") 
WRITE_RDATA_FILENAME <- paste(SUBJECTID, "_", Visit, "_1sec_BY_ACTIVITY_HOURS.RData", sep = "") 
#Note: output variables rounded to 3 digits after the decimal to read into REDCap
#---------------------------------------------------------------------


#Plot to check vector magnitudes, making sure two sides ~match
plot(V.MAG.R, type="l", col="blue", 
     ylab="Activity Counts", xlab = "Time(sec)")
matplot(V.MAG.L, type="l", col="red", add = 1)
Sys.sleep(2)

#remove(LEFT.DATA, RIGHT.DATA) #clean up workspace

#BREAK OUT HOURS BY ACTIVITY:
MAX.HOUR.R <- V.MAG.R[START.MAX:END.MAX]
MAX.HOUR.L <- V.MAG.L[START.MAX:END.MAX]

MED.HOUR.R <- V.MAG.R[START.MEDIAN:END.MEDIAN]
MED.HOUR.L <- V.MAG.L[START.MEDIAN:END.MEDIAN]

MIN.HOUR.R <- V.MAG.R[START.MIN:END.MIN]
MIN.HOUR.L <- V.MAG.L[START.MIN:END.MIN]


############## CALCULATIONS FOR MAX HOUR

HOUR.TYPE <- "MAX"


##############THRESHOLD#########################################
#Find values > THRESHOLD.  
THRESHOLD <- 2 
L.COUNT <- which(MAX.HOUR.L >= THRESHOLD)
R.COUNT <- which(MAX.HOUR.R >= THRESHOLD) 


######################################################## Variable calculations below
#Unused variable commented out:

####START WITH MAX HOUR ######################
#Duration variables
WEAR.TIME.MAX <- round((length(MAX.HOUR.L)/3600), 0) #Wearing times in hrs, as integer
LHrs.MAX <- round((length(L.COUNT)/3600), 3)
RHrs.MAX <- round((length(R.COUNT)/3600), 3)
NoMvt.MAX <- which(MAX.HOUR.L == 0 & MAX.HOUR.R == 0)

TotalNoMvtHrs.MAX <- round((length(NoMvt.MAX))/3600, 3)

TotalMvtHrs.MAX <- round(WEAR.TIME.MAX - TotalNoMvtHrs.MAX, 3)

TempData <- which(MAX.HOUR.L > THRESHOLD & MAX.HOUR.R > THRESHOLD)

SimultaneousHrs.MAX <- round(length(TempData)/3600, 3)

TempData <- which(MAX.HOUR.L > THRESHOLD & MAX.HOUR.R < THRESHOLD) 

LOnlyHrs.MAX <- round(length(TempData)/3600, 3)

TempData <- which(MAX.HOUR.L < THRESHOLD & MAX.HOUR.R > THRESHOLD) 

ROnlyHrs.MAX <- round(length(TempData)/3600, 3)

#Plot Duration variables
ChartData.MAX <- data.frame(Category = c("TotalNoMvtHrs.MAX","LOnlyHrs.MAX", "ROnlyHrs.MAX",
                                         "SimultaneousHrs.MAX"),
                        value = c(TotalNoMvtHrs.MAX, LOnlyHrs.MAX, ROnlyHrs.MAX, SimultaneousHrs.MAX))
b1 <- ggplot(ChartData.MAX, mapping = aes(x="", y=value, fill = Category))
b1 <- b1 + geom_bar(width = 1, stat = "identity") + scale_fill_manual(values = c("grey90", "grey70", "grey50", "black"))
b1 <- b1 + labs(x = "Time spent moving", y = "One Hour")
b1

#To make code compatible with previously written variables
LVMData <- MAX.HOUR.L
RVMData <- MAX.HOUR.R

#Magnitude variables, in activity counts
LMagnitude.MAX <- round(sum(LVMData, na.rm = TRUE), 3)
RMagnitude.MAX <- round(sum(RVMData, na.rm = TRUE), 3)

TotalMagnitude.MAX <- round(LMagnitude.MAX + RMagnitude.MAX, 2) 

LMagnitudeSD.MAX <- round(sd(LVMData[L.COUNT]), 3) #SDs when moving
RMagnitudeSD.MAX <- round(sd(RVMData[R.COUNT]), 3)


#Ratios 
if (Handedness == "R") {
  UseRatio.MAX <- round(LHrs.MAX/RHrs.MAX, 3) 
  VariationRatio.MAX <- round(LMagnitudeSD.MAX/RMagnitudeSD.MAX, 3)
  SimpleMagnitudeRatio.MAX <- round(LMagnitude.MAX/RMagnitude.MAX, 3)  
} else {
  UseRatio.MAX <- round(RHrs.MAX/LHrs.MAX, 3) 
  VariationRatio.MAX <- round(RMagnitudeSD.MAX/LMagnitudeSD.MAX, 3)
  SimpleMagnitudeRatio.MAX <- round(RMagnitude.MAX/LMagnitude.MAX, 3)
}


SampEn.R.MAX <- round(sample_entropy(MAX.HOUR.R, edim = 2, 0.2*sd(MAX.HOUR.R), tau = 1) ,3)
SampEn.L.MAX <- round(sample_entropy(MAX.HOUR.L, edim = 2, 0.2*sd(MAX.HOUR.L), tau = 1) ,3)

RESULTS.MAX <- c(SUBJECTID, Age, GROUP, Visit, HOUR.TYPE,
                 WEAR.TIME.MAX, TotalNoMvtHrs.MAX, TotalMvtHrs.MAX, 
                 LHrs.MAX, LOnlyHrs.MAX, RHrs.MAX, ROnlyHrs.MAX, SimultaneousHrs.MAX,
                 LMagnitude.MAX, RMagnitude.MAX, TotalMagnitude.MAX,
                 LMagnitudeSD.MAX, RMagnitudeSD.MAX,
                 UseRatio.MAX, VariationRatio.MAX, SimpleMagnitudeRatio.MAX, 
                 VMCorrCoef.MAX, SampEn.R.MAX, SampEn.L.MAX)

##### START AGAIN FOR MEDIAN HOUR ####
L.COUNT <- which(MED.HOUR.L >= THRESHOLD)
R.COUNT <- which(MED.HOUR.R >= THRESHOLD) 



HOUR.TYPE <- "MED"

#Duration variables
WEAR.TIME.MED <- round((length(MED.HOUR.L)/3600), 0) #Wearing times in hrs, as integer
LHrs.MED <- round((length(L.COUNT)/3600), 3)
RHrs.MED <- round((length(R.COUNT)/3600), 3)

NoMvt.MED <- which(MED.HOUR.L == 0 & MED.HOUR.R == 0)

TotalNoMvtHrs.MED <- round((length(NoMvt.MED))/3600, 3)

TotalMvtHrs.MED <- round(WEAR.TIME.MED - TotalNoMvtHrs.MED, 3)

TempData <- which(MED.HOUR.L > THRESHOLD & MED.HOUR.R > THRESHOLD)

SimultaneousHrs.MED <- round(length(TempData)/3600, 3)

TempData <- which(MED.HOUR.L > THRESHOLD & MED.HOUR.R < THRESHOLD) 

LOnlyHrs.MED <- round(length(TempData)/3600, 3)

TempData <- which(MED.HOUR.L < THRESHOLD & MED.HOUR.R > THRESHOLD) 

ROnlyHrs.MED <- round(length(TempData)/3600, 3)

#Plot Duration variables
ChartData.MED <- data.frame(Category = c("TotalNoMvtHrs.MED","LOnlyHrs.MED", "ROnlyHrs.MED",
                                         "SimultaneousHrs.MED"),
                            value = c(TotalNoMvtHrs.MED, LOnlyHrs.MED, ROnlyHrs.MED, SimultaneousHrs.MED))
b2 <- ggplot(ChartData.MED, mapping = aes(x="", y=value, fill = Category))
b2 <- b2 + geom_bar(width = 1, stat = "identity") + scale_fill_manual(values = c("grey90", "grey70", "grey50", "black"))
b2 <- b2 + labs(x = "Time spent moving", y = "One Hour")
b2

#To make code compatible with previously written variables
LVMData <- MED.HOUR.L
RVMData <- MED.HOUR.R

#Magnitude variables, in activity counts
LMagnitude.MED <- round(sum(LVMData, na.rm = TRUE), 3)
RMagnitude.MED <- round(sum(RVMData, na.rm = TRUE), 3)

TotalMagnitude.MED <- round(LMagnitude.MED + RMagnitude.MED, 2) 

LMagnitudeSD.MED <- round(sd(LVMData[L.COUNT]), 3) #SDs when moving
RMagnitudeSD.MED <- round(sd(RVMData[R.COUNT]), 3)


#Ratios 
if (Handedness == "R") {
  UseRatio.MED <- round(LHrs.MED/RHrs.MED, 3) 
  VariationRatio.MED <- round(LMagnitudeSD.MED/RMagnitudeSD.MED, 3)
  SimpleMagnitudeRatio.MED <- round(LMagnitude.MED/RMagnitude.MED, 3)  
} else {
  UseRatio.MED <- round(RHrs.MED/LHrs.MED, 3) 
  VariationRatio.MED <- round(RMagnitudeSD.MED/LMagnitudeSD.MED, 3)
  SimpleMagnitudeRatio.MED <- round(RMagnitude.MED/LMagnitude.MED, 3)
}


#For output to RData file
#Generate dataframe to house variables

SampEn.R.MED <- round(sample_entropy(MED.HOUR.R, edim = 2, 0.2*sd(MED.HOUR.R), tau = 1),3)
SampEn.L.MED <- round(sample_entropy(MED.HOUR.L, edim = 2, 0.2*sd(MED.HOUR.L), tau = 1),3)

RESULTS.MED <- c(SUBJECTID, Age, GROUP, Visit, HOUR.TYPE,
                   WEAR.TIME.MED, TotalNoMvtHrs.MED, TotalMvtHrs.MED, 
                   LHrs.MED, LOnlyHrs.MED, RHrs.MED, ROnlyHrs.MED, SimultaneousHrs.MED,
                   LMagnitude.MED, RMagnitude.MED, TotalMagnitude.MED,
                   LMagnitudeSD.MED, RMagnitudeSD.MED,
                   UseRatio.MED, VariationRatio.MED, SimpleMagnitudeRatio.MED, 
                   VMCorrCoef.MED, SampEn.R.MED, SampEn.L.MED)


####### START AGAIN CALCULATING MINIMUM HOUR ############

L.COUNT <- which(MIN.HOUR.L >= THRESHOLD)
R.COUNT <- which(MIN.HOUR.R >= THRESHOLD) 


HOUR.TYPE <- "MIN"

#Duration variables
WEAR.TIME.MIN <- round((length(MIN.HOUR.L)/3600), 0) #Wearing times in hrs, as integer
LHrs.MIN <- round((length(L.COUNT)/3600), 3)
RHrs.MIN <- round((length(R.COUNT)/3600), 3)

NoMvt.MIN <- which(MIN.HOUR.L == 0 & MIN.HOUR.R == 0)

TotalNoMvtHrs.MIN <- round((length(NoMvt.MIN))/3600, 3)

TotalMvtHrs.MIN <- round(WEAR.TIME.MIN - TotalNoMvtHrs.MIN, 3)

TempData <- which(MIN.HOUR.L > THRESHOLD & MIN.HOUR.R > THRESHOLD)

SimultaneousHrs.MIN <- round(length(TempData)/3600, 3)

TempData <- which(MIN.HOUR.L > THRESHOLD & MIN.HOUR.R < THRESHOLD) 

LOnlyHrs.MIN <- round(length(TempData)/3600, 3)

TempData <- which(MIN.HOUR.L < THRESHOLD & MIN.HOUR.R > THRESHOLD) 

ROnlyHrs.MIN <- round(length(TempData)/3600, 3)

#Plot Duration variables
ChartData.MIN <- data.frame(Category = c("TotalNoMvtHrs.MIN","LOnlyHrs.MIN", "ROnlyHrs.MIN",
                                         "SimultaneousHrs.MIN"),
                            value = c(TotalNoMvtHrs.MIN, LOnlyHrs.MIN, ROnlyHrs.MIN, SimultaneousHrs.MIN))
b3 <- ggplot(ChartData.MIN, mapping = aes(x="", y=value, fill = Category))
b3 <- b3 + geom_bar(width = 1, stat = "identity") + scale_fill_manual(values = c("grey90", "grey70", "grey50", "black"))
b3 <- b3 + labs(x = "Time spent moving", y = "Hours")
b3

#To make code compatible with previously written variables
LVMData <- MIN.HOUR.L
RVMData <- MIN.HOUR.R

#Magnitude variables, in activity counts
LMagnitude.MIN <- round(sum(LVMData, na.rm = TRUE), 3)
RMagnitude.MIN <- round(sum(RVMData, na.rm = TRUE), 3)

TotalMagnitude.MIN <- round(LMagnitude.MIN + RMagnitude.MIN, 2) 


LMagnitudeSD.MIN <- round(sd(LVMData[L.COUNT]), 3) #SDs when moving
RMagnitudeSD.MIN <- round(sd(RVMData[R.COUNT]), 3)



#Ratios 
if (Handedness == "R") {
  UseRatio.MIN <- round(LHrs.MIN/RHrs.MIN, 3) 
  VariationRatio.MIN <- round(LMagnitudeSD.MIN/RMagnitudeSD.MIN, 3)
  SimpleMagnitudeRatio.MIN <- round(LMagnitude.MIN/RMagnitude.MIN, 3)  
} else {
  UseRatio.MIN <- round(RHrs.MIN/LHrs.MIN, 3) 
  VariationRatio.MIN <- round(RMagnitudeSD.MIN/LMagnitudeSD.MIN, 3)
  SimpleMagnitudeRatio.MIN <- round(RMagnitude.MIN/LMagnitude.MIN, 3)
}



#For output to RData file
#Generate dataframe to house variables

SampEn.R.MIN <- round(sample_entropy(MIN.HOUR.R, edim = 2, 0.2*sd(MIN.HOUR.R), tau = 1),3) 
SampEn.L.MIN <- round(sample_entropy(MIN.HOUR.L, edim = 2, 0.2*sd(MIN.HOUR.L), tau = 1),3)

record_id <- SUBJECTID
RESULTS.MIN <- c(record_id, Age, GROUP, Visit, HOUR.TYPE,
                 WEAR.TIME.MIN, TotalNoMvtHrs.MIN, TotalMvtHrs.MIN, 
                 LHrs.MIN, LOnlyHrs.MIN, RHrs.MIN, ROnlyHrs.MIN, SimultaneousHrs.MIN,
                 LMagnitude.MIN, RMagnitude.MIN, TotalMagnitude.MIN,
                 LMagnitudeSD.MIN, RMagnitudeSD.MIN,
                 UseRatio.MIN, VariationRatio.MIN, SimpleMagnitudeRatio.MIN, 
                 VMCorrCoef.MIN, SampEn.R.MIN, SampEn.L.MIN)



############################### Start again calculating for 12 hours########


L.COUNT <- which(V.MAG.L >= THRESHOLD)
R.COUNT <- which(V.MAG.R >= THRESHOLD) 


HOUR.TYPE <- "12hrs"

#Duration variables
WEAR.TIME.12 <- round((length(V.MAG.L)/3600), 0) #Wearing times in hrs, as integer
LHrs.12 <- round((length(L.COUNT)/3600), 3)
RHrs.12 <- round((length(R.COUNT)/3600), 3)

NoMvt.12 <- which(V.MAG.L == 0 & V.MAG.R == 0)

TotalNoMvtHrs.12 <- round((length(NoMvt.12))/3600, 3)

TotalMvtHrs.12 <- round(WEAR.TIME.12 - TotalNoMvtHrs.12, 3)

TempData <- which(V.MAG.L > THRESHOLD & V.MAG.R > THRESHOLD)

SimultaneousHrs.12 <- round(length(TempData)/3600, 3)

TempData <- which(V.MAG.L > THRESHOLD & V.MAG.R < THRESHOLD) 

LOnlyHrs.12 <- round(length(TempData)/3600, 3)

TempData <- which(V.MAG.L < THRESHOLD & V.MAG.R > THRESHOLD) 

ROnlyHrs.12 <- round(length(TempData)/3600, 3)

#Plot Duration variables
ChartData.12 <- data.frame(Category = c("TotalNoMvtHrs.12","LOnlyHrs.12", "ROnlyHrs.12",
                                         "SimultaneousHrs.12"),
                            value = c(TotalNoMvtHrs.12, LOnlyHrs.12, ROnlyHrs.12, SimultaneousHrs.12))
b4 <- ggplot(ChartData.12, mapping = aes(x="", y=value, fill = Category))
b4 <- b4 + geom_bar(width = 1, stat = "identity") + scale_fill_manual(values = c("grey90", "grey70", "grey50", "black"))
b4 <- b4 + labs(x = "Time spent moving", y = "Hours")
b4

#To make code compatible with previously written variables
LVMData <- V.MAG.L
RVMData <- V.MAG.R

#Magnitude variables, in activity counts
LMagnitude.12 <- round(sum(LVMData, na.rm = TRUE), 3)
RMagnitude.12 <- round(sum(RVMData, na.rm = TRUE), 3)

TotalMagnitude.12 <- round(LMagnitude.12 + RMagnitude.12, 2) 


LMagnitudeSD.12 <- round(sd(LVMData[L.COUNT]), 3) #SDs when moving
RMagnitudeSD.12 <- round(sd(RVMData[R.COUNT]), 3)



#Ratios 
if (Handedness == "R") {
  UseRatio.12 <- round(LHrs.12/RHrs.12, 3) 
  VariationRatio.12 <- round(LMagnitudeSD.12/RMagnitudeSD.12, 3)
  SimpleMagnitudeRatio.12 <- round(LMagnitude.12/RMagnitude.12, 3)  
} else {
  UseRatio.12 <- round(RHrs.12/LHrs.12, 3) 
  VariationRatio.12 <- round(RMagnitudeSD.12/LMagnitudeSD.12, 3)
  SimpleMagnitudeRatio.12 <- round(RMagnitude.12/LMagnitude.12, 3)
}



#For output to RData file
#Generate dataframe to house variables

#FOR NOW NOT CALCULATING THE 12 hour sample entropy because it takes ~2 hours
SampEn.R.12 <- "NA"
  #round(sample_entropy(V.MAG.R, edim = 2, 0.2*sd(V.MAG.R), tau = 1),3) 
SampEn.L.12 <- "NA"
  #round(sample_entropy(V.MAG.L, edim = 2, 0.2*sd(V.MAG.L), tau = 1),3)

RESULTS.12 <- c(record_id, Age, GROUP, Visit, HOUR.TYPE,
                 WEAR.TIME.12, TotalNoMvtHrs.12, TotalMvtHrs.12, 
                 LHrs.12, LOnlyHrs.12, RHrs.12, ROnlyHrs.12, SimultaneousHrs.12,
                 LMagnitude.12, RMagnitude.12, TotalMagnitude.12,
                 LMagnitudeSD.12, RMagnitudeSD.12,
                 UseRatio.12, VariationRatio.12, SimpleMagnitudeRatio.12, 
                 VMCorrCoef.12, SampEn.R.12, SampEn.L.12)




#Put the results from each hour (Max, Med, Min) into one data frame
Acc1SecData <- data.frame(rbind(RESULTS.12,RESULTS.MAX,RESULTS.MED,RESULTS.MIN),GENDER)




#For output to csv file
#Rename output columns to exactly match REDCap variable names
colnames(Acc1SecData) <- c("record_id", "Age", "Group", "Visit", "Hour.Type",
                       "Wear.Time", "TotalNoMvtHrs", "TotalMvtHrs", 
                        "LHrs", "LOnlyHrs", "RHrs", "ROnlyHrs", "SimultaneousHrs",
                        "LMagnitude", "RMagnitude", "TotalMagnitude",
                        "LMagnitudeSD", "RMagnitudeSD",
                        "UseRatio", "VariationRatio", "SimpleMagnitudeRatio", 
                        "VMCorrCoef", "SampleEntropy.R","SampleEntropy.L", "Gender")
#Changing values back to numeric from characters
Acc1SecData$TotalMagnitude<-as.numeric(Acc1SecData$TotalMagnitude)
Acc1SecData$Age<-as.numeric(Acc1SecData$Age)
Acc1SecData$Wear.Time<-as.numeric(Acc1SecData$Wear.Time)
Acc1SecData$TotalNoMvtHrs<-as.numeric(Acc1SecData$TotalNoMvtHrs)
Acc1SecData$TotalMvtHrs<-as.numeric(Acc1SecData$TotalMvtHrs)
Acc1SecData$LHrs<-as.numeric(Acc1SecData$LHrs)
Acc1SecData$RHrs<-as.numeric(Acc1SecData$RHrs)
Acc1SecData$LOnlyHrs<-as.numeric(Acc1SecData$LOnlyHrs)
Acc1SecData$ROnlyHrs<-as.numeric(Acc1SecData$ROnlyHrs)
Acc1SecData$SimultaneousHrs<-as.numeric(Acc1SecData$SimultaneousHrs)
Acc1SecData$LMagnitude<-as.numeric(Acc1SecData$LMagnitude)
Acc1SecData$LMagnitudeSD<-as.numeric(Acc1SecData$LMagnitudeSD)
Acc1SecData$RMagnitude<-as.numeric(Acc1SecData$RMagnitude)
Acc1SecData$RMagnitudeSD<-as.numeric(Acc1SecData$RMagnitudeSD)
Acc1SecData$UseRatio<-as.numeric(Acc1SecData$UseRatio)
Acc1SecData$VariationRatio<-as.numeric(Acc1SecData$VariationRatio)
Acc1SecData$SimpleMagnitudeRatio<-as.numeric(Acc1SecData$SimpleMagnitudeRatio)
Acc1SecData$VMCorrCoef<-as.numeric(Acc1SecData$VMCorrCoef)
Acc1SecData$SampleEntropy.R<-as.numeric(Acc1SecData$SampleEntropy.R)
Acc1SecData$SampleEntropy.L<-as.numeric(Acc1SecData$SampleEntropy.L)



Sample.Entropy.Data <- data.frame("Subject" = SUBJECTID,
                                  "SampEn" = c(SampEn.R.MAX, SampEn.L.MAX, SampEn.R.MED, 
                                               SampEn.L.MED,SampEn.R.MIN, SampEn.L.MIN),
                                  "Side" = c('R','L','R','L','R','L'),
                                  "Hour.Type" = c('Max', 'Max', 'Med', 'Med', 'Min', 'Min'), 
                                  "Gender" = GENDER, "Visit" = VISIT)


Sys.sleep(3)
#Save data to R data file
save(Acc1SecData, file = WRITE_RDATA_FILENAME, row.names = F) #Save the main data frame

#Save the data frame from each hour (for both sides) all into one RData file:
save(MAX.HOUR.R, MAX.HOUR.L, MED.HOUR.R, MED.HOUR.L, MIN.HOUR.R, MIN.HOUR.L, RIGHT_DATA, LEFT_DATA,
     file = paste("Activity_Hours_", Visit,".RData", sep = ""))

#Save the main data frame as csv
write.csv(Acc1SecData, file = WRITE_CSV_FILENAME, row.names = F) 

#Save the sample entropy data as csv
write.csv(Sample.Entropy.Data,
          file = paste(SUBJECTID, "_",Visit,"_Sample_Entropy.csv", sep = ""), row.names = F)

jpeg(paste("Duration.max.",SUBJECTID,"_",VISIT,".jpg", sep = ""))
b1
dev.off()

jpeg(paste("Duration.med.",SUBJECTID,"_",VISIT,".jpg", sep = ""))
b2
dev.off()

jpeg(paste("Duration.min.",SUBJECTID,"_",VISIT,".jpg", sep = ""))
b3
dev.off()

jpeg(paste("Duration.12.",SUBJECTID,"_",VISIT,".jpg", sep = ""))
b4
dev.off()
}
