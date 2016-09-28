########################################################################################################
# Process to generate a NetCDF file with an analysis of the vertical data (L1.2 product)
# Project: SOCIB Seaturtle
# Author: Chloé Dalleau
# Date: 08/08/2016
# Based on the script of David March
#
# Description :
# This script is the third step of a global processing about the analysis of the turtles' trajectory.
# The following steps create a NetCDF file with the analysis of the vertical data:
# - Import NetCDF L1 product
# - Use dive Move
# - Export in a NetCDF file
#
# WARNING:
# - zero-offset correction (ZOC) provides 3 methods: a) visual (as current version), b) setting a global obset, c) using smoothing filter. Eddy seems to present
# an offset that varies with the time. So, method b) should not be valid. Method a) requires using a GUI that works no so good. So, explore method c) with Eddy, and even contact the author of the package.
#
#
# References:
# Diving  Behaviour Analysis in R. Luque 2007
########################################################################################################

### Remove the previous data
rm(list=ls()) 

### Import libraries and sources
source( 'ext/lib/ncParse.R')
source("src/utils.R")
library(RNetCDF)
library(ncdf4)
library(stats4)
library(diveMove)
library(lubridate)
library(dplyr) # intersect
library(sp) # merge

###################################### Step 1: Import data ############################################

### Set the turtle ID
tag_id <-  00000 # modify the tag_id according to your turtle 

### Import and convert the NetCDF L1 to a list
data_type <- "netcdf"
level <- "L1_1"
file <- paste("data/output/",tag_id,"/",level,"/",tag_id,"-",level,".nc", sep="")
dataset <- ncParse(file)

### Export interesting variables in a data frame
time <- dataset$dimensions$time$data[which(is.na(dataset$variables$depth$data)==F)]
time <- parse_date_time(time, c("HMS dbY", "Ymd HMS"), locale=Sys.setlocale("LC_TIME", "English"), tz="GMT")
depth <- dataset$variables$depth$data[which(is.na(dataset$variables$depth$data)==F)]
depth <- as.numeric(depth)
quality_d <- dataset$variables$depth_range$data[which(is.na(dataset$variables$depth$data)==F)]
csvfile <- dataset$variables$source_series$data[which(is.na(dataset$variables$depth$data)==F)]
error <- dataset$variables$errorD$data[which(is.na(dataset$variables$depth$data)==F)]
data <- data.frame(quality_d, csvfile, error)



#######################################################################################################

###################################### Step 2: Analysis of data #######################################

### Creation of TDR (time-series of depth readings) file 
TDR <- createTDR(time = time, depth = depth, concurrentData = data, dtime = 300, file = file, speed=FALSE ) # for dtime(sampling interval): 300 seconds = 5 min

### Identification of activies at various scales
## Method of zero-offset correction:
# (1) in the diveMove tab select a time window with a unique surface depth value
# (2) click on "Zero Offset Correct a range"
# (3) on the plot select the start time
# (4) on the plot select the end time
# This procedure can be repeated as many times as needed (if there are overlapping time windows, the last one prevails)
# (5) click on the button "Quit" on the diveMove tab
dcalib <-calibrateDepth(TDR) 
#dcalib <-calibrateDepth(TDR, zoc.method="filter", k=c(3, 5760), probs=c(0.5, 0.02)) 
## all the parameters for calibrateDepth :
# calibrateDepth(TDR,offset =3,wet.thr=70,dry.thr=3610,dive.thr=4,descent.crit.q=0.1,ascent.crit.q=0.1,wiggle.tol=0.8) # offset in meter

## To view of dcalib
show(dcalib) # give such as the number of dry phases and aquatic phases

### Extract data from decalib
## dive activity: L:dry, W:wet, U:underwater, D:diving, Z:brief wet 
dive <- getDAct(dcalib, "dive.activity")
dive_act <- as.character(dive)
legend_dive_act <- "legend = L: dry, W:wet, U:underwater, D:diving, Z: brief wet"
## gross dive activity L: dry, W:wet, U: underwater, D: diving, Z: brief wet periods
gross_dive_act <- dcalib@gross.activity$activity
gross_dive_act <- as.character(gross_dive_act)
## dive phases: D:descent, DB:descent/bottom, B:bottom, BA:bottom/ascent, A:ascent, X:surface
dive_phases <- dcalib@dive.phases
dive_phases <- as.character(dive_phases)
legend_dive_phases <- "legend = D:descent, DB:descent/bottom, B:bottom, BA:bottom/ascent, A:ascent, X:surface"
### create a data frame with all the information from dcalib 
## add NA for the other time observations 
extract_data <- cbind(time,gross_dive_act,dive_act,dive_phases)
ini_time <- data.frame(dataset$dimensions$time$data)
colnames(ini_time) <- c( "time")
ini_time$time <- parse_date_time(ini_time$time, c("HMS dbY", "Ymd HMS"), locale=Sys.setlocale("LC_TIME", "English"), tz="GMT")
decalib_data <- merge(ini_time,extract_data, all=TRUE)

### Data frame with the final summary for each dive
tdrXSumm1 <- diveStats(dcalib)
tBudget<-timeBudget(dcalib, ignoreZ = FALSE)
tBudget$phase.no <- as.numeric(tBudget$phase.no)
dimph <- dim(tBudget)[1]

### Data frame for wet data
trip.labs <- stampDive(dcalib,ignoreZ = FALSE)
if (dim(trip.labs)[1]!=dim(tdrXSumm1)[1]){
  size <- dim(trip.labs)[1]-dim(tdrXSumm1)[1]
  if (size>0){
    trip.labs <- trip.labs[1:dim(tdrXSumm1)[1],]
  }else {
    tdrXSumm1 <- tdrXSumm1[1:dim(trip.labs)[1],]
  }
}
tdrXSumm2 <- data.frame(trip.labs,tdrXSumm1,c(1:dim(tdrXSumm1)[1]))
colnames(tdrXSumm2) <- c(colnames(trip.labs),colnames(tdrXSumm1),"dive")
dimdiv <- dim(tdrXSumm2)[1]

# ## Export in a csv file
level <- "L1_2"
output.file1 <- paste("data/output/",tag_id,"/",level,"/",tag_id,"-",level,"-activity.csv", sep= "")
output.file2 <- paste("data/output/",tag_id,"/",level,"/",tag_id,"-",level,"-dive.csv", sep= "")
write.table(decalib_data, output.file1, row.names=FALSE, sep=",", dec=".")
write.table(tdrXSumm2, output.file2, row.names=FALSE, sep=",", dec=".")
#######################################################################################################

###################################### Step 3: Prepare the data for the NetCDF document ########################
## To add a new variable with attributes (long_name,units,details) in the NetCDF: modify ONLY the CSV 
## file parameter_netcdf
## To add other attributes such as maximum, minimum, add them manually (see also maximum and minimum for the temperature)
## The dimensions and the details for the NC_Global are inserted manually

### Import parameters from parameter_netcdf.csv
## WARNING the names of the data created in the previous steps have to be the same as the var_name in 
## allparameter (or long_name for NC_GLOBAL).
allparameter <- "data/parameter_netcdf.csv"
allparameter <- read.csv(allparameter, sep=",", dec=".", header=TRUE, fill=TRUE)
if (dataset$metadata$otherTypeLoc != "GPS") { 
  allparameter <- allparameter[-which(allparameter$var_name == "satellites"),]
  allparameter <- allparameter[-which(allparameter$var_name == "residual"),]
}

product <- paste(level,"product",sep = "")

variables <- allparameter[which(allparameter[[product]] == "x"),]
for (i in 1:length(variables)) { variables[,i] <- as.character(variables[,i])}

timevar <- variables[which(variables$units=="seconds since 1970-01-01 00:00:00"),] # select the variables time according to the unit for var.put.nc
dimvar <- dim(variables)[1]
dimtime <- dim(timevar)[1]
#######################################################################################################

###################################### Step 4: Export data in a new NetCDF ############################

### Export the L0 product in a new NetCDF
## Copy and paste the L0 file in another place
newfile <- paste("data/output/",tag_id,"/",level,"/",tag_id,"-",level,".nc", sep="")
file.copy(from = file,to = newfile,overwrite = TRUE)
# Open the new file
ncL1 <- open.nc(newfile, write = TRUE)


### Creation of new NetCDF dimensions, the dimension is defined by the length of the variables
## A NetCDF file may only contain one unlim dimension. We choose the variable with the most observations (time)
dim.def.nc(ncL1, "dive", dimdiv)


### Creation of new NetCDF variables
## Definition of the variables in the NetCDF file with the format 
## var.def.nc(netcdf_file, "variable_name", "value_type","dimension"), such as: 
## - var.def.nc(dataset, "time", "NC_INT","time")
## - var.def.nc(dataset, "lat", "NC_DOUBLE","time")
## - var.def.nc(dataset, "source_loc", "NC_CHAR",c("max_string_32","time")): for the character, 
##   the UNLIM dimension has to be at last in the dimension vector
for ( i in 1:dimvar ) { 
  if (variables$value_type[i]=="NC_CHAR"){
    var.def.nc(ncL1, variables$var_name[i], variables$value_type[i], c(variables$dimCHAR[i],variables$dimL1_2product[i]))
  } else {
    var.def.nc(ncL1, variables$var_name[i], variables$value_type[i], variables$dimL1_2product[i]) 
  }
}


### Put attributes in the variables or in the NC_GLOBAL
## The attributes are either the meta data of the variables or the global meta data of the NetCDF (NC_GLOBAL)
## the format is : att.put.nc(netcdf_file, "variables_name-or-NC_GLOBAL", "meta_data_name", "value_type", data), such as:
## - att.put.nc(dataset, "NC_GLOBAL", "title", "NC_CHAR", title)
## - att.put.nc(dataset, "time", "long_name" , "NC_CHAR", name_time)
## - att.put.nc(dataset, "temp", "_FillValue", "NC_DOUBLE", -99999.9), _FillValue has to be added for the creation of the figures
#
# For variables 
for ( i in 1:dimvar ) { 
  if (variables$standard_name[i] != "") {
    att.put.nc(dataset, variables$var_name[i], "standard_name", "NC_CHAR", variables$standard_name[i])  # add standard name
  }
  if (variables$long_name[i] != "") {
    att.put.nc(ncL1, variables$var_name[i], "long_name", "NC_CHAR", variables$long_name[i])  # add a long name 
  } 
  if (variables$units[i] != "") {
    att.put.nc(ncL1, variables$var_name[i], "units", "NC_CHAR", variables$units[i]) # add the unit for the variables having a unit definition
  }
  if (variables$value_type[i] == "NC_DOUBLE"){
    att.put.nc(ncL1, variables$var_name[i], "_FillValue", "NC_DOUBLE", -99999.9) # -99999.9 to see Michna, P. & Milton Woods. RNetCDF - A Package for Reading and Writing NetCDF Datasets. (2013).
  }
  if (variables$details[i] != "") {
    att.put.nc(ncL1, variables$var_name[i], "details", "NC_CHAR", variables$details[i])  # add details
  }
}
# Other attributes for NC_GLOBAL
detail_1 <- "L1_2 product, analysis of vertical data using diveMove."
att.put.nc(ncL1, "NC_GLOBAL", "detail_1", "NC_CHAR", detail_1 )


### Write the contents of a NetCDF variable.
## format: var.put.nc(netcdf_file, varialable_name, data), such as: var.put.nc(dataset, "lon", fusion$lon)
## the time variable data must be temporarily converted to a UTC referenced date, format of the convertion: dataconvert <- utinvcal.nc(units, data)
## for CHAR the NA must be replaced by ""
## Warning: the var.put.nc will not work if the format of the data is different from the format given in var.def.nc
for (i in 1 : dimvar){
  if (variables$units[i] =="seconds since 1970-01-01 00:00:00"){
      if (variables$dimL1_2product[i] == "dive") { datatime <- tdrXSumm2 
      } else { warning("No data frame found") }
      id.time <- as.numeric(which(colnames(datatime) == variables$var_name[i])) # if problem, check if id is empty
      mytime <- utinvcal.nc("seconds since 1970-01-01 00:00:00", datatime[,id.time]) # conversion of time
      var.put.nc(ncL1, variables$var_name[i], mytime)
  } else  if (variables$value_type[i] == "NC_CHAR"){ # all the variables with character use the dimension phase.no or time 
      if (variables$dimL1_2product[i] == "time") { datachar <- decalib_data 
      } else { warning("No data frame found") }
      id.char <- as.numeric(which(colnames(datachar) == variables$var_name[i])) # if problem, check if id is empty
      mydata <- datachar[,id.char] # select the data
      mydata <- as.character(mydata) #  warning : the "as.character" have to be before ' replace NA by "" '
      mydata[is.na(mydata)] <-"" # replace NA by ""
      var.put.nc(ncL1, variables$var_name[i], mydata) 
  } else {
      if (variables$dimL1_2product[i] == "dive") { datavar <- tdrXSumm2 
      } else { warning("No data frame found") }
      id.var <- as.numeric(which(colnames(datavar) == variables$var_name[i])) # if problem, check if id is empty
      datavar[,id.var] <- as.numeric(datavar[,id.var])
      var.put.nc(ncL1, variables$var_name[i], datavar[,id.var])
  }
}

### View the NetCDF
# print.nc(ncL1)
# var.get.nc(ncL1, "ascdist")

### Close the opened NetCDF file
close.nc(ncL1)
##########################################      END      ##############################################
