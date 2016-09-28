########################################################################################################
# Process to generate a NetCDF file with an analysis of the horizontal data (L1.1 product)
# Project: SOCIB Seaturtle
# Author: Chloé Dalleau
# Date: 08/08/2016
# Based on the script of David March
#
# Description:
# This script is the second step of a global processing about the analysis of the turtles' trajectory.
# The following steps create a NetCDF file with the analysis of the horizontal data :
# - Import Netcdf L0 product
# - Filter location data according to the position on land
# - Filter location data using the algorithm described in Freitas et al 2008
# - Interpolation of location and battery voltage
# - Calcul of the speed and the cumulative distance using either all the locations or only the location
#   filtered by Freitas et al 2008 algorithm
# - Tests if the temperature and the depth are in the defined range according to the defined metadata
# - Export in a NetCDF file
#
# References:
# Freitas, C., Lydersen, C., Ims, R.A., Fedak, M.A. and Kovacs, K.M. (2008) A simple new algorithm
# to filter marine mammal Argos locations Marine Mammal Science 24:315-325.
########################################################################################################

### Remove the previous data
rm(list=ls()) 

### Import libraries and sources
source( 'ext/lib/ncParse.R')
source("src/utils.R")
library(ncdf4)
library(RNetCDF)
library(argosfilter)
library(stats4)
library(sp)
library(maptools)
library(akima)
library(geosphere)
library(dplyr)
library(diveMove)
library(lubridate)
library(zoo)

###################################### Step 1: Import data ###########################################

### Set the turtle ID and parameters
tag_id <-  00000 # modify the tag_id according to your turtle
ref_max_gps <- 100 # a value over ref_max_gps may indicate an erroneous location
ref_max_sea_temp <- 30 # reference for the maximum of sea water temperature
ref_min_sea_temp <- 10 # reference for the minimum of sea water temperature
ref_max_depth <- 100 # reference for the maximum of depth
ref_min_depth <- 0 # reference for the minimum of depth
sda_vmax <- 1.39  # value of the maximum of velocity using in sdafilter in m/s
sda_ang <- c(15, 25) # value of the angle using in sdafilter, no spikes are removed if ang=-1
sda_distlim <- c(2500, 5000) # value of the limite distance using in sdafilter, no spikes are removed if ang=-1

## Define parameters data frame
sda_parameter <- paste("vmax: ", sda_vmax, ", ang: c(",sda_ang[1],",",sda_ang[2],"), distlim: c(" , sda_distlim[1] , ",", sda_distlim[2], ")", sep="")
ref_data <- matrix( data = c("lc","ref_max_gps", "NC_DOUBLE", ref_max_gps,"temp_range","ref_max_sea_temp","NC_DOUBLE",ref_max_sea_temp,
                             "temp_range","ref_min_sea_temp","NC_DOUBLE",ref_min_sea_temp,"depth_range","ref_max_depth","NC_DOUBLE",ref_max_depth,
                             "depth_range","ref_min_depth","NC_DOUBLE",ref_min_depth),
                    nrow = 5, ncol = 4, byrow=T, dimnames = NULL)
ref_data <- as.data.frame.matrix(ref_data)
colnames(ref_data) <- c("var_name","att_name","value_type","data")
for (i in 1:3) ref_data[,i] <- as.character(ref_data[,i])
ref_data$data <- as.numeric(as.character(ref_data$data))



### Import and convert the NetCDF L0 to a list
data_type <- "netcdf"
level <- "L0"
file <- paste("data/output/",tag_id,"/",level,"/",tag_id,"-",level,".nc", sep="")
dataset <- ncParse(file)

#######################################################################################################

###################################### Step 2: Analysis of data #######################################

### Filter location using Argosfilter
## Preparation of data for sdafilter (the observations are selected according to the source)
## Warning: don't work with NA entries in the location csv file
lat <- dataset$variables$lat$data[is.na(dataset$variables$lat$data) == FALSE] 
lon <- dataset$variables$lon$data[is.na(dataset$variables$lat$data) == FALSE]
time <- dataset$dimensions$time$data[is.na(dataset$variables$lat$data) == FALSE]
lc <- dataset$variables$lc$data[is.na(dataset$variables$lat$data) == FALSE]
if (dataset$metadata$otherTypeLoc == "GPS") residual <- dataset$variables$residual$data[is.na(dataset$variables$lat$data) == FALSE]
csvfile <- dataset$variables$source_loc$data[is.na(dataset$variables$lat$data) == FALSE]

## Put "-9" in lc for the bad values from GPS
if (dataset$metadata$otherTypeLoc == "GPS") {
  location_data <- data.frame(time,lc,residual,lat,lon,csvfile)
  for (i in 1:dim(location_data)[1]){
    if (is.na(location_data$residual[i])==F && location_data$residual[i]> ref_max_gps) location_data$lc[i] <- -9
  }
} else {
  location_data <- data.frame(time,lc,lat,lon,csvfile)
}


### Positions on land
landmask <- readShapePoly("ext/data/gshhs/landmask_med", proj4string=CRS("+proj=longlat +ellps=WGS84"))
location_data$onland <- point.on.land(lat = lat, lon = lon, landmask) # don't work if NA value

### Filter Argos locations
location_data$argosfilter <- sdafilter(lat = lat, lon = lon, dtime = time, lc = lc,
                              vmax = sda_vmax, # in m/s
                              ang = sda_ang, # No spikes are removed if ang=-1
                              distlim = sda_distlim)

### Insert the data from Argosfilter and onland in the data
dim_time <- length(dataset$dimensions$time$data)
j=0
for (i in 1:dim_time){
  if ( is.na(dataset$variables$lat$data[i]) == FALSE){
    j=j+1
    dataset$variables$argosfilter[i] <- location_data$argosfilter[j]  # WARNING: order of time should be the same between netCDF and filters!
    dataset$variables$onland[i] <- location_data$onland[j]  # WARNING: order of time should be the same between netCDF and filters!
  } else {
    dataset$variables$argosfilter[i] <- ""
    dataset$variables$onland[i] <- ""
  }
}
## To view the data
# viz <- data.frame(dataset$dimensions$time$data,dataset$variables$lat$data,dataset$variables$lon$data,dataset$variables$argosfilter,dataset$variables$onland)
# colnames(viz) <- c("time","lat","lon","argosfilter","onland")


### Interpolation of latitude and longitude using the filtered locations
## Selection of the filtered data 
lat_sda <- location_data$lat[which((location_data$argosfilter == "not") == "TRUE")]
lon_sda <- location_data$lon[which((location_data$argosfilter == "not") == "TRUE")]
time_sda <- location_data$time[which((location_data$argosfilter == "not") == "TRUE")]
## Interpolation of data
## Note : the function aspline from akima can produce NaN values, the function na.spline from zoo gives bad interpolations
## Note 2 : the function na.approx from zoo, performs an interpolation between the first location and the last location with 
##          values (non NA). So the final lat and lon can have NA for the times before the first location data and for the  
##          times after the last location data. No extrapolation is performed.
loc <- data.frame(dataset$variables$lat$data,dataset$variables$lon$data)
colnames(loc)<-c("lat","lon")
object <- zoo(loc, order.by = ,as.POSIXct(dataset$dimensions$time$data, "%Y-%m-%d %H:%M:%S", tz="GMT"))
interp_loc_zoo <- na.approx(object = object, na.rm = F)
lat_zoo <- round(coredata(interp_loc_zoo[,"lat"]), 4) #coredata : extract only the data from time serie
lon_zoo <- round(coredata(interp_loc_zoo[,"lon"]), 4)
## To view the interpolated location
# viz <- data.frame(dataset$dimensions$time$data,lat_zoo, lon_zoo)
# colnames(viz) <- c("time" ,"lat","lon")


### Interpolation of voltage battery: function aspline (akima) with "original" method, WARNING erroneous data can be used
## Selection of the data
batt <- dataset$variables$batt$data[which((is.na(dataset$variables$batt$data)) == "FALSE")]
time_batt <- dataset$dimensions$time$data[which((is.na(dataset$variables$batt$data)) == "FALSE")]
## Interpolation
interp_batt <- aspline(x=time_batt, y=batt,xout = dataset$dimensions$time$data , method="original")$y # aspline can create Nan values


### Save initials lat, lon and batt (without interpolation) to check the results. This data will not be exported in the NetCDF file
dataset$variables$lat$dataini <- dataset$variables$lat$data 
dataset$variables$lon$dataini <- dataset$variables$lon$data 
dataset$variables$batt$dataini <- dataset$variables$batt$data 

### Insert the interpolation of lat, lon and batt if NA in the data
result_lat <- insertinterp(interp_data=lat_zoo,data=dataset$variables$lat$data,source=dataset$variables$source_loc$data)
dataset$variables$lat$data <- result_lat$data
dataset$variables$lat$source_lat <- result_lat$source
result_lon <- insertinterp(interp_data=lon_zoo,data=dataset$variables$lon$data,source=dataset$variables$source_loc$data)
dataset$variables$lon$data <- result_lon$data
dataset$variables$lon$source_lon <- result_lon$source
result_batt <- insertinterp(interp_data=interp_batt,data=dataset$variables$batt$data,source=dataset$variables$source_status$data)
dataset$variables$batt$data <- result_batt$data
dataset$variables$batt$source_batt <- result_batt$source

### Calculate the cumulative distance and speed of all data or if Argosfilter==not 
result_dist_argos <- distandSpeed(dtime = dataset$dimensions$time$data,lat=dataset$variables$lat$data, lon=dataset$variables$lon$data, dataremoved=dataset$variables$argosfilter, crit="not")
dataset$variables$argoscumdist <- result_dist_argos$cumdist # km
dataset$variables$argosspeed <-  round(result_dist_argos$speed,6)  # m/s 
result_dist_all <- distandSpeed(dtime = dataset$dimensions$time$data,lat=dataset$variables$lat$data, lon=dataset$variables$lon$data, dataremoved=NA, crit=NA)
dataset$variables$cumdist <- result_dist_all$cumdist # km
dataset$variables$speed <- round(result_dist_all$speed,6) # m/s 


### Test if the temperature and the depth are in the defined ranges as defined in the metadata
dataset$variables$temp$quality <- testrange(data=dataset$variables$temp$data, min=ref_min_sea_temp ,max=ref_max_sea_temp)
dataset$variables$depth$quality <- testrange(data=dataset$variables$depth$data, min=ref_min_depth ,max=ref_max_depth)

### Data export in a data frame
if (dataset$metadata$otherTypeLoc == "GPS"){
  alldata <- data.frame(parse_date_time(dataset$dimensions$time$data, c("HMS dbY", "Ymd HMS"), locale=Sys.setlocale("LC_TIME", "English"), tz="GMT"),
              dataset$variables$lc$data, dataset$variables$satellites$data,dataset$variables$residual$data, dataset$variables$lat$dataini, dataset$variables$lat$data,
              dataset$variables$lon$dataini, dataset$variables$lon$data, dataset$variables$lon$source_lon, dataset$variables$batt$dataini,
              dataset$variables$batt$data, dataset$variables$batt$source_batt, dataset$variables$temp$data, dataset$variables$temp$quality,
              dataset$variables$depth$data, dataset$variables$depth$quality, dataset$variables$argosfilter, dataset$variables$onland,
              dataset$variables$cumdist, dataset$variables$argoscumdist, dataset$variables$speed, dataset$variables$argosspeed)
  colnames(alldata) <- c("time","lc","satellites","residual","latini","lat","lonini","lon","source_loc","battini","batt","source_status","temp",
                         "temp_range","depth","depth_range","argosfilter","onland","cumdist","argoscumdist","speed","argosspeed")
} else {
  alldata <- data.frame(parse_date_time(dataset$dimensions$time$data, c("HMS dbY", "Ymd HMS"), locale=Sys.setlocale("LC_TIME", "English"), tz="GMT"),
                        dataset$variables$lc$data, dataset$variables$lat$dataini, dataset$variables$lat$data, dataset$variables$lon$dataini,
                        dataset$variables$lon$data, dataset$variables$lon$source_lon, dataset$variables$batt$dataini, dataset$variables$batt$data,
                        dataset$variables$batt$source_batt, dataset$variables$temp$data, dataset$variables$temp$quality, dataset$variables$depth$data,
                        dataset$variables$depth$quality, dataset$variables$argosfilter, dataset$variables$onland,dataset$variables$cumdist, 
                        dataset$variables$argoscumdist,dataset$variables$speed,dataset$variables$argosspeed)
  colnames(alldata) <- c("time","lc","latini","lat","lonini","lon","source_loc","battini","batt","source_status","temp","temp_range","depth",
                         "depth_range","argosfilter","onland","cumdist","argoscumdist","speed","argosspeed")
  
}


### Export in a csv file
level <- "L1_1"
output.file <- paste("data/output/",tag_id,"/",level,"/",tag_id,"-",level,".csv", sep= "")
if (dataset$metadata$otherTypeLoc == "GPS"){
  data <- alldata[c("time","lc","satellites","residual","lat","lon","source_loc","batt","source_status","temp",
                  "temp_range","depth","depth_range","argosfilter","onland","cumdist","argoscumdist","speed","argosspeed")]
} else {
  data <- alldata[c("time","lc","lat","lon","source_loc","batt","source_status","temp",
                    "temp_range","depth","depth_range","argosfilter","onland","cumdist","argoscumdist","speed","argosspeed")]
} 
write.table(alldata, output.file, row.names=FALSE, sep=",", dec=".")

#######################################################################################################

###################################### Step 3: Prepare the data for the NetCDF document ########################
## To add a new variable with attributes (long_name,units,details) in the NetCDF:  modify the previous step and the CSV 
## file parameter_netcdf
## To add other attributes such as maximum, minimum, add them manually (see also maximum and minimum for the temperature)
## The dimensions and the details for the NC_Global must also be inserted manually

### Import parameters from parameter_netcdf.csv
## WARNING the name of the data created in the previous step have to be the same as the var_name in allparameter
allparameter <- "data/parameter_netcdf.csv"
allparameter <- read.csv(allparameter, sep=",", dec=".", header=TRUE, fill=TRUE)
if (dataset$metadata$otherTypeLoc != "GPS") { 
  allparameter <- allparameter[-which(allparameter$var_name == "satellites"),]
  allparameter <- allparameter[-which(allparameter$var_name == "residual"),]
}

product <- paste(level,"product",sep = "")

variables <- allparameter[which(allparameter[[product]] == "x"),]
for (i in 1:length(variables)) { variables[,i] <- as.character(variables[,i])}

newvar <- variables[-which(variables$L0product=="x"),] # select the new variables for var.def.nc
dimvar <- dim(variables)[1]
dimnew <- dim(newvar)[1]
#######################################################################################################

###################################### Step 4: Export data in a new NetCDF file #######################

### Export the L0 product in a new NetCDF file
## Copy and paste the L0 file in an other place
newfile <- paste("data/output/",tag_id,"/",level,"/",tag_id,"-",level,".nc", sep="")
file.copy(from = file,to = newfile,overwrite = TRUE)
# Open the new file
ncL1 <- open.nc(newfile, write = TRUE)


### Creation of new NetCDF variables
## Definition of the variables in the NetCDF file with the format 
## var.def.nc(netcdf_file, "variable_name", "value_type","dimension"), such as: 
## - var.def.nc(dataset, "time", "NC_INT","time")
## - var.def.nc(dataset, "lat", "NC_DOUBLE","time")
## - var.def.nc(dataset, "source_loc", "NC_CHAR",c("max_string_32","time")): for the character, 
##   the UNLIM dimension has to be at last in the dimension vector
for ( i in 1:dimnew ) { 
  if (newvar$value_type[i]=="NC_CHAR"){
    var.def.nc(ncL1, newvar$var_name[i], newvar$value_type[i], c(newvar$dimCHAR[i],newvar$dimL1_1product[i]))
  } else {
    var.def.nc(ncL1, newvar$var_name[i], newvar$value_type[i], newvar$dimL1_1product[i]) 
  }
}



### Put attributes in the variables or in the NC_GLOBAL
## The attributes are either the meta data of the variables or the global meta data of the NetCDF (NC_GLOBAL)
## the format is : att.put.nc(netcdf_file, "variables_name-or-NC_GLOBAL", "meta_data_name", "value_type", data), such as:
## - att.put.nc(dataset, "NC_GLOBAL", "title", "NC_CHAR", title)
## - att.put.nc(dataset, "time", "long_name" , "NC_CHAR", name_time)
## - att.put.nc(dataset, "temp", "_FillValue", "NC_DOUBLE", -99999.9), _FillValue has to be added for the creation of the figures
#
## For new variables
for ( i in 1:dimnew ) { 
  if (newvar$standard_name[i] != "") {
    att.put.nc(dataset, newvar$var_name[i], "standard name", "NC_CHAR", newvar$standard_name[i])  # add standard name
  }
  if (newvar$long_name[i] != "") {
    att.put.nc(ncL1, newvar$var_name[i], "long_name", "NC_CHAR", newvar$long_name[i])  # add a long name 
  } 
  if (newvar$units[i] != "") {
    att.put.nc(ncL1, newvar$var_name[i], "units", "NC_CHAR", newvar$units[i]) # add the unit for the variables having an unit 
  }
  if (newvar$value_type[i] == "NC_DOUBLE"){
    att.put.nc(ncL1, newvar$var_name[i], "_FillValue", "NC_DOUBLE", -99999.9) # -99999.9 to see Michna, P. & Milton Woods. RNetCDF - A Package for Reading and Writing NetCDF Datasets. (2013).
  }
  if (newvar$details[i] != "") {
    att.put.nc(ncL1, newvar$var_name[i], "details", "NC_CHAR", newvar$details[i])  # add details
  }
}
## Other attributes for the variables
# WARNING the names in column "var_name" from ref_data have to be the same as the variables defined previously
for ( i in 1 : dim(ref_data)[1]) {
  att.put.nc(ncL1, ref_data$var_name[i], ref_data$att_name[i], ref_data$value_type[i], ref_data$data[i])
}
att.put.nc(ncL1, "argosfilter", "parameters", "NC_CHAR", sda_parameter)
## Other attributes for NC_GLOBAL
detail_1 <- "L1_1 product (analysis of horizontal data) : (1) filtration of location ; (2) interpolation of lat,lon and batt ; (3) Position on land ; (4) calcul of cumulative distance and the speed, (5) quality of temperature and depth."
att.put.nc(ncL1, "NC_GLOBAL", "detail_1", "NC_CHAR", detail_1 )



### Write the contents of a NetCDF variable.
## format: var.put.nc(netcdf_file, varialable_name, data), such as: var.put.nc(dataset, "lon", fusion$lon)
## the time variable data must be temporarily converted to a UTC referenced date, format of the convertion: dataconvert <- utinvcal.nc(units, data)
## for CHAR the NA must be replaced by ""
## Warning: the var.put.nc will not work if the format of the data is diffrent from the format given in var.def.nc
for (i in 1 :dimvar ){ 
  if (variables$value_type[i] == "NC_CHAR"){
    id.char <- as.numeric(which(colnames(data) == variables$var_name[i])) # select the id of variables using character in fusion 
    mydata <- data[,id.char] # select the data
    mydata <- as.character(mydata) #  warning : the "as.character" have to be before ' remplace NA by "" '
    mydata[is.na(mydata)] <-"" # remplace NA by ""
    var.put.nc(ncL1, variables$var_name[i], mydata) 
  } else {
    id.var <- as.numeric(which(colnames(data) == variables$var_name[i])) #select the other id
    var.put.nc(ncL1, variables$var_name[i], data[,id.var])
  }
}


### View the NetCDF
# print.nc(ncL1)
# var.get.nc(ncL1, "depth_range")

### Close the opened NetCDF file
close.nc(ncL1)
##########################################      END      ##############################################
