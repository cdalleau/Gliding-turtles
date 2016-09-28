########################################################################################################
# Process to generate a NetCDF file with an interpolated (L2) product
# Project: SOCIB Seaturtle
# Author: David March
# Co author : Chloé Dalleau
# Date: 08/08/2016
#
#
# Description
# - Import L1 product
# - Interpolate positions using a state-space model
# - Create derived variables: bearing, distance (using spherical trigonometry) 
#
# References:
# Jonsen et al 2005
########################################################################################################

### Remove the previous data
rm(list=ls()) 

### Import libraries and sources
library(PBSmapping)
library(coda)
library(rjags)
library(bsam)
library(ncdf4)
source( 'ext/lib/ncParse.R')
source("src/utils.R")

###################################### Step 1 : Import data ###########################################

### Set the turtle ID and parameters for the function fitSSM
tag_id <-  00000 # modify the tag_id according to your turtle
ssm_model="DCRWS"
ssm_tstep=0.25
ssm_adapt=30000
ssm_samples=10000
ssm_thin=10
ssm_chains=2
parameters_ssm <- paste("model=",ssm_model, ", tstep=",ssm_tstep, ", adapt=",ssm_adapt, ", samples=",ssm_samples, 
                        ", thin=",ssm_thin, ", chains=",ssm_chains, sep="")

### Meta data import
meta_data <- "data/turtles_metadata.csv"
meta_data <- read.csv(meta_data, sep=",", dec=".", header=TRUE, fill=TRUE)
colnames(meta_data)<-c("argosid", "name", "dateDeployment","refMaxSeaTemp","refMinSeaTemp","refMaxDepth","refMinDepth", "title", "author", "publisher","fileVersion","otherTypeLoc")
meta_data$argosid<-as.character(meta_data$argosid)
meta_data$name<-as.character(meta_data$name)
meta_data$dateDeployment <- as.POSIXct(meta_data$dateDeployment, "%Y-%m-%d %H:%M:%S", tz="GMT")
meta_data$refMaxSeaTemp<-as.numeric(as.character(meta_data$refMaxSeaTemp))
meta_data$refMinSeaTemp<-as.numeric(as.character(meta_data$refMinSeaTemp))
meta_data$refMaxDepth<-as.numeric(as.character(meta_data$refMaxDepth))
meta_data$refMinDepth<-as.numeric(as.character(meta_data$refMinDepth))
meta_data$title<-as.character(meta_data$title)
meta_data$author<-as.character(meta_data$author)
meta_data$publisher<-as.character(meta_data$publisher)
meta_data$fileVersion<-as.character(meta_data$fileVersion)
meta_data$otherTypeLoc<-as.character(meta_data$otherTypeLoc)

### Meta data selection using the turtle ID
select_data <- meta_data[which((meta_data$argosid == tag_id) == "TRUE"),]

### Import and convert the netcdf L1_2 to a list
data_type <- "netcdf"
level <- "L1_1"
file <- paste("data/output/",tag_id,"/",level,"/",tag_id,"-",level,".nc", sep="")
dataset <- ncParse(file)

### Filter data according to QC
lat <- dataset$variables$lat$data[which(dataset$variables$argosfilter$data == "not" & dataset$variables$onland$data == FALSE)]
lon <- dataset$variables$lon$data[which(dataset$variables$argosfilter$data == "not" & dataset$variables$onland$data == FALSE)]
time <-dataset$dimensions$time$data[which(dataset$variables$argosfilter$data == "not" & dataset$variables$onland$data == FALSE)]
time <- parse_date_time(time, c("HMS dbY", "Ymd HMS"), locale=Sys.setlocale("LC_TIME", "English"), tz="GMT")
lc <- dataset$variables$lc$data[which(dataset$variables$argosfilter$data == "not" & dataset$variables$onland$data == FALSE)]
for ( i in 1:length(lc)){ if (lc[i]=="gps") lc[i] <- 3}
argosfilter <- dataset$variables$argosfilter$data[which(dataset$variables$argosfilter$data == "not" & dataset$variables$onland$data == FALSE)]
onland <- dataset$variables$onland$data[which(dataset$variables$argosfilter$data == "not" & dataset$variables$onland$data == FALSE)]
locationdata <- data.frame(time,lat,lon,lc,argosfilter,onland)



### Convert to bsam data format (check for changes in variables)
indata <- data.frame(id = tag_id, date = time, lc = lc, lon = lon, lat = lat)
indata$id<-as.character(indata$id)
indata$lc<-as.character(indata$lc)
indata$lat<-as.numeric(indata$lat)
indata$lon<-as.numeric(indata$lon)

#######################################################################################################

###################################### Step 2: process Of the SSM  ####################################

### Fit DCRWS model for state filtering, regularization and behavioural state estimation
fit = fitSSM(indata, model="DCRWS", tstep=0.25, adapt=30000, samples=10000, thin=10, chains=2)
plotSSM(fit, save.to.pdf=TRUE)
diagSSM(fit, save.to.pdf=TRUE)

### Save SSM results & move pdf plots to output/ssm
level <- "L2"
output.file <- paste("data/output/",tag_id,"/",level,"/",tag_id,".RData",sep="")
save(fit, file = output.file)

file.plot <- paste("ssm", tag_id,".pdf", sep = "")
file.diag <- paste("diag", tag_id, "_params.pdf", sep = "")
filestocopy <- c(file.plot, file.diag)
file.copy(from=filestocopy, to=paste("data/output/",tag_id,"/",level,"/",filestocopy,sep=""))
file.remove(filestocopy)

### Define product level 2
newdata <- fit[[1]]$summary
vars <- c("id", "date", "lon", "lat", "b") # keep mean values
newdata <- newdata[,vars]
colnames(newdata) <- c("tag_id","time","lon","lat","behaviour")
newdata$time <- parse_date_time(newdata$time, c("HMS dbY", "Ymd HMS"), locale=Sys.setlocale("LC_TIME", "English"), tz="GMT") 

# ## Export in a csv file
level <- "L2"
output.file <- paste("data/output/",tag_id,"/",level,"/",tag_id,"-",level,".csv", sep= "")
write.table(newdata, output.file, row.names=FALSE, sep=",", dec=".")
#######################################################################################################

###################################### Step 3: Prepare the data for the NetCDF document ########################
## To add a new variable with attributes (long_name,units,details) in the NetCDF: modify ONLY the CSV 
## file parameter_netcdf
## To add other attributes such as maximum, minimum, add them manually (see also maximum and minimum for the temperature)
## The dimensions and the details for the NC_Global are inserted manually

### Import parameters from parameter_netcdf.csv
## WARNING the names of the data created in the previous steps have to be the same as the var_name in allparameter
allparameter <- "data/parameter_netcdf.csv"
allparameter <- read.csv(allparameter, sep=",", dec=".", header=TRUE, fill=TRUE)
if (dataset$metadata$otherTypeLoc != "GPS") {
  allparameter <- allparameter[-which(allparameter$var_name == "satellites"),]
  allparameter <- allparameter[-which(allparameter$var_name == "residual"),]
}

product <- paste(level,"product",sep = "")

variables <- allparameter[which(allparameter[[product]] == "x"),]
for (i in 1:length(variables)) { variables[,i] <- as.character(variables[,i])}
dimvar <- dim(variables)[1]

glob_att<- "data/nc_global_att.csv"
glob_att <- read.csv(glob_att, sep=",", dec=".", header=TRUE, fill=TRUE)
glob_att <- glob_att[which(glob_att[[product]] == "x"),]
for (i in 1:length(glob_att)) { glob_att[,i] <- as.character(glob_att[,i])}
dimglob <- dim(glob_att)[1]

###################################### Step 4: Export data in a new NetCDF ############################


### Creation of new NetCDF file
data_type <- "netcdf"
filename <- paste("data/output/",tag_id,"/",level,"/",tag_id,"-",level,".nc", sep="")
ncL1 <- create.nc(filename)

### Creation of new NetCDF dimensions, the dimension is defined by the length of the variables
## A NetCDF file may only contain one unlim dimension. We choose the variable with the most observations (time)
dim.def.nc(ncL1, "time", unlim=TRUE)

### Creation of new NetCDF variables
## Definition of the variables in the NetCDF file with the format 
## var.def.nc(netcdf_file, "variable_name", "value_type","dimension"), such as: 
## - var.def.nc(dataset, "time", "NC_INT","time")
## - var.def.nc(dataset, "lat", "NC_DOUBLE","time")
## - var.def.nc(dataset, "source_loc", "NC_CHAR",c("max_string_32","time")): for the character, 
##   the UNLIM dimension has to be at last in the dimension vector
for ( i in 1:dimvar ) { 
  if (variables$value_type[i]=="NC_CHAR"){
    var.def.nc(ncL1, variables$var_name[i], variables$value_type[i], c(variables$dimCHAR[i],variables$dimL2product[i]))
  } else {
    var.def.nc(ncL1, variables$var_name[i], variables$value_type[i], variables$dimL2product[i]) 
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
    att.put.nc(ncL1, variables$var_name[i], "standard name", "NC_CHAR", variables$standard_name[i])  # add standard name
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
## For NC_GLOBAL
## WARNING the names of the data in select_data have to be the same as the att_name in glob_att
for ( i in 1:dimglob ) { 
  if ( length(intersect(colnames(select_data),glob_att$att_name[i])) == 1)  {
    id.glob_att <- which(colnames(select_data) == glob_att$att_name[i]) 
    id.glob_att <- as.numeric(id.glob_att)
    att.put.nc(ncL1, "NC_GLOBAL", glob_att$att_name[i], glob_att$value_type[i],  format(select_data[,id.glob_att]) ) # format is to keep the format of the select_data
  }
}
## Other attributes for NC_GLOBAL
att.put.nc(ncL1, "NC_GLOBAL", "ssm parameters", "NC_CHAR", parameters_ssm )
detail_1 <- "L2 product, interpolate positions using the state-space model from Jonsen et all 2005."
att.put.nc(ncL1, "NC_GLOBAL", "detail_1", "NC_CHAR", detail_1 )
if ( is.null(dataset$metadata$detail_2) == FALSE ){
  detail_2 <- dataset$metadata$detail_2
  att.put.nc(ncL1, "NC_GLOBAL", "detail_2", "NC_CHAR", detail_2 )
}


### Write the contents of a NetCDF variable.
## format: var.put.nc(netcdf_file, varialable_name, data), such as: var.put.nc(dataset, "lon", fusion$lon)
## the time variable data must be temporarily converted to a UTC referenced date, format of the convertion: dataconvert <- utinvcal.nc(units, data)
## for CHAR the NA must be replaced by ""
## Warning: the var.put.nc will not work if the format of the data is diffrent from the format given in var.def.nc
## WARNING the names of the data in the data frames have to be the same as the var_name in variables
for (i in 1 : dimvar){
  if (variables$var_name[i] =="time"){
    mytime <- utinvcal.nc(variables$units[which(variables$var_name == "time")], newdata$time) #conversion of time
    var.put.nc(ncL1, "time", mytime)
  } else  if (variables$value_type[i] == "NC_CHAR"){
    id.char <- as.numeric(which(colnames(newdata) == variables$var_name[i])) # select the id of variables using character in newdata 
    mydata <- newdata[,id.char] # select the data
    mydata <- as.character(mydata) #  warning : the "as.character" have to be before ' remplace NA by "" '
    mydata[is.na(mydata)] <-"" # remplace NA by ""
    var.put.nc(ncL1, variables$var_name[i], mydata) 
  } else {
    id.var <- as.numeric(which(colnames(newdata) == variables$var_name[i])) #select the other id
    var.put.nc(ncL1, variables$var_name[i], newdata[,id.var])
  }
}


### View the NetCDF
# print.nc(ncL1)
# var.get.nc(ncL1, "lat")

### Close the opened NetCDF file
close.nc(ncL1)
##########################################      END      ##############################################
