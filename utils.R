########################################################################################################
# Unpackaged set of custom functions developed for the project
# Project: SOCIB Seaturtle, OASIS
# Author: David March & Chloe Dalleau
########################################################################################################

#############################################
# LIST OF FUNCTIONS
#
# ncdf2df           netCDF to data frame
# point.on.land     Check if location overlap with landmask
# insertinterp      Insert interpolated data keeping the raw data in the final data
# distandSpeed      Calculate the distance and speed with latitude and longitude (with argosfilter=not)
# testrange         Make a classification ("O.K." , "out range") of the data according to a range of reference data.


#----------------------------------------------------------------------------------------------
# ncdf2df     netCDF to data frame
#----------------------------------------------------------------------------------------------
ncdf2df <- function(nc){
  #
  # Description:
  # This function converts a netCDF into a data.frame. Current version only works for trajectory feature type
  #
  # Arguments:
  # nc = netcdf file
  #
  # Author:
  # David March, SOCIB
  #
  # Example: 
  # nc <- "data/locations/L1/dep0001_turtle-manda_scb-ttrk009_L1_2009-08-07.nc"
  # nc.df <- ncdf2df(nc)
  
  # Load dependencies
  require(ncdf4)
  source("ext/lib/ncParse.R") 
  
  # Parse netCDF
  mync <- ncParse(nc)
  
  # check that netCDF is a trajectory
  if (mync$metadata$featureType != "trajectory") stop ("feature type is not a trajectory")
  
  # Extract time dimension and variables
  time <- mync$dimensions$time$data
  vars <- data.frame(sapply(mync$variables, '[[', "data"))  # extract elements from sublist
  
  # Create data.frame
  df <- data.frame(time, vars)
  return(df)
}
#----------------------------------------------------------------------------------------------


#----------------------------------------------------------------------------------------------
# point.on.land    Check if location overlap with landmask
#----------------------------------------------------------------------------------------------
point.on.land <- function(lon, lat, landmask){
  
  require(sp)
  
  xy <- cbind(lon,lat)
  pts <- SpatialPoints(xy, proj4string=CRS("+proj=longlat +ellps=WGS84"))
  ov <- over(pts, landmask)
  onland <- !is.na(ov$id)
  return(onland)
}
#----------------------------------------------------------------------------------------------

#----------------------------------------------------------------------------------------------

#----------------------------------------------------------------------------------------------
# insertinterp         insert interpolation of data 
#----------------------------------------------------------------------------------------------
insertinterp <- function(interp_data,data,source){
  #
  # Description:
  # insert interpolated data  keeping the raw data in the final data
  #
  # Arguments:
  # interp_data = interpolated data
  # data = data where the interpolate data have to add
  # source = the initial source for each observations
  #
  # Values:
  # data frame with a list of values with both interpolated and raw data and a list of the source of this data.
  #
  # Details:
  # This function adds interpolated data where there are NA and the source of the data (interpolation or the initial source)
  # 
  #
  # Author:
  # Chloe Dalleau
  #
  # Example: 
  # result_lat <- insertinterp(interp_data=lat_interp,data=dataset$variables$lat$data,source=dataset$variables$source_loc$data)
  
  dim_data <- length(data)
  result <- c()
  for (i in 1:dim_data){
    if (is.na(data[i])){
      result$data[i] <- interp_data[i]
      result$source[i] <- "interpolation"
    } else {
      result$data[i] <- data[i]
      result$source[i] <- source[i]
    }
  }
  return(result)
}
#----------------------------------------------------------------------------------------------

#----------------------------------------------------------------------------------------------
# distandSpeed          Calculate the distance and speed with latitude and longitude (with argosfilter=not or argosfilter=end_location)
#----------------------------------------------------------------------------------------------
distandSpeed <- function(dtime,lat,lon,dataremoved,crit){
  #
  # Description:
  # Calcul of the cumulative distance and the Speed for locations according to a parameter
  #
  # Arguments:
  # dtime = date with day and time with Posixct format
  # lat = latitude in degree (can have NA values)
  # lon = longitude in degree (can have NA values)
  # dataremoved = qualitative table indicating what data to use
  # crit = criterion in dataremoved indicating what data to use
  #
  # Values:
  # dataframe with 4 columns, the three first column are the distance (dist), the time difference (timediff) and the speed (speed) between
  # two good locations. The last column is the cumulative distance.
  #
  # Details:
  # This function (1) prepares the data, (2) uses the package diveMove to caculate the distance and the speed
  # for all locations and (3) calculates the cumulative distance using the distance from the diveMove package.
  # The data can have NA values, for these observations the speed and the cumulative distance are not calculated
  # 
  #
  # Author:
  # Chloe Dalleau
  #
  # Example: 
  # distandSpeed(dtime = dataset$dimensions$time$data,lat=dataset$variables$lat$data, lon=dataset$variables$lon$data, dataremoved=dataset$variables$argosfilter, crit="not")
  
  require(diveMove)
  
  #preparation of data
  loc <- data.frame(as.POSIXct(dtime, "%Y-%m-%d %H:%M:%S", tz="GMT"),lon,lat) #dataframe with the date, lon and lat
  colnames(loc) <- c("dtime","lon","lat")
  id.na <- c()
  id.obs <- c()
  for ( i in 1:dim(loc)[1]){
    if (is.na(loc$lat)[i] == T){
      id.na <- c(id.na,i)
    } else { id.obs <- c(id.obs,i)}
  }
  size <- length(lat)
  result <- data.frame(NA,NA,NA) 
  result[id.obs[1],] <- c(0,0,0) # put 0 for the first location 
  colnames(result) <- c("dist","timedif","speed") # distance (km), time difference (s), and speed (m/s). according to the package diveMove
  cumdist <- data.frame(NA) # data frame for the cumulative distance
  cumdist[id.obs[1],] <- 0
  colnames(cumdist) <- c("cumulativeDist")
  
  
   if ( length(dataremoved) == 1) { # if the function uses all the observations 
       for (i in 2:length(id.obs)){ # calcul of the speed and the cumulative distance for the observations without missing value.
         result[id.obs[i],] <- distSpeed(loc[id.obs[i-1],],loc[id.obs[i],], method = "VincentyEllipsoid") # A matrix with three columns: distance (km), time difference (s), and speed (m/s).
         cumdist[id.obs[i],] <- cumdist[id.obs[i-1],] + result$dist[id.obs[i]]
      }
   } else { # if the function uses only certain observations
     dataremoved[id.obs[1]] <- "start"
     for (i in 2:length(id.obs)){
         a=i-1
         if (dataremoved[id.obs[i]]== crit) {
           while (dataremoved[id.obs[a]] != crit && dataremoved[id.obs[a]] != "start") {a=a-1} # if for the previous location Argosfilter=removed, the last good argosfilter is took
           result[id.obs[i],] <- distSpeed(loc[id.obs[a],],loc[id.obs[i],], method = "VincentyEllipsoid") # A matrix with three columns: distance (km), time difference (s), and speed (m/s).
           cumdist[id.obs[i],] <- cumdist[id.obs[a],] + result$dist[id.obs[i]]
         } else {
           result[id.obs[i],] <- NA # if Argosfilter=removed or "", NA is put
           cumdist[id.obs[i],] <- NA # if Argosfilter=removed or "", NA is put
         }
       }
    }
  
  
  result$cumdist <- cumdist$cumulativeDist
  for (i in 1:length(id.na)) result[id.na[i],]<-NA # put NA for all the observations with location = NA 
  
  return(result)
}
#----------------------------------------------------------------------------------------------

#----------------------------------------------------------------------------------------------
# testrange          check if the observation are between two numbers
#----------------------------------------------------------------------------------------------
testrange <- function(data,min,max){
  #
  # Description:
  # make a classification ("O.K." , "out range") of the data according to a range of reference data.
  #
  # Arguments:
  # data = data to test
  # min = minimum value of reference
  # max  maximum value of reference
  #
  # Values:
  # list of qualitative value : "O.K" or "outrange"
  #
  # Details:
  # This function adds a qualitative value where there aren't NA
  # 
  #
  # Author:
  # Chloé Dalleau
  #
  # Example: 
  # dataset$variables$depth$quality <- testrange(data=dataset$variables$depth$data, min=select_data$refMinDepth ,max=select_data$refMaxDepth)
  
  dim <- length(data)
  quality <- c()
  for ( i in 1:dim){
    
    if (is.na(data[i]) == "FALSE") {
        if (data[i] >= min && data[i] <= max) {
          quality[i] <- "O.K." 
        } else {
          quality[i] <- "outrange"
        }
    } else {
      quality[i] <- ""
    }
  }
  return(quality)
}

#----------------------------------------------------------------------------------------------
