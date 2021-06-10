require(raster)
require(sp)
require(tidyverse)
require(ggspatial)
require(sf)


################################################################
### Bathymetric Depth
################################################################

load(file="data/depth-data/code-output/combinedDive.Rd")
benthos <- raster::raster("Data/Mapping/Benthic/D4_2018.asc/benthos.tif")

coordinates(combined.dive) <- ~lon+lat

combined.dive$bathymetry <- abs(raster::extract(benthos, combined.dive))
combined.dive$bathymetry <- ifelse(combined.dive$bathymetry < 0, 0, combined.dive$bathymetry)

combined.dive <- as.data.frame(combined.dive) %>%
  drop_na(bathymetry) %>% #some dives will be effectively on land given the reslution of the depth data
  filter(bathymetry>0) %>%
  mutate(prop_benthos = max_depth/bathymetry) %>%
  mutate(seabed_dist = bathymetry-max_depth) %>%
  filter(prop_benthos<1.51)


################################################################
### Tidal State
################################################################

Tides11 <- read_csv("Data/GillsBay2011_2012.csv") 
Tides16 <- read_csv("Data/GillsBay2016_2018.csv") 
Tides <- rbind(Tides11, Tides16)

Tides$Tide <- "NA"                                                  
Tides$Tide[1] <- ifelse(Tides$Level_m[1] < Tides$Level_m[2], "LW", "HW")  
for (i in 2:NROW(Tides)){                                           # loop for remaining rows
  Tides$Tide[i] <- ifelse(Tides$Level_m[i] < Tides$Level_m[i-1], "LW", "HW")  # if height lower than row above = "LW", else "HW"
}


nearest.time <- function(x, y, key = "mid_dt"){y[which.min(abs(difftime(x, y$mid_dt))),key]}

# subset tide data to HW only
HW2016.2018 <- Tides[Tides$Tide == "HW",]

# create posixct DateTime columns
HW2016.2018$mid_dt <- as.POSIXct(paste(HW2016.2018$Date, HW2016.2018$Time), format="%d/%m/%Y %H:%M", origin="1970-01-01", tz = "UTC")

# new column of time of nearest high water
combined.dive$NearestHW <- as.POSIXct(unlist(mapply(function(x) nearest.time(x, HW2016.2018), combined.dive$mid_dt)), origin="1970-01-01", tz = "UTC")

# column of time around nearest high water
combined.dive$TimeAroundHW <- difftime(combined.dive$mid_dt, combined.dive$NearestHW, units = "hours")
combined.dive$TimeAroundHW <- as.numeric(combined.dive$TimeAroundHW)


ggplot()+
  geom_density(data=combined.dive, aes(x=seabed_dist, y=..count..))

################################################################
### Proportion of time at risk depths 
################################################################

diveriskFun <- function(dive_time, dive_depths, bathymetry, shallowrisk=5, deeprisk=25 ) {
 
   
    #temp = df[i,]
    diveSecs = seq(from=0, to=dive_time, by=1) # create sequence of one second intervals for the dive
   
    
    if(is.na(dive_depths[11])){
      diveT = seq(from = 0, to=dive_time, by=dive_time/10) 
      diveD = data.frame(c(0, dive_depths[1:9], 0)) # I have a feeling that these should be bookended by 1.5 rather than 0, but I don;t think it makes a difference for this analysis
      
      # interpolate ldepths onto the one second sequence: outputs is d.reg (Depth)
      d.reg = approx(diveT, diveD, xout=diveSecs)$y
      
      # calculate the depths as distances from the seabed
      d.regSeabed = bathymetry-d.reg
      
      # calculate proportion of dive time in risk zone 
      riskdive = length(d.regSeabed[d.regSeabed>shallowrisk & d.regSeabed<deeprisk ])/dive_time
      
    } else {
      diveT = seq(from = 0, to=dive_time, by=dive_time/12) 
      diveD = data.frame(c(0, dive_depths[1:11], 0)) 
      
      # interpolate ldepths onto the one second sequence: outputs is d.reg (Depth)
      d.reg = approx(diveT, diveD, xout=diveSecs)$y
      
      # calculate the depths as distances from the seabed
      d.regSeabed = bathymetry-d.reg
      
      # calculate proportion of dive time in risk zone 
      riskdive = length(d.regSeabed[d.regSeabed>shallowrisk & d.regSeabed<deeprisk ])/dive_time
      
     
  }
   
  return(riskdive)
    print("done")
}
  

diveriskFun(combined.dive$dive_duration[1], combined.dive[1,7:17],combined.dive$bathymetry[1])

combined.dive$prop_risk <- tibble::tibble(
  dive_time = combined.dive$dive_duration, 
  bathymetry = combined.dive$bathymetry, 
  dive_depths = lapply(1:nrow(combined.dive), function(x) combined.dive[x, 7:17]), 
  shallowrisk = 5, 
  deeprisk = 25
) %>% 
  purrr::pmap(diveriskFun) %>%
  unlist() %>%
  round(digits=3)


################################################################
### Dive angle based on TDR alone
################################################################

## this function will calculate dive angle relative to the horizontal plane based on diving speed alone 
## This is not at all perfect as x dimension is time not spatial so absolute vertical transit (i.e. 90 degrees) is virtualy impossible given time step.
## The angles get more accurate as animal approaches 0 degrees as can assume no depth change is likely an animal on a more or less flat plane if propulsion is assumed


diveanglesFUN <- function(dive_time, dive_depths, bathymetry, shallowrisk=5, deeprisk=25 ) {
  
  ## just to try it on one dive
  #dive_time = combined.dive$dive_duration[2]
  #bathymetry = combined.dive$bathymetry[2] 
  #dive_depths = combined.dive[2, 7:17]
  #shallowrisk = 5
  #deeprisk = 25
  
  
 
  diveSecs = seq(from=0, to=dive_time, by=1) # create sequence of one second intervals for the dive
  
  if(is.na(dive_depths[11])){
    diveT = seq(from = 0, to=dive_time, by=dive_time/10) 
    diveD = data.frame(c(0, dive_depths[1:9], 0)) #
    
    d.reg = approx(diveT, diveD, xout=diveSecs)$y
    
    # calculate the depths as distances from the seabed
    d.regSeabed = bathymetry-d.reg
    
    tdr <- tibble(t=diveSecs, d=d.regSeabed)
    
    # calculate proportion of dive time in risk zone 
    riskdepths = tdr[which(between(tdr$d, shallowrisk, deeprisk)),]
    
    dive_angles <- numeric()
    for(i in 1:nrow(riskdepths)){
      
      dive_angles[i] <- atan2(riskdepths$d[i+1]-riskdepths$d[i], riskdepths$t[i+1]-riskdepths$t[i])*180/pi
      
    }
      
  } else {
    diveT = seq(from = 0, to=dive_time, by=dive_time/12) 
    diveD = data.frame(c(0, dive_depths[1:11], 0)) #
    
    d.reg = approx(diveT, diveD, xout=diveSecs)$y
    
    # calculate the depths as distances from the seabed
    d.regSeabed = bathymetry-d.reg
    
    tdr <- tibble(t=diveSecs, d=d.regSeabed)
    
    # calculate proportion of dive time in risk zone 
    riskdepths = tdr[which(between(tdr$d, shallowrisk, deeprisk)),]
    
    dive_angles <- numeric()
    for(i in 1:nrow(riskdepths)){
      
      dive_angles[i] <- atan2(riskdepths$d[i+1]-riskdepths$d[i], riskdepths$t[i+1]-riskdepths$t[i])*180/pi
    
    
  }}
  
    
  dive_angles_all <- na.omit(dive_angles) 

      des_angle_mean = mean(dive_angles_all[which(dive_angles_all<0)])
      des_angle_sd = sd(dive_angles_all[which(dive_angles_all<0)])
      asc_angle_mean = mean(dive_angles_all[which(dive_angles_all>0)])
      asc_angle_sd = sd(dive_angles_all[which(dive_angles_all>0)])
      
  angle_list <- list(des_angle_mean, des_angle_sd, asc_angle_mean, asc_angle_sd)    
    
  return(angle_list)
  print("done")
  
}

########## try it

#march_to_kill <- tibble::tibble(
#  dive_time = combined.dive$dive_duration[1:10], 
#  bathymetry = combined.dive$bathymetry[1:10], 
#  dive_depths = lapply(1:10, function(x) combined.dive[x, 7:17]), 
 # shallowrisk = 5, 
 # deeprisk = 25
#) %>% 
#  purrr::pmap(diveanglesFUN)


#~~~~~~~~~~~~~~~~~~~~~~~~


dive_angles_output <- tibble::tibble(
  dive_time = combined.dive$dive_duration, 
  bathymetry = combined.dive$bathymetry, 
  dive_depths = lapply(1:nrow(combined.dive), function(x) combined.dive[x, 7:17]), 
  shallowrisk = 5, 
  deeprisk = 25
) %>% 
  purrr::pmap(diveanglesFUN) 

injector <- data.frame(matrix(unlist(dive_angles_output), nrow=length(dive_angles_output), byrow=TRUE),stringsAsFactors=FALSE) # turn list into data frame where every list element is one dive with 4 variables of ascent and descent mean angles and sds

combined.dive <- mutate(combined.dive, des_angle_mean=injector$X1, des_angle_sd=injector$X2, asc_angle_mean=injector$X3, asc_angle_sd=injector$X4); rm(injector)


save(combined.dive, file="data/depth-data/code-output/mod-dive-dat.Rd")
