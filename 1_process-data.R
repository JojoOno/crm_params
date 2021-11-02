load('./Data/GPS-COMBINED-FLOW.Rd')
packages <- c("tidyverse", "fossil", "dplyr", "wesanderson", "sp", "REdaS", "geosphere", "rgdal", "gridExtra", "raster", "ggplot2", "sf", "ggspatial", "data.table")
lapply(packages, require, character.only=TRUE)

gps <- gps_flow_comb %>%
  janitor::clean_names() %>%
  dplyr::select("lat", "lon", "tag", "datetime", "current_speed", "current_dir", "trip_count", "hoflag", "interv") %>%
  mutate(ref=gps_flow_comb$tag)
  
#######################################################################
spDistsPr <- function(pts, longlat=TRUE){
  #	Function to get distances (km) between pairs of coordinates
  #	pts: matrix of coordinates, where each row is a pair of coordinates
  #		- columns 1:2 are LON,LAT of 1st point
  #		- columns 3:4 are LON,LAT of 2nd point
  #	longlat: whether coordinates are in longitude & latitude (default=TRUE)
  if (!is.matrix(pts)) pts <- as.matrix(pts)
  if (nrow(pts)==1){
    return(spDists(pts, longlat=longlat))
  } else {
    dists <- rep(NA, nrow(pts))
    index <- !apply(pts,1,function(x)any(is.na(x)))
    dists[index] <- apply( pts[index,], 1, function(x){spDistsN1(pts=matrix(x[1:2], ncol=2), 
                                                                 pt=matrix(x[3:4], ncol=2), longlat=longlat)} )
    return(dists)
  }
}
#######################################################################

Geo.Hydro.Metrics <- function(gps, refs=FALSE, flow=FALSE){ # gps = your data frame of locations
  
  if(class(gps$datetime)[1]!="POSIXct" ) stop ("Timestamp not in POSIXct")
  
  if(refs==FALSE) {gps$ref <- 1} # includes dummy ref variable so extra if statements not necessary
  ## Calculate time difference in seconds for simpler calculation
  
  
  tmp <- tapply(gps$datetime, gps$ref, function(x){
    c(NA,difftime(tail(x,-1),head(x,-1), tz="GMT", units="secs"))
  })
  
  gps$interv <- unlist(tmp) # create vector of time intervals between each location
  
  gps$interv[is.na(gps$interv)] <- 0
  
  ################ Calculating Distances from one location to another  #######################
  
  
  gps$distance.KM <- gps %>% 
    group_by(ref) %>%
    mutate(lon1=c(head(lon,-1),0), lat1=c(head(lat,-1),0), lon2=c(lon[-1],0), lat2=c(lat[-1],0)) %>%
    subset(select = c(lon1, lat1, lon2, lat2)) %>%
    spDistsPr(longlat = TRUE)
  
  gps$distance.KM <- c(0, gps$distance.KM[-nrow(gps)])  # ensure first location is 0 as cannot have distance between nothing. Just shifts every distance down 1.
  gps$distance.M <- gps$distance.KM*1000
  
  
  ############## Calculate speed over ground (pure point to point (distance / time) calculation in km/s)
  
  gps <- gps %>%
    group_by(ref) %>%
    mutate(speed.MS = distance.M / interv) 
  
  ############## Calculate geocentric orientation
  
  gps <- gps %>% 
    group_by(ref) %>%
    arrange(ref, datetime) %>%
    mutate("LonNext" = c(lon[-1],NA), "LatNext" = c(lat[-1],NA)) %>%
    mutate(swim.dir = earth.bear(lon, lat, LonNext, LatNext) )
  
  
  ############## Calculate current offset
  
  if(flow==TRUE) {
    if(!"current_speed" %in% names(gps)) stop ("Current speed not available")
    if(!"current_dir" %in% names(gps)) stop ("Current direction not available")
    
    Vc<-mutate(gps,
               complex.geo = complex(modulus = speed.MS, argument = swim.dir/360*2*pi),
               complex.flow = complex(modulus = current_speed, argument = (current_dir/360*2*pi)),
               complex.relat = complex.geo-complex.flow, 
               speed.relat = Mod(complex.relat),
               dir.relat = ifelse(Arg(complex.relat)*360/(2*pi)<0,
                                  Arg(complex.relat)*360/(2*pi)+360,
                                  Arg(complex.relat)*360/(2*pi))) 
    
    gps$active.swim.speed <- Vc$speed.relat
    gps$active.swim.dir <- Vc$dir.relat
    
    print("Calculated geospatial and hydrospatial speed and direction")
    
  } else {print("Calculated geospatial speed and direction only")}
  
  
  return(gps)
}


gps <- Geo.Hydro.Metrics(gps_data, flow=TRUE, refs=TRUE)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#


##########################################

at.sea <- subset(gps, gps$hoflag!=1)
seals <- unique(at.sea$tag)
at.sea <- subset(at.sea, at.sea$active.swim.speed<5)

library(data.table);setDT(at.sea)[, bin:=+(distance.KM-shift(distance.KM, fill=distance.KM[1L]) < 2)]
at.sea[, cumul := rleid(bin)]
at.sea[bin == 0, cumul := 0]                  
at.sea[bin == 1, cumul := rleid(cumul)] 
table(at.sea$cumul)

trips <- unique(at.sea$cumul)


for(i in 1:length(trips)){
  
  trip <- subset(at.sea, at.sea$cumul==trips[i])
  
  if(nrow(trip)<2) next
  
  for(j in 1:nrow(trip)-1){
    if(j==0) next # 1:nrow(trip)-1 causes a 0 to be added the beginning of the interger string - this skips the first row as row 0 is meaningless and creates one too many rows to append to the data at the end
    if(j==1){dest <- destPoint(p=data.frame(trip$lon[1], trip$lat[1]), b=trip$active.swim.dir[2], 
                               d=(trip$active.swim.speed[2]*trip$interv[2])) } else {
                                 
                                 dest <- destPoint(p=dest, 
                                                   b=trip$active.swim.dir[j+1], 
                                                   d=(trip$active.swim.speed[j+1]*trip$interv[j+1]))}
    
    if(j==1){ adj.track <- dest} else {adj.track <- rbind(adj.track, dest)}
    
  }
  
  trip$lon_hydr <- c(trip$lon[1], adj.track[,1])
  trip$lat_hydr <- c(trip$lat[1], adj.track[,2])
  
  
  if(i==1){finish <- trip}else{ # i==3 as that's the first trip with >1 locations - think of more transferable way than this 
    finish <- rbind(finish,trip)
  }
  print(i)
}

seal <- finish
seal$secs  <- as.numeric(as.POSIXct(seal$datetime,origin = "1970-01-01"))



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Plotting
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Projection string to use for entire analysis
projUTM <- "+proj=utm +ellps=WGS84 +datum=WGS84 +zone=30 +north +units=km"
projLATLON <- "+proj=longlat +ellps=WGS84 +datum=WGS84"


pent_hr <- st_read("C://users/eo01jo/OneDrive/PhD/Mapping/Scotland/MaximumResolution/ScotLatLon.shp")
pent_lr <- st_read("C://users/eo01jo/OneDrive/PhD/Mapping/Scotland/MaximumResolution/GSSHS_British_Isles poly.shp")


pent_hr <- pent_hr %>% 
  st_set_crs(projLATLON) %>%
  st_transform(projUTM)
pent_lr <- pent_lr %>% 
  st_set_crs(projLATLON) %>%
  st_transform(projUTM)


seal$diffsecs2 <- c(seal$interv[-1], 0)
dest <- destPoint(p=data.frame(seal$lon, seal$lat), b=seal$current_dir, 
                  d=(seal$current_speed*seal$diffsecs2)) 

seal$dest_cur_lon <- dest[,1]
seal$dest_cur_lat <- dest[,2]

coords <- cbind(seal$dest_cur_lon, seal$dest_cur_lat)

sp.points <- SpatialPoints(coords, proj4string = CRS(projLATLON))
sp.points <- spTransform(sp.points, CRS(projUTM))

seal$dest_cur_lonUTM <-  sp.points@coords[,1]
seal$dest_cur_latUTM <-  sp.points@coords[,2]


unique(seal$tag)

tag.sub <- "pv24-580-11"

temp_reg <- filter(seal, tag==tag.sub)

temp_sf_geo <- temp_reg %>%
  st_as_sf(coords=c("lon", "lat"), crs=projLATLON) %>%
  st_transform(crs=projUTM)  

temp_sf_hydro <- temp_reg %>%
  st_as_sf(coords=c("lon_hydr", "lat_hydr"), crs=projLATLON) %>%
  st_transform(crs=projUTM)  

combinetracks <- ggplot() +
  layer_spatial(data=temp_sf_geo, colour=wes_palette("Zissou1")[4], size=0.7) +
  layer_spatial(data=temp_sf_hydro, colour=wes_palette("Zissou1")[1], size=0.5, alpha=0.2) +
  annotation_spatial(pent_lr)+
  theme_bw()

combinetracks

save(seal, file="proc_data.Rd")
