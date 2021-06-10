require(DBI)
require(odbc)
require(dplyr)
require(dbplyr)
require(sf)
require(tidyverse)
require(data.table)
require(ggspatial)
require(ggplot2)

### Map to the database using odbc
db <- dbConnect(odbc(), .connection_string = "Driver={Microsoft Access Driver (*.mdb, *.accdb)};
                        DBQ=D:/OneDrive/RCode/Git Clones - UHI Cpu/CRMParam/Data/pv24.accdb")

######################
### Load Dive data ###
######################

dive <- tbl(db, "dive") %>%
  dplyr::select("ref", "DE_DATE", "SURF_DUR", "DIVE_DUR", "MAX_DEP", "DEPTH_STR","lat", "lon", "start_lat", "start_lon", "N_DEPTHS") %>%
  collect() %>%
  mutate(ds_date = DE_DATE-DIVE_DUR) %>% # for some reason dive start time throws up an error due to formatting (which looks fine in the database) so just manually re-create it using dive end timestammp and dive duration 
  drop_na(lon) %>%
  drop_na(lat) %>%
  mutate(D_DATE = as.POSIXct((as.numeric(DE_DATE)+as.numeric(ds_date))/2, origin='1970-01-01 00:00:00', tz="UTC"))

temp <- do.call("rbind", strsplit(dive$DEPTH_STR, ",")) # make TDR records separate columns for ease of plitting and manipulation 
temp <- data.frame(apply(temp, 2, as.numeric))
names(temp) <- paste("depth", 1:9, sep="")
dive <- cbind(dive, temp); rm(temp)

## check with plots
dive_num <- which(dive$MAX_DEP>110)

plot(as.numeric(dive[dive_num[2],14:22])~seq(dive$ds_date[dive_num[2]], dive$DE_DATE[dive_num[2]], length.out=9), 
     type="l", ylab="Depth (m)", xlab="Time", ylim=rev(range(as.numeric(dive[dive_num[2],14:22]))))

######################
#### Load GPS Data ###
######################

battlecross <- function(df){
  df$days_diff <- c(0, as.numeric(diff(df$D_DATE)/60/60/24))
  df <- slice(df, c(which(df$days_diff<1)))
  df <- tibble(df)
  return(df)
  } #function to remove locations if there is >1day between them


gps <- tbl(db, "gps") %>%
  dplyr::select(ref, D_DATE, NSATS_TRANSMITTED, LAT, LON, RESIDUAL, EST_SPEED) %>%
  collect() %>%
  filter(RESIDUAL>0 & RESIDUAL<25 & NSATS_TRANSMITTED >5) %>% #can probably get away with not using sattelite threshold but benthic data is high res so useful to ensure as high accuracy as possible (Russell et al. (2011))
  mutate(tag=ref) %>% #here as ref will be superceded in the function as will create a nested list and then revert back to a tibble
  group_by(ref) %>%
  group_map(~battlecross(.x)) %>%
  rbindlist(., fill=TRUE, idcol="ref") %>%
  dplyr::select(-ref) %>%
  rename(ref=tag) 


gps_sf <- st_as_sf(gps, coords=c("LON", "LAT"), crs=4326) %>%
  st_transform(crs=32630)

ggplot()+
 geom_sf(data=gps_sf, aes(colour=ref))


###################################################################
#### Clean dives by comparing times to nearest cleaned GPS fix ####
###################################################################

z <- lapply(intersect(dive$ref,gps$ref),function(id) { #function to match the mid-point of the dive times to the closest cleaned GPS location
  d1 <- subset(dive,ref==id)
  d2 <- subset(gps,ref==id)
  
  d1$indices <- sapply(d1$D_DATE,function(d) which.min(abs(d2$D_DATE - d)))
  d2$indices <- 1:nrow(d2)
  
  merge(d1,d2,by=c('ref','indices'))
})

dive <- do.call(rbind,z) %>%
  mutate(time_diff = as.numeric(abs(D_DATE.y - D_DATE.x))) %>%
  filter(time_diff<1050) # difftime is in seconds so <1800 filters out any dives greater than 30 minutes from the nearest GPS fix

dive_sf <- st_as_sf(dive, coords=c("LON", "LAT"), crs=4326) %>%
  st_transform(crs=32630)

ggplot()+
  geom_sf(data=dive_sf, aes(colour=ref))

save(dive, file="Data/Depth Data/Code Output/GSMDive_cleaned.Rd")


