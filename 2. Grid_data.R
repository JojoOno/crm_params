packages <- c("tidyverse", "fossil", "dplyr", "wesanderson", "sp", "REdaS", "geosphere", "rgdal", "gridExtra", "raster", "ggplot2", "sf", "ggspatial", "data.table")
lapply(packages, require, character.only=TRUE)
projUTM <- "+proj=utm +ellps=WGS84 +datum=WGS84 +zone=30 +north +units=km"
projLATLON <- "+proj=longlat +ellps=WGS84 +datum=WGS84"


### read data
load("proc_data.Rd")

### read shapefile
pent_hr <- st_read("C://users/eo01jo/OneDrive/PhD/Mapping/Scotland/MaximumResolution/ScotLatLon.shp") %>% 
  st_set_crs(projLATLON) %>%
  st_transform(projUTM)



### define coordinates and convert to SpatialPointsDataFrame
points <- seal
coordinates(points) <- ~ lon+lat
proj4string(points) <- projLATLON
points <- spTransform(points, projUTM)
### define SpatialGrid object
bb <- bbox(pent_hr)
cs <- c(3.28084, 3.28084)*6000  # cell size 6km x 6km (for illustration)
# 1 ft = 3.28084 m
cc <- bb[, 1] + (cs/2)  # cell offset
cd <- ceiling(diff(t(bb))/cs)  # number of cells per direction
grd <- GridTopology(cellcentre.offset=cc, cellsize=cs, cells.dim=cd)
grd
# cellcentre.offset 923018 129964
# cellsize           19685  19685
# cells.dim              8      8

sp_grd <- SpatialGridDataFrame(grd,
                               data=data.frame(id=1:prod(cd)),
                               proj4string=CRS(proj4string(shp)))
summary(sp_grd)
# Object of class SpatialGridDataFrame
# Coordinates:
#      min     max
# x 913175 1070655
# y 120122  277602
# Is projected: TRUE
# ...