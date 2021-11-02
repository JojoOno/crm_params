packages <- c("tidyverse", "move", "fossil", "dplyr", "wesanderson", "sp", "REdaS", "ggsubplot", "geosphere", "rgdal", "gridExtra", "raster", "ggplot2", "sf", "ggspatial", "data.table")
lapply(packages, require, character.only=TRUE)

### read data
load("data/proc_data.Rd")

### read shapefile
pent_hr <- st_read("data/mapping/scotland/maximum-resolution/ScotLatLon.shp") %>% 
  st_set_crs(4326) %>%
  st_transform(32630)

pent_lr <- st_read("data/mapping/scotland/maximum-resolution/GSSHS_British_Isles poly.shp") %>% 
  st_set_crs(4326) %>%
  st_transform(32630)

lease_site <- st_read("data/mapping/wave-and-tidal-proposed-sites/TCE_Lease_Tide_20160919.shp") %>%
  st_set_crs(4326) %>%
  st_crop(xmin=-4, xmax=-3.1, ymin=57, ymax=58.7) %>% # just keep meygen
  st_transform(32630)

flow_rate_shape <- st_read("data/mapping/Tidal model/Voronoi Polygons - Flow.shp") %>% 
  st_set_crs(4326) %>%
  st_transform(32630)

ggplot()+
  layer_spatial(data=flow_rate_shape, aes(fill=spring, colour=spring))+
  scale_fill_viridis_c(name="Mean Spring Flow Rate (m/s)")+
  scale_colour_viridis_c()+
  #guides(color=guide_legend("spring"), fill = FALSE)+
  theme(axis.text = element_text(size=30, family="serif"),
        legend.text = element_text(family="serif"),
        legend.title = element_text(family="serif"))+
  theme_bw()+
  annotation_spatial(data=pent_lr)+
    NULL

################


################
## calculate angles

seal <- seal %>%
  arrange(tag, datetime)

### Geospace
seal_move <- move(x=seal$lon, y=seal$lat, time=as.POSIXct(seal$datetime, tz="UTC"), 
                  proj=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84"), data=seal, animal=seal$tag)

angles <- angle(seal_move) %>%
  lapply(function(x){c(0,x)}) # add 0 to first location of each tag as no orientation with 1st location

seal$angle_geo <- as.numeric(unlist(angles)) 
  seal$angle_geo <- ifelse(seal$angle_geo<0, seal$angle_geo+360, seal$angle_geo)

### Hydrospace
  seal_move_hydro <- move(x=seal$lon_hydr, y=seal$lat_hydr, time=as.POSIXct(seal$datetime, tz="UTC"), 
                          proj=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84"), data=seal, animal=seal$tag)
  
  angles_hydro <- angle(seal_move_hydro) %>%
    lapply(function(x){c(0,x)}) # add 0 to first location of each tag as no orientation with 1st location
  
  seal$angle_hydro <- as.numeric(unlist(angles_hydro)) 
  seal$angle_hydro <- ifelse(seal$angle_hydro<0, seal$angle_hydro+360, seal$angle_hydro)

  
###################
seal_sf <- st_as_sf(seal, coords=c("lon", "lat"), crs=4326) %>%
    st_transform(32630)

  seals_in_meygen_sf <- st_join(seal_sf, lease_site, join=st_within)
  seals_in_meygen_df <- seal[which(!is.na(seals_in_meygen_sf$Lease_Star), arr.ind=TRUE),] # have done all this fannying around as st_geomtery seems to round off the locations when retreiving from sf so subset the data frame before creating the snipped sf object
  seals_in_meygen_sf <- filter(seals_in_meygen_sf, !is.na(Lease_Star))
  
lease_grid <- seals_in_meygen_sf %>%
    st_make_grid(cellsize = 500) %>%
    st_sf(grid_id = 1:length(.)) 

seals_in_meygen_df$grid_id <- st_join(seals_in_meygen_sf, lease_grid, join = st_intersects)$grid_id

ggplot()+
  annotation_spatial(data=lease_grid) +
  layer_spatial(seals_in_meygen_sf) +
  NULL

##### next to do - extract the grid_points so we have something to map the annotations to 

grid_points <-  lease_grid %>%
  st_centroid() %>%
  st_coordinates() %>%
  as_tibble() %>%
  mutate(grid_id=1:nrow(.)) %>%
  st_as_sf(coords=c("X", "Y"), crs=32630) 
  # add locations of centroids so can map plots to specific locations

ggplot()+
  geom_sf(data=lease_grid)+
  geom_sf(data=grid_points, aes(colour=grid_id))+
   NULL
  
###################
# Plot geo-angles on grid  
###################

geo_angle_plots <- seals_in_meygen_df %>%
  nest(-grid_id) %>%
  mutate(plot = map2(data, grid_id, 
                     ~ ggplot(.x) +  
                       theme_bw() +
                       geom_histogram(aes(x=angle_geo))+
                       coord_polar(start=0, direction = 1)+
                       scale_x_continuous(breaks=seq(0, 360, by=90), expand=c(0,0), lim=c(0, 360))+
                       xlab("")+
                       ylab("")+
                       theme(axis.text.y = element_blank(),
                         axis.ticks = element_blank(),
                         plot.background = element_rect(fill = "transparent")) 
                     ))

angle_plot_annotations <- grid_points %>%
 bind_cols(as_tibble(st_coordinates(.))) %>%
  st_drop_geometry() %>%
  #select(grid_id, X, Y) %>%
  left_join(geo_angle_plots, by = "grid_id") 
  
  angle_plot_annotations <- angle_plot_annotations[!angle_plot_annotations$plot=="NULL",]

  angle_plot_annotations <- angle_plot_annotations %>%
    mutate(annotation = pmap(list(X, Y, plot),
                           ~ annotation_custom(ggplotGrob(..3),
                                               xmin = ..1 - 250, xmax = ..1 + 250,
                                               ymin = ..2 - 250, ymax = ..2 + 250))) %>%
  pull(annotation)

coords_for_plot <- st_coordinates(seals_in_meygen_sf) # need a buffer around the edges rather than just using ggspatial so this helps us call it easily inside the plot function

ggplot(grid_points)+
  annotation_spatial(data=flow_rate_shape, aes(fill=spring, colour=spring))+
  scale_fill_viridis_c()+
  scale_colour_viridis_c()+
  labs(fill = expression(paste("Peak Spring Flow Speed ", (m.s^-1))))+
  guides(colour=FALSE) +
  theme(legend.position = "bottom") +
  xlim(c(min(coords_for_plot[,1])-500, max(coords_for_plot[,1]+500)))+
  ylim(c(min(coords_for_plot[,2])-500, max(coords_for_plot[,2]+500)))+
  geom_sf(alpha=0)+
  angle_plot_annotations +
  geom_sf(data=pent_hr)+
  labs(title = "Movement angle", subtitle = "Geospace")+
  theme_classic()+
  theme(text = element_text(family="serif", size=24), legend.position = "bottom")+
  NULL
  
###################
# Plot hydro-angles on grid  
###################

hydro_angle_plots <- seals_in_meygen_df %>%
  nest(-grid_id) %>%
  mutate(plot = map2(data, grid_id, 
                     ~ ggplot(.x) +  
                       theme_bw() +
                       geom_histogram(aes(x=angle_hydro))+
                       coord_polar(start=0, direction = 1)+
                       scale_x_continuous(breaks=seq(0, 360, by=90), expand=c(0,0), lim=c(0, 360))+
                       xlab("")+
                       ylab("")+
                       theme(axis.text.y = element_blank(),
                             axis.ticks = element_blank(),
                             plot.background = element_rect(fill = "transparent")) 
  ))

angle_plot_annotations <- grid_points %>%
  bind_cols(as_tibble(st_coordinates(.))) %>%
  st_drop_geometry() %>%
  #select(grid_id, X, Y) %>%
  left_join(hydro_angle_plots, by = "grid_id") 

angle_plot_annotations <- angle_plot_annotations[!angle_plot_annotations$plot=="NULL",]

angle_plot_annotations <- angle_plot_annotations %>%
  mutate(annotation = pmap(list(X, Y, plot),
                           ~ annotation_custom(ggplotGrob(..3),
                                               xmin = ..1 - 250, xmax = ..1 + 250,
                                               ymin = ..2 - 250, ymax = ..2 + 250))) %>%
  pull(annotation)

coords_for_plot <- st_coordinates(seals_in_meygen_sf) # need a buffer around the edges rather than just using ggspatial so this helps us call it easily inside the plot function

ggplot(grid_points)+
  annotation_spatial(data=flow_rate_shape, aes(fill=spring, colour=spring))+
  scale_fill_viridis_c()+
  scale_colour_viridis_c()+
  labs(fill = expression(paste("Peak Spring Flow Speed ", (m.s^-1))))+
  guides(colour=FALSE) +
  xlim(c(min(coords_for_plot[,1])-500, max(coords_for_plot[,1]+500)))+
  ylim(c(min(coords_for_plot[,2])-500, max(coords_for_plot[,2]+500)))+
  geom_sf(alpha=0)+
  angle_plot_annotations +
  geom_sf(data=pent_hr)+
  labs(title = "Movement angle", subtitle = "Hydrospace")+
  theme_classic()+
  theme(text = element_text(family="serif", size=24), legend.position = "bottom")+
  NULL

###########################################################
# Extra Processing

seals_in_meygen_df <- seals_in_meygen_df %>%
  mutate(angle_180=ifelse(.$angle_hydro>180, .$angle_hydro-180, .$angle_hydro),
         current_180=ifelse(.$current_di>180, .$current_di-180, .$current_di)) %>%
    mutate(swimdiff = abs(current_180 - angle_180))


Tides11 <- read_csv("data/GillsBay2011_2012.csv") 
Tides16 <- read_csv("data/GillsBay2016_2018.csv") 
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
seals_in_meygen_df$NearestHW <- as.POSIXct(unlist(mapply(function(x) nearest.time(x, HW2016.2018), as.POSIXct(seals_in_meygen_df$datetime))), origin="1970-01-01", tz = "UTC")

# column of time around nearest high water
seals_in_meygen_df$TimeAroundHW <- difftime(as.POSIXct(seals_in_meygen_df$datetime), seals_in_meygen_df$NearestHW, units = "hours")
seals_in_meygen_df$TimeAroundHW <- as.numeric(seals_in_meygen_df$TimeAroundHW)


benthos <- raster::raster("data/Mapping/Benthic/D4_2018.asc/benthos.tif")

coordinates(seals_in_meygen_df) <- ~lon+lat

seals_in_meygen_df$bathymetry <- abs(raster::extract(benthos, seals_in_meygen_df))
seals_in_meygen_df$bathymetry <- ifelse(seals_in_meygen_df$bathymetry < 0, 0, seals_in_meygen_df$bathymetry)

seals_in_meygen_df <- as.data.frame(seals_in_meygen_df) %>%
  drop_na(bathymetry) %>% #some dives will be effectively on land given the reslution of the depth data
  filter(bathymetry>0)%>%
  filter(TimeAroundHW<250)

save(seals_in_meygen_df, file="data/mod-move-dat.Rd")
