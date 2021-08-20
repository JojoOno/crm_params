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
                         axis.ticks = element_blank())))

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
  xlim(c(min(coords_for_plot[,1])-500, max(coords_for_plot[,1]+500)))+
  ylim(c(min(coords_for_plot[,2])-500, max(coords_for_plot[,2]+500)))+
  geom_sf(alpha=0)+
  angle_plot_annotations +
  geom_sf(data=pent_hr)+
  theme_classic()+
  NULL
  
