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
  
pent_grid <- seal_sf %>%
    st_make_grid(cellsize = 10000) %>%
    st_sf(grid_id = 1:length(.))

#####next to do - extract the grid_points so we have something to map the annotations to 

grid_points <-  pent_grid %>%
  st_centroid() %>%
  st_coordinates() %>%
  st_sfc() %>%
  st_sf()# add locations of centroids so can map plots to specific locations

ggplot()+
    layer_spatial(data=st_sf(points=st_sfc(st_multipoint(matrix(as.numeric(pent_grid$x), as.numeric(pent_grid$y), ,2)), crs=32630), aes(fill=grid_id)) +
    annotation_spatial(data=pent_lr)
  
grid_id <- seal_sf %>% st_join(pent_grid, join = st_intersects) %>% as.data.frame()
  seal$grid_id <- grid_id$grid_id
  seal_sf$grid_id <- grid_id$grid_id
  
###################
  
  seals_in_meygen_sf <- st_join(seal_sf, lease_site, join=st_within)
  seals_in_meygen_df <- seal[which(!is.na(seals_in_meygen_sf$Lease_Star), arr.ind=TRUE),] # have done all this fannying around as st_geomtery seems to round off the locations when retreiving from sf so subset the data frame before creating the snipped sf object
  seals_in_meygen_sf <- filter(seals_in_meygen_sf, !is.na(Lease_Star))
  

geo_angle_plots <- seals_in_meygen_df %>%
  nest(-grid_id) %>%
  mutate(plot = map2(data, grid_id, 
                     ~ ggplot(.x) + 
                       ggtitle(.y) + 
                       theme_bw() +
                       geom_histogram(aes(x=angle_geo))+
                       coord_polar()))


plot_locations <- seals_in_meygen_sf %>%
  st_make_grid(cellsize = 10000, what="centers") %>%
  st_sf(grid_id = 1:length(.))# get plot locations which is just the coordinates of the unique grid cells the plots lie in

ggplot()+geom_sf(data=plot_locations)

angle_plot_anotations <- plot_locations %>%
 bind_cols(as_tibble(st_coordinates(.))) %>%
  #st_drop_geometry() %>%
  select(grid_id, X, Y) %>%
  left_join(geo_angle_plots, by = "grid_id") %>%
  mutate(annotation = pmap(list(X, Y, plot),
                           ~ annotation_custom(ggplotGrob(..3),
                                               xmin = ..1 - 2000, xmax = ..1 + 2000,
                                               ymin = ..2 - 1000, ymax = ..2 + 1000))) %>%
  pull(annotation)

ggplot()+
    geom_histogram(data=seals_in_meygen_df, aes(x=angle_geo))+
    coord_polar()+
    facet_wrap(~grid_id)+
      theme_bw()



  