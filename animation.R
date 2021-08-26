require(moveVis)
require(gganimate)
require(lubridate)
require(ggspatial)

################################## 
#Trying with MoveVis
################################## 
seal_move <- seal %>%
  mutate(datetime=as.POSIXct(datetime)) %>%
  df2move(proj="+init=epsg:4326 +proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0", x="lon", y="lat", time="datetime", track_id="tag") %>%
  align_move(res=10, unit="mins")


frames <- frames_spatial(seal_move[1:1000], map_service = "carto", map_type = "voyager_no_labels", alpha = 0.5)%>%  # path
  add_labels(x = "Longitude", y = "Latitude") %>% # add some customizations
  add_northarrow(colour = "black", position = "bottomright") %>% 
  add_scalebar(colour = "black", position = "bottomleft") %>% 
  add_progress(size = 2)

animate_frames(frames[1:100], out_file = "plots/animation.gif", overwrite = TRUE)


##################################
#try with custom
################################## 

seal_move <- seal[1:100,] %>%
  dplyr::mutate(date_time = ymd_hms(datetime)) %>%
  dplyr::mutate(hour = lubridate::yday(date_time)) %>%
  sf::st_as_sf(crs = 4326, coords=c("lon", "lat")) %>% 
  sf::st_transform(crs=32630) %>%
  dplyr::arrange(tag,hour) %>% 
  dplyr::group_by(tag,hour) %>% 
  dplyr::summarise(do_union = FALSE)

anim_data <- seal_move

animated_plot <- ggplot() +
  annotation_spatial(pent_lr, fill = "grey", lwd = 0) +
  layer_spatial(seal_move, aes(color = tag)) +
  scale_color_viridis_d() + 
  theme(legend.position = "none") +
  gganimate::transition_time(hour) +
  shadow_wake(wake_length = 0.1)

gganimate::animate(animated_plot)
range(seal_move$hour)
