require(tidyverse)
require(data.table)
require(janitor)
require(ggspatial)


load("Data/Depth Data/Code Output/UHFDive_cleaned_locs.Rd")
load("Data/Depth Data/Code Output/GSMDive_cleaned.Rd")

colnames(TDR.D.locs)
colnames(dive)

gsm <- dive %>%
  dplyr::select(ref, DIVE_DUR, MAX_DEP, lon, lat, D_DATE.x, depth1, depth2, depth3, depth4, depth5, depth6, depth7, 
         depth8, depth9) %>%
  setnames(old=c("D_DATE.x", "MAX_DEP", "DIVE_DUR"), 
           new=c("mid.dt", "max.depth", "dive.duration")) %>%
  dplyr::mutate(depth10 = NA, depth11 = NA)

uhf <- TDR.D.locs %>%
  dplyr::select(tag, x, y, max.depth, dive.duration, mid.dt, depth1, depth2, depth3, depth4, depth5, depth6, depth7,
         depth8, depth9, depth10, depth11) %>%
  setnames(old=c("tag", "x", "y"), new=c("ref", "lon", "lat"))

combined.dive <- rbind(gsm, uhf) %>%
  clean_names()

## check ##
pent_lr <- st_read("Data/Mapping/Scotland/MaximumResolution/GSSHS_British_Isles poly.shp")
pent_lr <- pent_lr %>% 
  st_set_crs(4326) 

ggplot()+
  geom_point(data=combined.dive, aes(x=lon, y=lat, col=max_depth)) + #UTM is in Km in the gps data set
  scale_colour_viridis_c()+
  annotation_spatial(pent_lr)+
  facet_wrap(~ref)

save(combined.dive, file="Data/Depth Data/Code Output/combinedDive.Rd")
