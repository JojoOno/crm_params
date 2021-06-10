require(tidyverse)
require(data.table)
require(sf)
require(ggspatial)


load(file="Data/Depth Data/TDR2016.Rd")
load(file="Data/Depth Data/TDR2017.Rd")
load(file="Data/Depth Data/TDR2018.Rd")

################################################################################
# Make one data frame from individual dive records
# Identify discrete dives based on depth thresholds ( a la GSM tags )
################################################################################
seals <- unique(atsea.flow$Tag)    
for (i in 1:length(seals)) {
  d <- subset(atsea.flow, atsea.flow$Tag==seals[i])
  dsn <- paste("Tag", seals[i], sep="")
  assign(dsn, d)
}

tdrs16 <- unique(TDR2016$tag)
for (i in 1:length(tdrs16)) {
  d <- subset(TDR2016, TDR2016$tag==tdrs16[i])
  dsn <- paste("tdr", tdrs16[i], sep="")
  assign(dsn, d)
}

tdrs17 <- unique(TDR2017$tag)
for (i in 1:length(tdrs17)) {
  d <- subset(TDR2017, TDR2017$tag==tdrs17[i])
  dsn <- paste("tdr", tdrs17[i], sep="")
  assign(dsn, d)
}

tdrs18 <- unique(TDR2018$tag)
for (i in 1:length(tdrs18)) {
  d <- subset(TDR2018, TDR2018$tag==tdrs18[i])
  dsn <- paste("tdr", tdrs18[i], sep="")
  assign(dsn, d)
}


TDR2016$`Dive Number`<- NULL
TDR2017$Voltage <- NULL
TDR2018$Voltage <- NULL
TDR2016$datetime2 <- NULL

TDR <- rbind(TDR2016, TDR2017, TDR2018)

dive <- vector(mode="numeric", length=nrow(TDR))
dive <- ifelse(TDR$Depth > 8, 1, 0) # 8m threshold for defining dives
head(dive, n=600)
tail(dive, n=600)

dive.surf <- vector(mode="character", length=nrow(TDR))
dive.start <- vector(mode="character", length=nrow(TDR)-1)
dive.end <- vector(mode="character", length=nrow(TDR)-1)

dive1 <- dive[1:(length(dive)-1)]; length(dive1)
dive2 <- dive[2:length(dive)]; length(dive2)
dive.lag <- cbind(dive1,dive2)

dive.start <- ifelse(dive.lag[,1]==0 & dive.lag[,2]==1, 1, 0) # pick out the beginning of dives
dive.end <- ifelse(dive.lag[,1]==1 & dive.lag[,2]==0, 1, 0) # pick out the end of dives
dive.start.index <- which(dive.start==1); dive.start.index <- dive.start.index+1; length(dive.start.index)
dive.end.index <- which(dive.end==1); length(dive.end.index)

if(dive[1]==1){ tt <- c(1, dive.start.index) } else { tt <- dive.start.index } ; tt
if(dive[length(dive)]==1){ ss <- c(dive.end.index, length(dive)) } else { ss <- dive.end.index } ; ss
dive.index <- cbind(tt, ss); dive.index

surf.start.index <- which(dive.end==1)+1
surf.end.index <- which(dive.start==1) 
if(dive[1]==0){ kk <- c(1, surf.start.index) } else { kk <- surf.start.index } ; kk
if(dive[length(dive)]==0){ jj <- c(surf.end.index, length(dive)) } else { jj <- surf.end.index } ; jj
surf.index <- cbind(kk, jj); surf.index

for (i in 1:nrow(dive.index))
{
  dive.surf[c(dive.index[i,1]:dive.index[i,2])] <- paste("dive",i, sep="")
}

for (i in 1:nrow(surf.index)){
  
  dive.surf[c(surf.index[i,1]:surf.index[i,2])] <- paste("surf",i+1, sep="")
}

TDR$DIVE_SURF <- dive.surf

save(TDR, file=paste("Data/Depth Data/TDRFULL.Rd"))

TDR <- na.omit(TDR)
rm(list=ls(pattern="Tag"))
rm(list=ls(pattern="tdr"))

################################################################################
# Broken stick model to allocate inflection points to continuous dive data 
# Photopolou et al. (2015) Methods in Ecology and Evolution 
################################################################################

load(file=paste("Data/Depth Data/TDRFULL.Rd"))

dives <- unique(TDR$DIVE_SURF)

detailed.dive.data <- TDR

detailed.dive.data$DIVE_SURF <- as.factor(detailed.dive.data$DIVE_SURF)
detailed.dive.data <- detailed.dive.data[grep("dive", detailed.dive.data$DIVE_SURF),]

#detailed.dive.data <- detailed.dive.data[1:1000,]
## Remove short dives

remove.shorts <- function(x){
  xdf <- data.frame(x)
  if(nrow(xdf)>9){return(x)}else(return(NA))
}

templist <- tapply(detailed.dive.data$Depth, detailed.dive.data$DIVE_SURF, remove.shorts)
templist <- which(is.na(templist))

null.dives <- names(templist)

detailed.dive.data <- detailed.dive.data %>%
  filter(!(DIVE_SURF%in%null.dives))


get.BSMdive2 <- purrr::possibly(get.BSMdive, otherwise=NA) ##alternative to tryCatch using purrr -- if an error is 

BSM.list <- tapply(detailed.dive.data$Depth, detailed.dive.data$DIVE_SURF, FUN=get.BSMdive2, 
                   BS.dive=FALSE, res=9, breakpoints=9, dev.new=FALSE, plot.truth=FALSE, plot.bsm=FALSE, 
                   draw.numbers=FALSE, draw.lines=FALSE)


BSM.list <- BSM.list[lapply(BSM.list,length)>0]; BSM.list <- BSM.list[!is.na(BSM.list)]

d.nu <- "dive1000"

plot(BSM.list[[d.nu]]$time.depth[,2]*-1~BSM.list[[d.nu]]$time.depth[,1], type='l')
points(BSM.list[[d.nu]]$BS.time.depth[,2]*-1~BSM.list[[d.nu]]$BS.time.depth[,1], col="red")



d.nest.list <- vector(mode = "list")
for(i in 1:length(BSM.list)) {
  
  d.list <- find.Xzone.Joe(output=BSM.list[[i]], BS.dive=F, plot.truth=FALSE, plot = FALSE)
  d.list$diveRecord <- names(BSM.list)[i]
  d.nest.list[[i]] <- d.list  
  print(i)
}

d.nu <- runif(1,1,length(d.nest.list))
plot(-d.nest.list[[d.nu]]$bs.points[,2]~d.nest.list[[d.nu]]$bs.points[,1], type="l")

save(d.nest.list, file="Data/Depth Data/Code Output/UHFDive_nestedlist.Rd")

################################################################################
# Create data frame where one observation = one dive so directly comparable with GSM format
################################################################################

load(file="Data/Depth Data/Code Output/GSMDive_cleaned.Rd")
dive.gsm <- dive; rm(dive)
colnames(dive.gsm)
  
TDR.D <- setNames(data.frame(matrix(ncol = 6, nrow=length(d.nest.list))), 
                  c("ref","dive.id", "start.dt", "end.dt", "dive.duration", 
                    "max.depth")) 
tds <- paste0("depth", 1:11, sep="")
TDR.D[,tds] <- NA

for (j in 1:length(d.nest.list)){
  
  dive <- d.nest.list[[j]]
  
    TDR.D$ref[j] <- detailed.dive.data$tag[which(detailed.dive.data$DIVE_SURF==dive$diveRecord)[1]]
    TDR.D$dive.id[j] <- dive$diveRecord
    
    
      
      TDR.D$start.dt[j] <- detailed.dive.data$datetime[which(detailed.dive.data$DIVE_SURF==dive$diveRecord)] %>% min() %>% as.POSIXct(as.character(), origin = "1970-01-01", tz="UTC")
      TDR.D$end.dt[j] <- detailed.dive.data$datetime[which(detailed.dive.data$DIVE_SURF==dive$diveRecord)] %>% max() %>% as.POSIXct(as.character(), origin = "1970-01-01", tz="UTC")
      TDR.D$dive.duration[j] <- max(dive$bs.points[,1])*60
      TDR.D$max.depth[j] <- max(dive$bs.points[,2])
      
      
      TDR.D[j,7:17] <- d.nest.list[[j]]$bs.points$depth
     
      print(j)
   
 } 


TDR.D$start.dt <- as.POSIXct(TDR.D$start.dt, format="%Y-%m-%d %H:%M:%S", origin = "1970-01-01 00:00:00") %>%
  format(tz="UTC") %>%
  as.POSIXct()

TDR.D$end.dt <- as.POSIXct(TDR.D$end.dt, format="%Y-%m-%d %H:%M:%S", origin = "1970-01-01 00:00:00") %>%
  format(tz="UTC") %>%
  as.POSIXct()

TDR.D$mid.dt <- TDR.D$end.dt - (TDR.D$dive.duration/2) 

TDR.D <- TDR.D[!is.na(TDR.D$mid.dt), ]


save(TDR.D, file="Data/Depth Data/Code Output/UHFDive_cleaned.Rd")
save(detailed.dive.data, file="Data/Depth Data/Code Output/UHFDive_detailed.Rd")

################################################################################
# Allocate locations to dives
################################################################################

load(file="Data/GPS-COMBINED-FLOW.Rd")
load(file="Data/Depth Data/Code Output/UHFDive_cleaned.Rd")
captures <- read_csv(file="Data/Basic Capture Info.csv") # need the capture data base so we can lookup the names of the gps tags that go nwith the TDRs and match the correct tags and times to each other to get locations

uhf_gps <- filter(gps_flow_comb, datetime>as.POSIXct("2015-01-01 00:00:00")) #only want UHF data for this processing

colnames(captures)[5:6] <- c("tag", "ref") 
TDR.D$tag <- NA
TDR.D <- merge(TDR.D, captures[, c("ref", "tag")], by="ref", all.x = TRUE) %>%
  select(-tag.x) %>%
  rename(tag = tag.y)

dive.table <- data.table(TDR.D) %>%
  mutate(time = mid.dt) %>%
  setkey(tag, time) %>%
  setDT


gps.table <- data.table(uhf_gps) %>%
  mutate(time=datetime) %>%
  setkey(tag, time) %>%
  setDT()

TDR.D.locs <- gps.table[dive.table, roll="nearest"] %>%
  drop_na(lat) %>%
  mutate(difference = as.numeric(abs((datetime-mid.dt)/60))) %>%
  filter(difference<31)

save(TDR.D.locs, file="Data/Depth Data/Code Output/UHFDive_cleaned_locs.Rd")

pent_lr <- st_read("Data/Mapping/Scotland/MaximumResolution/GSSHS_British_Isles poly.shp")
pent_lr <- pent_lr %>% 
  st_set_crs(4326) %>%
  st_transform(32630)

ggplot()+
  geom_point(data=TDR.D.locs, aes(x=lonUTM*1000, y=latUTM*1000, col=max.depth)) + #UTM is in Km in the gps data set
  scale_colour_viridis_c()+
  annotation_spatial(pent_lr)+
  facet_wrap(~tag)


