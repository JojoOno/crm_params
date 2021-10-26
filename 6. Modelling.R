require(ggplot2)
require(tidyverse)
require(ggspatial)
require(sf)
require(geepack)
require(mvtnorm)  
require(splines)
require(mgcv)
require(boot)
require(sp)
require(gratia)
require(grid)
require(gridExtra)

##load model data

load(file="data/depth-data/code-output/mod-dive-dat.Rd")
load(file="data/mod-move-dat.Rd")

################################################################
### Lease site specification 
################################################################

### Make meygen polygon and clip data to polygon
bbox_meygen <- c(xmin=-3.17, ymin=55.5, xmax=-3.08, ymax=58.9) # clip it so just using Meygen as the shapefile contains all tidal lease sites in Scotland
bbox_pentland <- c(xmin=-3.4, ymin=55.5, xmax=-2, ymax=58.9)

meygen_lease <- st_read("data/mapping/wave-and-tidal-proposed-sites/TCE_Lease_Tide_20160919.shp") %>%
  st_set_crs(4326) %>%
  st_crop(bbox_meygen)

pentland_lease <- st_read("data/mapping/wave-and-tidal-proposed-sites/TCE_Lease_Tide_20160919.shp") %>%
  st_set_crs(4326) %>%
  st_crop(bbox_pentland)

pent_lr <- st_read("data/mapping/scotland/maximum-resolution/GSSHS_British_Isles poly.shp") %>%
  st_set_crs(4326)
pent_hr <- st_read("data/mapping/scotland/maximum-resolution/ScotLatLon.shp") %>%
  st_set_crs(4326)

## check you've just got meygen
ggplot()+
  annotation_spatial(data=pent_lr) +
  geom_sf(data=meygen_lease, fill="red", alpha=0.5) +
  NULL

## check you've just got pentland
ggplot()+
  annotation_spatial(data=pent_lr) +
  geom_sf(data=pentland_lease, fill="red", alpha=0.5) +
  NULL

##m make sf object of points and select points inside meygen lease site

dive_sf <- combined.dive %>%
  st_as_sf(coords=c("lon", "lat"), crs=4326)

dives_in_meygen_sf <- st_join(dive_sf, meygen_lease, join=st_within)
dives_in_meygen_df <- combined.dive[which(!is.na(dives_in_meygen_sf$Lease_Star), arr.ind=TRUE),] # have done all this fannying around as st_geomtery seems to round off the locations when retreiving from sf so subset the data frame before creating the snipped sf object
dives_in_meygen_sf <- filter(dives_in_meygen_sf, !is.na(Lease_Star))

dives_in_pentland_sf <- st_join(dive_sf, pentland_lease, join=st_within)
dives_in_pentland_df <- combined.dive[which(!is.na(dives_in_pentland_sf$Lease_Star), arr.ind=TRUE),] # have done all this fannying around as st_geomtery seems to round off the locations when retreiving from sf so subset the data frame before creating the snipped sf object
dives_in_pentland_sf <- filter(dives_in_pentland_sf, !is.na(Lease_Star))


################################################################
### Simple plots to confirm all is well with data limits
################################################################

## all dives
ggplot()+
  annotation_spatial(data=pent_lr) +
  geom_sf(data=meygen,  alpha=0.5) +
  layer_spatial(data=dive_sf, aes(colour=max_depth))+
  scale_colour_viridis_c(name="Dive Depth (m)")+
  theme_bw()+
  theme(axis.text = element_text(size=16, family="serif"),
        legend.text = element_text(family="serif"),
        legend.title = element_text(family="serif"))+
  NULL

## meygen dives
ggplot()+
  annotation_spatial(data=pent_lr) +
  geom_sf(data=meygen_lease,  alpha=0.5, fill="red") +
  layer_spatial(data=dives_in_meygen_sf, aes(colour=max_depth))+
  scale_colour_viridis_c(name="Dive Depth (m)")+
  theme_bw()+
  theme(axis.text = element_text(size=16, family="serif"),
        legend.text = element_text(family="serif"),
        legend.title = element_text(family="serif"))+
  NULL

## pentland dives
ggplot()+
  annotation_spatial(data=pent_lr) +
  geom_sf(data=pentland_lease,  alpha=0.5, fill="red") +
  layer_spatial(data=dives_in_pentland_sf, aes(colour=max_depth))+
  scale_colour_viridis_c(name="Dive Depth (m)")+
  theme_bw()+
  theme(axis.text = element_text(size=16, family="serif"),
        legend.text = element_text(family="serif"),
        legend.title = element_text(family="serif"))+
  NULL

##################################################################################################################################
##################################################################################################################################
# Modeling  Modeling Modeling Modeling Modeling Modeling Modeling Modeling Modeling Modeling Modeling Modeling Modeling Modeling #
##################################################################################################################################
##################################################################################################################################


##################################################################################################################################
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Proportion of dives within the risk zone - Meygen
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##################################################################################################################################

fit.gam.1 <- mgcv::gam(prop_risk ~ s(TimeAroundHW, bs="cc") +
                         s(lon,lat, by=TimeAroundHW)+
                         te(lon,lat, by=TimeAroundHW)+
                         s(lon,lat),
          family = "binomial",
          data=dives_in_meygen_df,
          method="ML")

summary(fit.gam.1)
plot(fit.gam.1)

appraise(fit.gam.1)

## Cyclic spline for tidal smooth

#~~~~~~~~~~~~~~~~~~~~~~~~~
spl_bs <- cSplineDes(x = dives_in_meygen_df$TimeAroundHW, 
                     knots = seq(min(dives_in_meygen_df$TimeAroundHW),
                                 max(dives_in_meygen_df$TimeAroundHW),
                                 length.out = 8), 
                     ord = 4, derivs = 0)

colnames(spl_bs) <- paste0("cyclic.", 1:ncol(spl_bs))
df <- cbind(spl_bs, dives_in_meygen_df)

model <- reformulate(c(-1,
                       paste0("cyclic.", 1:(ncol(spl_bs))),
                       "splines::bs(bathymetry)", "splines::bs(lon+lat)"), 
                     response = "prop_risk")

fit1 <- geepack::geeglm(model,
                        family="binomial",
                        data=df,
                        id=as.factor(ref),
                        corstr = "independence")

acf(fit1$residuals)
summary(fit1)
aod::wald.test(coef(fit1))
## Partial effects
#~~~~~~~~~~~~~~~~~~~~~~~~~

quant.func<- function(x){quantile(x, probs=c(0.05,0.95))} # ci levels

coefs<-coefficients(fit1) # Grabs the coefficients from your model
  
varcov<- summary(fit1)$cov.unscaled # Grabs the variance-covariance matrix from your model
  varcov[lower.tri(varcov)] = t(varcov)[lower.tri(varcov)] # you will often get an error that your matrix is not symetric - it actually is but the the number of significant figures your covariance values go to is such that it sppears they are not. By definition the covariance matrix has to be symetric so this is just an R fault - this line ensures R knows the matrix is symetric.

BootstrapParameters1 <- rmvnorm(500, coefs, varcov)

par(mfrow=c(1,2))

### Benthos
start=8; finish=10; Variable=dives_in_meygen_df$bathymetry; xlabel="Benthos"; ylabel="Proportion of dive in risk zone"
  PlottingVar1<-seq(min(Variable), max(Variable), length=500)

CenterVar1<-model.matrix(fit1)[,c(1,start:finish)]*coef(fit1)[c(1,start:finish)]
  BootstrapCoefs1<-BootstrapParameters1[,c(1,start:finish)]

Basis1<-gam(rbinom(500,2,0.5)~s(PlottingVar1, k=4), fit=F, family=binomial)$X[,1:4]
  RealFit1<-Basis1%*%coef(fit1)[c(1,start:finish)]
  RealFitCenter1<-RealFit1-mean(CenterVar1)
  RealFitCenter1a<-exp(RealFitCenter1)/(1+exp(RealFitCenter1))

BootstrapFits1<-Basis1%*%t(BootstrapCoefs1)
  quant.func1<-function(x){quantile(x,probs=c(0.025, 0.975))}
  cis1<-apply(BootstrapFits1, 1, quant.func1)-mean(CenterVar1)
  cis1a<-inv.logit(cis1)

plot(PlottingVar1,(RealFitCenter1a), type="l", col=1,ylim=c(0, 1),xlab=xlabel, ylab=ylabel, xlim=c(min(PlottingVar1),max(PlottingVar1)), main="", cex.lab = 1.5, cex.axis=1.5)
  segments(PlottingVar1,(cis1a[1,]),PlottingVar1,(cis1a[2,]), col="grey", main = "Bathymetry relationship")
  lines(PlottingVar1,(RealFitCenter1a),lwd=2, col=1)
  rug(Variable)

benthos.pred <- ggplot()+
  geom_line(aes(x=PlottingVar1, y=RealFitCenter1a))+
  geom_ribbon(aes(x=PlottingVar1, ymin=cis1a[1,], ymax=cis1a[2,]), alpha=0.1)+
  theme_bw()+
  labs(x="Bathymetric Depth (m)", y="Propotion of Dive in Risk Zone")+
  theme(axis.text = element_text(size=24, family="serif"),
        axis.title = element_text(size=24, family="serif"),
        legend.text = element_text(family="serif"),
        legend.title = element_text(family="serif"))+
  geom_rug(aes(dives_in_meygen_df$bathymetry))+
  NULL

### Tide
start=1; finish=7; Variable=dives_in_meygen_df$TimeAroundHW; xlabel="Hours Around High Water"; ylabel="Proportion of dive in risk zone"
  PlottingVar2<-seq(min(Variable), max(Variable), length=500)

CenterVar2<-model.matrix(fit1)[,c(1,start:finish)]*coef(fit1)[c(1,start:finish)]
  BootstrapCoefs2<-BootstrapParameters1[,c(1,start:finish)]

Basis2<-gam(rbinom(500,2,0.5)~s(PlottingVar2, bs="cc", k=9), fit=F, family=binomial)$X
  RealFit2<-Basis2%*%coef(fit1)[c(1,start:finish)]
  RealFitCenter2<-RealFit2-mean(CenterVar2)
  RealFitCenter2a<-exp(RealFitCenter2)/(1+exp(RealFitCenter2))

BootstrapFits2<-Basis2%*%t(BootstrapCoefs2)
  quant.func1<-function(x){quantile(x,probs=c(0.025, 0.975))}
  cis2<-apply(BootstrapFits2, 1, quant.func1)-mean(CenterVar2)
  cis2a<-inv.logit(cis2)

plot(PlottingVar2,(RealFitCenter2a), type="l", col=1,ylim=c(0, 1),xlab=xlabel, ylab=ylabel, xlim=c(min(PlottingVar2),max(PlottingVar2)), main="", cex.lab = 1.5, cex.axis=1.5)
  segments(PlottingVar2,(cis2a[1,]),PlottingVar2,(cis2a[2,]), col="grey", main = "Bathymetry relationship")
  lines(PlottingVar1,(RealFitCenter2a),lwd=2, col=1)
  rug(Variable)


tide.pred <- ggplot()+
    geom_line(aes(x=PlottingVar2, y=RealFitCenter2a))+
    geom_ribbon(aes(x=PlottingVar2, ymin=cis2a[1,], ymax=cis2a[2,]), alpha=0.1)+
    theme_bw()+
    labs(x="Time Around High Water (hrs)", y="Propotion of Dive in Risk Zone")+
    theme(axis.text = element_text(size=24, family="serif"),
          axis.title.x = element_text(size=24, family="serif"),
          axis.title.y = element_blank(),
          legend.text = element_text(family="serif"),
          legend.title = element_text(family="serif"))+
  geom_rug(aes(dives_in_meygen_df$TimeAroundHW))+
    NULL

grid.arrange(benthos.pred, tide.pred, nrow=1)

### Location predictions

  predgrid <- data.frame(expand.grid(lat=seq(min(dives_in_meygen_df$lat-0.01), max(dives_in_meygen_df$lat+0.01), length.out = 50),
                                     lon=seq(min(dives_in_meygen_df$lon-0.01), max(dives_in_meygen_df$lon+0.01), length.out = 50),
                                     TimeAroundHW=seq(-6, 6, by=1)))
  
 predgrid_sf <- st_join(st_as_sf(predgrid, coords=c("lon", "lat"), crs=4326), meygen, join=st_within)
  predgrid_df <- predgrid[which(!is.na(predgrid_sf$Lease_Star), arr.ind=TRUE),] # have done all this fannying around as st_geomtery seems to round off the locations when retreiving from sf so subset the data frame before creating the snipped sf object

  benthos <- raster::raster("data/mapping/Benthic/D4_2018.asc/benthos.tif")
  
  coordinates(predgrid_df) <- ~ lon+lat
  
  predgrid_df$bathymetry <- abs(raster::extract(benthos, predgrid_df))
  predgrid_df$bathymetry <- ifelse(predgrid_df$bathymetry < 0, 0, predgrid_df$bathymetry)
  predgrid_df <- na.omit(as.data.frame(predgrid_df))
  
  new_spline <- cSplineDes(predgrid_df$TimeAroundHW, 
                           knots = seq(min(df$TimeAroundHW),
                                       max(df$TimeAroundHW),
                                       length.out = 8),
                           ord = 4, derivs=0)
  
  colnames(new_spline) <- paste0("cyclic.", 1:ncol(new_spline))
  
  newdf <- cbind(new_spline, predgrid_df)
  
  plot_tide_time <- c(-6, -3, 0, 3)
  plot_tide_names <- c("LowWater", "Flood", "HighWater", "Ebb")
  
  for(i in 1:length(plot_tide_time)) {
 
  tide_time <- plot_tide_time[i]
  
  tide_predgrid <- newdf[newdf$TimeAroundHW==tide_time,]
  tide_predgrid$preds <- predict(fit1, tide_predgrid, type="response")

  
 p <- ggplot()+
    geom_tile(data=tide_predgrid, aes(x=lon, y=lat, fill=preds)) +
    scale_fill_viridis_c(limits=c(0,1)) +
   # annotation_spatial(data=pent_hr)+
    theme_bw() +
    ggtitle(paste(plot_tide_names[i])) +
   labs(fill = "")+
   theme(axis.text = element_text(size=16, family="serif"),
         axis.title.x = element_text(size=16, family="serif"),
         axis.title.y = element_blank(),
         legend.text = element_text(family="serif"),
         legend.title = element_text(family="serif"))+
   NULL
   
 
 assign(paste(plot_tide_names[i]), p)
  
  }


gridExtra::grid.arrange(LowWater, Flood, HighWater, Ebb, top=textGrob("Predicted proportion of time in risk zone", gp = gpar(fontsize = 20, fontface = 'bold')))  
  


##################################################################################################################################
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Proportion of dives within the risk zone - Pentland
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##################################################################################################################################
dives_in_pentland_df <- arrange(dives_in_pentland_df, ref)

fit.gam.1 <- mgcv::gam(prop_risk ~ s(TimeAroundHW, bs="cc") +
                         s(lon,lat, by=TimeAroundHW)+
                         te(lon,lat, by=TimeAroundHW)+
                         s(lon,lat),
                       family = "binomial",
                       data=dives_in_pentland_df,
                       method="ML")

summary(fit.gam.1)
plot(fit.gam.1)

appraise(fit.gam.1)
acf(resid(fit.gam.1))

## Cyclic spline for tidal smooth

#~~~~~~~~~~~~~~~~~~~~~~~~~
spl_bs <- cSplineDes(x = dives_in_pentland_df$TimeAroundHW, 
                     knots = seq(min(dives_in_pentland_df$TimeAroundHW),
                                 max(dives_in_pentland_df$TimeAroundHW),
                                 length.out = 8), 
                     ord = 4, derivs = 0)

colnames(spl_bs) <- paste0("cyclic.", 1:ncol(spl_bs))
df <- cbind(spl_bs, dives_in_pentland_df)

model <- reformulate(c(-1,
                       paste0("cyclic.", 1:(ncol(spl_bs))),
                       "splines::bs(bathymetry)", "splines::bs(lon+lat)"), 
                     response = "prop_risk")

fit1 <- geepack::geeglm(model,
                        family="binomial",
                        data=df,
                        id=as.factor(ref),
                        corstr = "independence")

acf(fit1$residuals)
summary(fit1)
aod::wald.test(coef(fit1))
## Partial effects
#~~~~~~~~~~~~~~~~~~~~~~~~~

quant.func<- function(x){quantile(x, probs=c(0.05,0.95))} # ci levels

coefs<-coefficients(fit1) # Grabs the coefficients from your model

varcov<- summary(fit1)$cov.unscaled # Grabs the variance-covariance matrix from your model
varcov[lower.tri(varcov)] = t(varcov)[lower.tri(varcov)] # you will often get an error that your matrix is not symetric - it actually is but the the number of significant figures your covariance values go to is such that it sppears they are not. By definition the covariance matrix has to be symetric so this is just an R fault - this line ensures R knows the matrix is symetric.

BootstrapParameters1 <- rmvnorm(500, coefs, varcov)

par(mfrow=c(1,2))

### Benthos
start=8; finish=10; Variable=dives_in_pentland_df$bathymetry; xlabel="Benthos"; ylabel="Proportion of dive in risk zone"
PlottingVar1<-seq(min(Variable), max(Variable), length=500)

CenterVar1<-model.matrix(fit1)[,c(1,start:finish)]*coef(fit1)[c(1,start:finish)]
BootstrapCoefs1<-BootstrapParameters1[,c(1,start:finish)]

Basis1<-gam(rbinom(500,2,0.5)~s(PlottingVar1, k=4), fit=F, family=binomial)$X[,1:4]
RealFit1<-Basis1%*%coef(fit1)[c(1,start:finish)]
RealFitCenter1<-RealFit1-mean(CenterVar1)
RealFitCenter1a<-exp(RealFitCenter1)/(1+exp(RealFitCenter1))

BootstrapFits1<-Basis1%*%t(BootstrapCoefs1)
quant.func1<-function(x){quantile(x,probs=c(0.025, 0.975))}
cis1<-apply(BootstrapFits1, 1, quant.func1)-mean(CenterVar1)
cis1a<-inv.logit(cis1)

plot(PlottingVar1,(RealFitCenter1a), type="l", col=1,ylim=c(0, 1),xlab=xlabel, ylab=ylabel, xlim=c(min(PlottingVar1),max(PlottingVar1)), main="", cex.lab = 1.5, cex.axis=1.5)
segments(PlottingVar1,(cis1a[1,]),PlottingVar1,(cis1a[2,]), col="grey", main = "Bathymetry relationship")
lines(PlottingVar1,(RealFitCenter1a),lwd=2, col=1)
rug(Variable)

benthos.pred <- ggplot()+
  geom_line(aes(x=PlottingVar1, y=RealFitCenter1a))+
  geom_ribbon(aes(x=PlottingVar1, ymin=cis1a[1,], ymax=cis1a[2,]), alpha=0.1)+
  theme_bw()+
  labs(x="Bathymetric Depth (m)", y="Propotion of Dive in Risk Zone")+
  theme(axis.text = element_text(size=24, family="serif"),
        axis.title = element_text(size=24, family="serif"),
        legend.text = element_text(family="serif"),
        legend.title = element_text(family="serif"))+
  geom_rug(aes(dives_in_pentland_df$bathymetry))+
  NULL

### Tide
start=1; finish=7; Variable=dives_in_pentland_df$TimeAroundHW; xlabel="Hours Around High Water"; ylabel="Proportion of dive in risk zone"
PlottingVar2<-seq(min(Variable), max(Variable), length=500)

CenterVar2<-model.matrix(fit1)[,c(1,start:finish)]*coef(fit1)[c(1,start:finish)]
BootstrapCoefs2<-BootstrapParameters1[,c(1,start:finish)]

Basis2<-gam(rbinom(500,2,0.5)~s(PlottingVar2, bs="cc", k=9), fit=F, family=binomial)$X
RealFit2<-Basis2%*%coef(fit1)[c(1,start:finish)]
RealFitCenter2<-RealFit2-mean(CenterVar2)
RealFitCenter2a<-exp(RealFitCenter2)/(1+exp(RealFitCenter2))

BootstrapFits2<-Basis2%*%t(BootstrapCoefs2)
quant.func1<-function(x){quantile(x,probs=c(0.025, 0.975))}
cis2<-apply(BootstrapFits2, 1, quant.func1)-mean(CenterVar2)
cis2a<-inv.logit(cis2)

plot(PlottingVar2,(RealFitCenter2a), type="l", col=1,ylim=c(0, 1),xlab=xlabel, ylab=ylabel, xlim=c(min(PlottingVar2),max(PlottingVar2)), main="", cex.lab = 1.5, cex.axis=1.5)
segments(PlottingVar2,(cis2a[1,]),PlottingVar2,(cis2a[2,]), col="grey", main = "Bathymetry relationship")
lines(PlottingVar1,(RealFitCenter2a),lwd=2, col=1)
rug(Variable)


tide.pred <- ggplot()+
  geom_line(aes(x=PlottingVar2, y=RealFitCenter2a))+
  geom_ribbon(aes(x=PlottingVar2, ymin=cis2a[1,], ymax=cis2a[2,]), alpha=0.1)+
  theme_bw()+
  labs(x="Time Around High Water (hrs)", y="Propotion of Dive in Risk Zone")+
  theme(axis.text = element_text(size=24, family="serif"),
        axis.title.x = element_text(size=24, family="serif"),
        axis.title.y = element_blank(),
        legend.text = element_text(family="serif"),
        legend.title = element_text(family="serif"))+
  geom_rug(aes(dives_in_pentland_df$TimeAroundHW))+
  NULL

grid.arrange(benthos.pred, tide.pred, nrow=1)

### Location predictions

predgrid <- data.frame(expand.grid(lat=seq(min(dives_in_pentland_df$lat-0.01), max(dives_in_pentland_df$lat+0.01), length.out = 150),
                                   lon=seq(min(dives_in_pentland_df$lon-0.01), max(dives_in_pentland_df$lon+0.01), length.out = 150),
                                   TimeAroundHW=seq(-6, 6, by=1)))

predgrid_sf <- st_join(st_as_sf(predgrid, coords=c("lon", "lat"), crs=4326), meygen, join=st_within)
predgrid_df <- predgrid[which(!is.na(predgrid_sf$Lease_Star), arr.ind=TRUE),]
predgrid_df <-mutate(predgrid_df, name=predgrid_sf$Name_Ten[which(!is.na(predgrid_sf$Lease_Star), arr.ind=TRUE)])# have done all this fannying around as st_geomtery seems to round off the locations when retreiving from sf so subset the data frame before creating the snipped sf object

benthos <- raster::raster("data/mapping/Benthic/D4_2018.asc/benthos.tif")

coordinates(predgrid_df) <- ~ lon+lat

predgrid_df$bathymetry <- abs(raster::extract(benthos, predgrid_df))
predgrid_df$bathymetry <- ifelse(predgrid_df$bathymetry < 0, 0, predgrid_df$bathymetry)
predgrid_df <- na.omit(as.data.frame(predgrid_df))

new_spline <- cSplineDes(predgrid_df$TimeAroundHW, 
                         knots = seq(min(df$TimeAroundHW),
                                     max(df$TimeAroundHW),
                                     length.out = 8),
                         ord = 4, derivs=0)

colnames(new_spline) <- paste0("cyclic.", 1:ncol(new_spline))

newdf <- cbind(new_spline, predgrid_df)

plot_tide_time <- c(-6, -3, 0, 3)
plot_tide_names <- c("LowWater", "Flood", "HighWater", "Ebb")

newdf <- newdf %>%
  mutate(name_abr = case_when(
  name == unique(name)[1] ~ "Brims",
  name == unique(name)[2] ~ "Meygen",
  name == unique(name)[3] ~ "SPR",
  TRUE ~ as.character(name)
))

sites <- unique(newdf$name_abr)


for(i in 1:length(sites)) {
  
  site_spec <- sites[i]
  
  for(j in 1:length(plot_tide_time)){
  tide_time <- plot_tide_time[j]
  
  tide_predgrid <- newdf[newdf$TimeAroundHW==tide_time & newdf$name_abr==site_spec,]
  tide_predgrid$preds <- predict(fit1, tide_predgrid, type="response")
  
  
  p <- ggplot()+
    geom_tile(data=tide_predgrid, aes(x=lon, y=lat, fill=preds)) +
    scale_fill_viridis_c(limits=c(0,1), option="cividis") +
    # annotation_spatial(data=pent_hr)+
    theme_bw() +
    ggtitle(paste(plot_tide_names[j])) +
    labs(fill = "")+
    theme(axis.text = element_text(size=16, family="serif"),
          axis.title.x = element_text(size=16, family="serif"),
          axis.title.y = element_blank(),
          legend.text = element_text(family="serif"),
          legend.title = element_text(family="serif"))+
    NULL
  
  
  assign(paste(plot_tide_names[j],site_spec, sep=""), p)
  }
}


gridExtra::grid.arrange(LowWaterBrims, FloodBrims, HighWaterBrims, EbbBrims, top=textGrob("Predicted proportion of time in risk zone", gp = gpar(fontsize = 20, fontface = 'bold')))  
gridExtra::grid.arrange(LowWaterMeygen, FloodMeygen, HighWaterMeygen, EbbMeygen, top=textGrob("Predicted proportion of time in risk zone", gp = gpar(fontsize = 20, fontface = 'bold')))  
gridExtra::grid.arrange(LowWaterSPR, FloodSPR, HighWaterSPR, EbbSPR, top=textGrob("Predicted proportion of time in risk zone", gp = gpar(fontsize = 20, fontface = 'bold')))  



##################################################################################################################################
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Orientation
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##################################################################################################################################

fit.gam.1 <- mgcv::gam(swimdiff ~ s(TimeAroundHW, bs="cc") +
                         s(lon,lat, by=TimeAroundHW)+
                         te(lon,lat, by=TimeAroundHW)+
                         s(lon,lat),
                       family = "Gamma",
                       data=seals_in_meygen_df,
                       method="ML")

summary(fit.gam.1)
plot(fit.gam.1)

appraise(fit.gam.1)

## Cyclic spline for tidal smooth

#~~~~~~~~~~~~~~~~~~~~~~~~~
spl_bs <- cSplineDes(x = seals_in_meygen_df$TimeAroundHW, 
                     knots = seq(min(seals_in_meygen_df$TimeAroundHW),
                                 max(seals_in_meygen_df$TimeAroundHW),
                                 length.out = 8), 
                     ord = 4, derivs = 0)

colnames(spl_bs) <- paste0("cyclic.", 1:ncol(spl_bs))
df <- cbind(spl_bs, seals_in_meygen_df)

model <- reformulate(c(-1,
                       paste0("cyclic.", 1:(ncol(spl_bs))),
                       "splines::bs(bathymetry)", "splines::bs(lon+lat)"), 
                     response = "swimdiff")

fit1 <- geepack::geeglm(model,
                        family="Gamma",
                        data=df,
                        id=as.factor(tag),
                        corstr = "independence")


acf(fit1$residuals)
summary(fit1)
## Partial effects
#~~~~~~~~~~~~~~~~~~~~~~~~~

quant.func<- function(x){quantile(x, probs=c(0.05,0.95))} # ci levels

coefs<-coefficients(fit1) # Grabs the coefficients from your model

varcov<- summary(fit1)$cov.unscaled # Grabs the variance-covariance matrix from your model
varcov[lower.tri(varcov)] = t(varcov)[lower.tri(varcov)] # you will often get an error that your matrix is not symetric - it actually is but the the number of significant figures your covariance values go to is such that it sppears they are not. By definition the covariance matrix has to be symetric so this is just an R fault - this line ensures R knows the matrix is symetric.

BootstrapParameters1 <- rmvnorm(500, coefs, varcov)

par(mfrow=c(1,2))

### Tide (bathymetry was, unsurprisingly, no effect on movement direction)
start=1; finish=7; Variable=seals_in_meygen_df$TimeAroundHW; xlabel="Hours Around High Water"; ylabel="Difference in Swim Angle"
PlottingVar2<-seq(min(Variable), max(Variable), length=500)

CenterVar2<-model.matrix(fit1)[,c(1,start:finish)]*coef(fit1)[c(1,start:finish)]
BootstrapCoefs2<-BootstrapParameters1[,c(1,start:finish)]

Basis2<-gam(abs(rnorm(500,mean(seals_in_meygen_df$swimdiff),sd(seals_in_meygen_df$swimdiff)))~s(PlottingVar2, bs="cc", k=9), fit=F, family=Gamma)$X
RealFit2<-Basis2%*%coef(fit1)[c(1,start:finish)]
RealFitCenter2<-RealFit2-mean(CenterVar2)
RealFitCenter2a<-exp(RealFitCenter2)/(1+exp(RealFitCenter2))

BootstrapFits2<-Basis2%*%t(BootstrapCoefs2)
quant.func1<-function(x){quantile(x,probs=c(0.025, 0.975))}
cis2<-apply(BootstrapFits2, 1, quant.func1)-mean(CenterVar2)
cis2a<-inv.logit(cis2)

plot(PlottingVar2,(RealFitCenter2), type="l", col=1,xlab=xlabel, ylab=ylabel, xlim=c(min(PlottingVar2),max(PlottingVar2)), main="", cex.lab = 1.5, cex.axis=1.5)
segments(PlottingVar2,(cis2[1,]),PlottingVar2,(cis2[2,]), col="grey", main = "Bathymetry relationship")
lines(PlottingVar1,(RealFitCenter2),lwd=2, col=1)
rug(Variable)


tide.pred <- ggplot()+
  geom_line(aes(x=PlottingVar2, y=RealFitCenter2*1000+50))+
  geom_ribbon(aes(x=PlottingVar2, ymin=cis2[1,]*1000+50, ymax=cis2[2,]*1000+50), alpha=0.1)+
  theme_bw()+
  labs(x="Time Around High Water (hrs)", y="Difference in Swim Angle")+
  theme(axis.text = element_text(size=24, family="serif"),
        axis.title.x = element_text(size=24, family="serif"),
        axis.title.y = element_text(size=24, family="serif"),
        legend.text = element_text(family="serif"),
        legend.title = element_text(family="serif"))+
  geom_rug(aes(dives_in_meygen_df$TimeAroundHW))+
  NULL

tide.pred

### Location predictions

predgrid <- data.frame(expand.grid(lat=seq(min(dives_in_meygen_df$lat-0.01), max(dives_in_meygen_df$lat+0.01), length.out = 50),
                                   lon=seq(min(dives_in_meygen_df$lon-0.01), max(dives_in_meygen_df$lon+0.01), length.out = 50),
                                   TimeAroundHW=seq(-6, 6, by=1)))

predgrid_sf <- st_join(st_as_sf(predgrid, coords=c("lon", "lat"), crs=4326), meygen, join=st_within)
predgrid_df <- predgrid[which(!is.na(predgrid_sf$Lease_Star), arr.ind=TRUE),] # have done all this fannying around as st_geomtery seems to round off the locations when retreiving from sf so subset the data frame before creating the snipped sf object

benthos <- raster::raster("data/mapping/Benthic/D4_2018.asc/benthos.tif")

coordinates(predgrid_df) <- ~ lon+lat

predgrid_df$bathymetry <- abs(raster::extract(benthos, predgrid_df))
predgrid_df$bathymetry <- ifelse(predgrid_df$bathymetry < 0, 0, predgrid_df$bathymetry)
predgrid_df <- na.omit(as.data.frame(predgrid_df))

new_spline <- cSplineDes(predgrid_df$TimeAroundHW, 
                         knots = seq(min(df$TimeAroundHW),
                                     max(df$TimeAroundHW),
                                     length.out = 8),
                         ord = 4, derivs=0)

colnames(new_spline) <- paste0("cyclic.", 1:ncol(new_spline))

newdf <- cbind(new_spline, predgrid_df)

plot_tide_time <- c(-6, -3, 0, 3)
plot_tide_names <- c("LowWater", "Flood", "HighWater", "Ebb")

for(i in 1:length(plot_tide_time)) {
  
  tide_time <- plot_tide_time[i]
  
  tide_predgrid <- newdf[newdf$TimeAroundHW==tide_time,]
  tide_predgrid$preds <- predict(fit1, tide_predgrid, type="response")
  tide_predgrid$preds <- ifelse(tide_predgrid$preds>180, 180, tide_predgrid$preds)
  
  
  p <- ggplot()+
    geom_tile(data=tide_predgrid, aes(x=lon, y=lat, fill=abs(preds))) +
    scale_fill_viridis_c(limits=c(0,180)) +
    # annotation_spatial(data=pent_hr)+
    theme_bw() +
    ggtitle(paste(plot_tide_names[i])) +
    labs(fill = "")+
    theme(axis.text = element_text(size=16, family="serif"),
          axis.title.x = element_text(size=16, family="serif"),
          axis.title.y = element_blank(),
          legend.text = element_text(family="serif"),
          legend.title = element_text(family="serif"))+
    NULL
  
  
  assign(paste(plot_tide_names[i]), p)
  
}


gridExtra::grid.arrange(LowWater, Flood, HighWater, Ebb, top=textGrob("Difference in Swim Angle", gp = gpar(fontsize = 20, fontface = 'bold')))  
