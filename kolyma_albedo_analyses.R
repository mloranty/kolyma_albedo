#################################
#
# analyses of post-fire 
# vegetation ~ albedo dynamics
# for Kolyma River basis
#
# MML 07/24/2020
#################################

### load required packages
library(raster)
library(sf)
library(dplyr)


### set working directory
setwd("C:/Users/mloranty/Documents/GitHub/kolyma_albedo/")

### define functions
# standard error of the mean
se <- function(x, ...) sqrt(var(x, ...)/length(x))

### calculate stand level biomass
# read raw density gradient data
dg <- read.csv("data/DensityGradientTrees_raw_2010_2017.csv", header = T)

# read raw Polaris data; archived at Arctic Data Center https://arcticdata.io/catalog/view/doi:10.5065/D6NG4NP0
pol <- read.csv("data/Polaris_TS_2012_2013_Master_trees.csv", header = T)
  
# calculate larch aboveground biomass using equations from Alexander et al 2012
# divide by plot area and multiply by 0.47 to calculate g C/m2
dg$biomassC <- ifelse(is.na(dg$BD..cm.),
                      0.47*((179.2*dg$DBH..cm.^2.01)/dg$Area.sampled..m2.),
                      0.47*((39.46*dg$BD..cm.^2.26)/dg$Area.sampled..m2.))

pol$biomassC <- ifelse(is.na(pol$BD),
                      0.47*((179.2*pol$DBH^2.01)/pol$Plot_area),
                      0.47*((39.46*pol$BD^2.26)/pol$Plot_area))

# sum to cumulative biomass C for each plot
dg.plot <- aggregate(dg$biomassC~ dg$Site+dg$Plot,FUN = "sum")
pol.plot <- aggregate(pol$biomassC~ pol$Site+pol$Trans,FUN = "sum")
#change header names
names(dg.plot) <- c("site", "plot", "biomassC")
names(pol.plot) <- c("site", "plot", "biomassC")

# calculate stand-level mean and SE 
dg.site <- aggregate(dg.plot$biomassC~dg.plot$site,FUN = "mean")
dg.site$se <- aggregate(dg.plot$biomassC~ dg.plot$site,FUN = "se")[,2]

pol.site <- aggregate(pol.plot$biomassC~pol.plot$site,FUN = "mean")
pol.site$se <- aggregate(pol.plot$biomassC~ pol.plot$site,FUN = "se")[,2]

# append data sets and clean up workspace

### calculate stand-level percent canopy cover from densiometer measurements


pol.dens <- read.csv("data/Polaris_TS_2012_2013_Master_densiometer.csv", header = T)

# calculate plot means
pol.dens.plot <- aggregate(pol.dens$PctCover~pol.dens$Site+pol.dens$Trans,FUN = "mean")

# calculate site means
pol.dens.site <- aggregate(pol.dens.plot$`pol.dens$PctCover`~pol.dens.plot$`pol.dens$Site`,FUN = "mean")
pol.dens.site$se <- aggregate(pol.dens.plot$`pol.dens$PctCover`~pol.dens.plot$`pol.dens$Site`,FUN = "se")[,2]
# append data and clean up workspace


### end stand data calculations

# note this is manually collated from Polaris Project Data on ADC and unpublished data from H. Alexander
tree <- read.csv("C:/Users/mloranty/Documents/GitHub/kolyma_albedo/kolyma_stand_data_FINAL.csv", header = T)


# create spatial points
tree.sf <- st_as_sf(tree, coords = c("Long", "Lat"), crs = 4326) 

# read in biomass data from Berner et al 2012
agb <- raster("L:/data_repo/gis_data/Berner_2011_Kolyma_fire_biomass/kolyma_landsat5_larch_AGB_gm2_2007.tif")

# read in fire shape files from Berner et al 2012
fire <- st_read("L:/data_repo/gis_data/Berner_2011_Kolyma_fire_biomass/Kolyma_FireScars_LS1972-2011.shp")


