#################################
#
# analyses of post-fire 
# vegetation ~ albedo dynamics
# for Kolyma River basis
#
# MML 07/24/2020
#################################

library(raster)
library(sf)
library(dplyr)


# need to get canopy cover and biomass for each stand to compare with albedo data
# read in stand data
setwd("C:/Users/mloranty/Documents/GitHub/kolyma_albedo/")

dg <- read.csv("data/DensityGradientTrees_raw_2010_2017.csv", header = T)

# calculate cumulative biomass for each plot
dg.bio <- aggregate(dg$Tot.Tree.Carbon.CALC.ALLOM..g.C.m.2.~ dg$Site+dg$Plot,FUN = "sum")

#change header names
names(dg.bio) <- c("site", "plot", "biomassC")

#define function to calculate standard error of the mean
se <- function(x, ...) sqrt(var(x, ...)/length(x))

# calculate mean and SE for each stand based on plot sums
dg.tree <- aggregate(dg.bio$biomassC~dg.bio$site,FUN = "mean")
                 
dg.tree$se <- aggregate(dg.bio$biomassC~ dg.bio$site,FUN = "se")[,2]



# note this is manually collated from Polaris Project Data on ADC and unpublished data from H. Alexander
tree <- read.csv("C:/Users/mloranty/Documents/GitHub/kolyma_albedo/kolyma_stand_data_FINAL.csv", header = T)




# create spatial points
tree.sf <- st_as_sf(tree, coords = c("Long", "Lat"), crs = 4326) 

# read in biomass data from Berner et al 2012
agb <- raster("L:/data_repo/gis_data/Berner_2011_Kolyma_fire_biomass/kolyma_landsat5_larch_AGB_gm2_2007.tif")

# read in fire shape files from Berner et al 2012
fire <- st_read("L:/data_repo/gis_data/Berner_2011_Kolyma_fire_biomass/Kolyma_FireScars_LS1972-2011.shp")


