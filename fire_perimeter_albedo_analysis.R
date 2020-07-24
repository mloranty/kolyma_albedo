##############################
#
# read Kolyma fire perimeters   
# from Berner et al for  
# Landsat Albedo analysis
#
# MML 12/03/15
#
##############################
# nohup R CMD BATCH ./file.R &
## clear workspace and load packages ##
rm(list=ls())
require(sp)
require(rgdal)
require(raster)

setwd('/Users/mloranty/Google Drive/Documents/Research/kolyma_landsat_albedo/') #mac
setwd('C:/Users/mloranty/Desktop/Google Drive/Documents/Research/kolyma_landsat_albedo/') #pc

## read data from Berner et al 2011 - fire scar perimeters and agb ##
fire <- readOGR('/Users/mloranty/Google Drive/Documents/Research/GIS_datasets/Berner_2011_data/Kolyma_FireScars_LS1972-2011.shp',
                'Kolyma_FireScars_LS1972-2011')
fire <- readOGR('/Users/mloranty/Desktop/Google Drive/Documents/Research/GIS_datasets/Berner_2011_data/Kolyma_FireScars_LS1972-2011.shp',
                'Kolyma_FireScars_LS1972-2011')

agb <- raster('/Users/mloranty/Google Drive/Documents/Research/GIS_datasets/Berner_2011_data/kolyma_landsat5_larch_AGB_gm2_2007.tif')
agb <- raster('/Users/mloranty/Desktop/Google Drive/Documents/Research/GIS_datasets/Berner_2011_data/kolyma_landsat5_larch_AGB_gm2_2007.tif')

## create rasters of fire year, and fire polygon ID for zonal analyses ##
fire.year <- rasterize(fire,agb,'burnYr',filename='fire_data/Kolyma_FireScars_LS1972-2011_year.tif',overwrite=T,progress='text',datatype='INT2U')
fire.id <- rasterize(fire,agb,'Id',filename='fire_data/Kolyma_FireScars_LS1972-2011_year.tif',overwrite=T,progress='text',datatype='INT1U')

## read albedo from Schaaf group @ umb ##
umb.LE7 <- raster('umb_albedo/lndAlbedo_LE72014234EDC00_Kolyma_merge_kaea.tif')
umb.LC8 <- raster('umb_albedo/lndAlbedo_LC82014194_Kolyma_merge_kaea.tif')
## read and mosaic LC8 derived using narrow to boreadband conversion coefficients from Schaaf group
LC8.20140423 <- raster('LC8_albedo/ee.kolyma.20140423.SW.albedo.laea.tif')

## zonal using extract function
smr.alb <- extract(umb.LC8,fire,fun='mean',na.rm=T,progress='text')
smr.alb.sd <- extract(umb.LC8,fire,fun='mean',na.rm=T,progress='text')

spr.alb <- extract(LC8.20140423,fire,fun='mean',na.rm=T,progress='text')
spr.alb.sd <- extract(LC8.20140423,fire,fun='sd',na.rm=T,progress='text')

