##############################
#
# read and calculate Landsat
# albedo from scenes
# for AGU
#
# MML 12/09/15 
##############################
require(raster)
require(sp)
require(rgdal)
rm(list=ls())
 
wd <- ('/Users/mloranty/Google Drive/Documents/Research/kolyma_landsat_albedo/LC8_scenes') #mac
wd <- ('/Users/mloranty/Desktop/Google Drive/Documents/Research/Landsat_albedo/chersky_earth_engine_landsat') #pc
wd <- ('W:/kolyma_landsat_albedo/ee.swath/')
setwd(wd)

## read agb file to use for reprojection ##
agb <- raster('/Users/mloranty/Google Drive/Documents/Research/GIS_datasets/Berner_2011_data/kolyma_landsat5_larch_AGB_gm2_2007.tif')
agb <- raster('C:/Users/mloranty/Desktop/Google Drive/Documents/Research/GIS_datasets/Berner_2011_data/kolyma_landsat5_larch_AGB_gm2_2007.tif')

## make a list of directories to process ##
dirs <- c(list.files(pattern='2014114'),list.files(pattern='2014098'))

## snow albedo function - this is using coefficients from Crystal Schaaf's group ##
snow.alb <- function(a,b,c,d,e,f,na.rm=T){
  return(a*1.2242+b*-0.4318+c*-.3446+d*.3367+e*.1834+f*.2555+.0052)}## snow free albedo function - this is using coefficients from Crystal Schaaf's group ##
snow.free.alb <- function(a,b,c,d,e,f,na.rm=T){
  return(a*.2453+b*.0508+c*.1804+d*.3081+e*.1332+f*.0521+.0011)}
return(a*1.2242+b*-0.4318+c*-.3446+d*.3367+e*.1834+f*.2555+.0052)}
## process snow albedo files - not these were downloaded from earth engine in WGS 84 - maybe not the best solution, but good for now ##
## snow albedo ##
for(i in 1:length(dirs))
{
  # switch to a directory
  setwd(paste(wd,dirs[i],sep="/"))
  
  # read files used to calculate albedo as a stack
  lc8 <- stack(c(list.files(pattern='B2.tif'),list.files(pattern='B3.tif'),
                 list.files(pattern='B4.tif'),list.files(pattern='B5.tif'),
                 list.files(pattern='B6.tif'),list.files(pattern='B7.tif')))
  
  # calculate albedo - files are WGS84
  lc8.alb <- overlay(lc8,fun=snow.alb,unstack=TRUE,progress='text',overwrite=TRUE,
                     filename=paste(dirs[i],'SW.albedo.tif',sep='.'))
  
  # reproject to laea
  lc8.alb.laea <- projectRaster(lc8.alb,agb,method='bilinear',overwrite=TRUE,progress='text',
                                 filename=paste(dirs[i],'SW.albedo.laea.tif',sep='.') )
  
}

## snow free albedo ##
for(i in c(5:6,11:13))
{
  # switch to a directory
  setwd(paste(wd,dirs[i],sep="/"))
  
  # read files used to calculate albedo as a stack
  lc8 <- stack(c(list.files(pattern='.B2.tif'),list.files(pattern='.B3.tif'),
                 list.files(pattern='.B4.tif'),list.files(pattern='.B5.tif'),
                 list.files(pattern='.B6.tif'),list.files(pattern='.B7.tif')))
  
  # calculate albedo - files are WGS84
  lc8.alb <- overlay(lc8,fun=snow.free.alb,unstack=TRUE,progress='text',overwrite=TRUE,
                     filename=paste(dirs[i],'SW.albedo.tif',sep='.'))
  
  # reproject to laea
  lc8.alb.laea <- projectRaster(lc8.alb,agb,method='bilinear',overwrite=TRUE,progress='text',
                                filename=paste(dirs[i],'SW.albedo.laea.tif',sep='.') )
  
}
