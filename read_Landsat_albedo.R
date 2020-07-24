##############################
#
# read and analyse Landsat
# albedo & surface temp
# for AGU
#
# MML 11/22/15 
##############################
require(raster)
require(sp)
require(rgdal)
rm(list=ls())

setwd('/Users/mloranty/Google Drive/Documents/Research/kolyma_landsat_albedo/') #from MAC
setwd('/Users/mloranty/Desktop/Google Drive/Documents/Research/kolyma_landsat_albedo/') #from PC - look into changing this. 


## for nw simply read these files and then write them as tiffs, then load into ArcGIS
## or we can try merge!!
## note - cannot merge/mosaic in UTM - need to reproject using Kolyma AEA projection

## read Logan's biomass map, which should be in KAEA
agb <- raster('/Users/mloranty/Google Drive/Documents/Research/GIS_datasets/Berner_2011_data/kolyma_landsat5_larch_AGB_gm2_2007.tif')

#get Kolyma AEA projection as proj4 string
kaea <- projection(agb)
#rcl.val <- c(32000,32767,NA)

############# PROCESS UMB ALBEDO FROM SCHAAF GROUP ################################
row14 <- stack('umb_albedo/lndAlbedo_LC81040142014194LGN00.bin')
#writeRaster(row14,filename='umb_albedo/lndAlbedo_LC81040142014194LGN00.tif',overwrite=T)
NAvalue(row14) <- 32767
#reclassify(row14,rcl.val,filename='umb_albedo/lndAlbedo_LC81040142014194LGN00.tif',overwrite=T)
row14.kaea <- projectRaster(row14,agb,method='bilinear',progress='text',overwrite=T,
                            filename='umb_albedo/lndAlbedo_LC81040142014194LGN00_kaea.tif')

row13 <- stack('umb_albedo/lndAlbedo_LC81040132014194LGN00.bin')
#writeRaster(row13,filename='umb_albedo/lndAlbedo_LC81040132014194LGN00.tif')
NAvalue(row13) <- 32767
row13.kaea <- projectRaster(row13,agb,method='bilinear',progress='text',overwrite=T,
                            filename='umb_albedo/lndAlbedo_LC81040132014194LGN00_kaea.tif')

row12 <- stack('umb_albedo/lndAlbedo_LC81040122014194LGN00.bin')
#reclassify(row12,rcl.val,filename='umb_albedo/lndAlbedo_LC81040122014194LGN00.tif',overwrite=T)
NAvalue(row12) <- 32767
row12.kaea <- projectRaster(row12,agb,method='bilinear',progress='text',overwrite=T,
                            filename='umb_albedo/lndAlbedo_LC81040122014194LGN00_kaea.tif')

row11 <- stack('umb_albedo/lndAlbedo_LC81040112014194LGN00.bin')
NAvalue(row11) <- 32767
row11.kaea <- projectRaster(row11,agb,method='bilinear',progress='text',overwrite=T,
                            filename='umb_albedo/lndAlbedo_LC81040112014194LGN00_kaea.tif')

alb <- mosaic(row14.kaea, row13.kaea, row12.kaea, row11.kaea, fun=mean,progress='text',
              filename='lndAlbedo_LC82014194_Kolyma_merge_kaea.tif',overwrite=T)


## process the Landsat7 files

L7 <- list.files(path='umb_albedo',pattern=glob2rx("lnd*LE7*.bin"),full.name=T)

L7.out <- paste(substr(L7,1,32),"kaea.tif",sep="_")

for(i in 1:length(L7))
{
  r <- stack(L7[i])
  NAvalue(r) <- 32767
  r.kaea <- projectRaster(r,agb,method='bilinear',progress='text',overwrite=T,
                filename=L7.out[i])
}

r11 <- raster('umb_albedo/lndAlbedo_LE710401120_kaea.tif')
r12 <- raster('umb_albedo/lndAlbedo_LE710401220_kaea.tif')
r13 <- raster('umb_albedo/lndAlbedo_LE710401320_kaea.tif')
r14 <- raster('umb_albedo/lndAlbedo_LE710401420_kaea.tif')

alb2 <- mosaic(r14, r13, r12, r11, fun=mean,progress='text',
               filename='lndAlbedo_LE72014234EDC00_Kolyma_merge_kaea.tif',overwrite=T)

############# PROCESS LC8 ALBEDO USING PARAMETRS FROM SCHAAF GROUP ################################
r11 <- raster('LC8_albedo/LC81040112014114LGN00.SW.albedo.laea.tif')
r12 <- raster('LC8_albedo/LC81040122014114LGN00.SW.albedo.laea.tif')
r13 <- raster('LC8_albedo/LC81040132014114LGN00.SW.albedo.laea.tif')
r14 <- raster('LC8_albedo/LC81040142014114LGN00.SW.albedo.laea.tif')

alb <- mosaic(r14, r13, r12, r11, fun=mean,progress='text',
               filename='LC8_albedo/LC82014114.kolyma.SW.albedo.laea.tif',overwrite=T)

r11 <- raster('LC8_albedo/LC81040112014098LGN00.SW.albedo.laea.tif')
r12 <- raster('LC8_albedo/LC81040122014098LGN00.SW.albedo.laea.tif')
r13 <- raster('LC8_albedo/LC81040132014098LGN00.SW.albedo.laea.tif')
r14 <- raster('LC8_albedo/LC81040142014098LGN00.SW.albedo.laea.tif')

alb <- mosaic(r14, r13, r12, r11, fun=mean,progress='text',
              filename='LC8_albedo/LC82014098.kolyma.SW.albedo.laea.tif',overwrite=T)


