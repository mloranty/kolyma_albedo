##############################
#
# read Siberia field data   
# from Polaris & Fire project  
# for Landsat Albedo analysis
#
# MML 12/07/15
#
##############################
# nohup R CMD BATCH ./file.R &
## clear workspace and load packages ##
rm(list=ls())
require(sp)
require(rgdal)
require(raster)
require(xlsx)

setwd('/Users/mloranty/Google Drive/Documents/Research/manuscripts/Siberia_albedo/')

### NEVERMIND ALL OF THIS FOR NOW - FOR AGU MOST RECENT BOREAL FIRE DATA WAS AGGREGATED MANUALLY WITH POLARIS DATA ###
### READ AS A SINGLE FILE BELOW ###
# ## read TS data - note need to delete 'NA' as xlsx interperets blanks cells as NA ##
# polaris.ts <- read.xlsx(file='/Users/mloranty/Google Drive/Documents/Research/polaris_data/Final_Terrestrial_Data/Polaris_Terr_Survey_Data_2015.xlsx',
#                         sheetIndex=1,header=TRUE,startRow=1,endRow=35,colIndex=1:57)
# ## read boreal data from 2014
# boreal.ts <- read.xlsx(file='field_data//Boreal_TS_Analysis_11.07.14.xlsx',
#                        sheetIndex=1,header=TRUE,startRow=1,endRow=9,colIndex=1:57)
# 
# ## read & aggregate density gradient data from Heather
# boreal.dg <- read.xlsx(file='field_data//DenGradAllSites Mike AGU 2015.xlsx',
#                        sheetIndex=1,header=TRUE)
# boreal.dg <- merge(aggregate(boreal.dg[,3:7],by=list(boreal.dg$Site),FUN='mean',na.rm=T),
#                   aggregate(boreal.dg[,3:7],by=list(boreal.dg$Site),FUN=function(x){return=sd(x,na.rm=T)/sqrt(length(na.omit(x)))}),
#                   by="Group.1",suffixes=c('mean','sd'))
# write.xlsx(boreal.dg,file='field_data/DensityGradient_aggregate_AGU2015.xlsx',col.names=T,row.names=F,sheetName='site_data')

boreal.dat <- read.xlsx(file='field_data//AGU_all_vegC_data_12Dec15.xlsx',
                        sheetIndex=1,header=T,startRow=1,endRow=57,colIndex=1:19)

## convert to spatial points data frame
## note that coords must be ordered long/lat
# polaris.ts.spatial <- SpatialPointsDataFrame(polaris.ts[,6:5],polaris.ts,coords.nrs=6:5,
#                                              proj4string=CRS('+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs'))
# 
# boreal.ts.spatial <- SpatialPointsDataFrame(as.data.frame(boreal.ts[,6:5]),boreal.ts,coords.nrs=6:5,
#                                              proj4string=CRS('+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs'))

boreal.WGS <- SpatialPointsDataFrame(as.data.frame(boreal.dat[,6:5]),boreal.dat,coords.nrs=6:5,
                                     proj4string=CRS('+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs'))


## read in earth engine landsat rasters and reproject the points data frame ##
## SW albedo ##
ee.alb <- list.files(path='/Users/mloranty/Google Drive/Documents/Research/kolyma_landsat_albedo/chersky_earth_engine_landsat/',
                     pattern='SW.albedo.laea.tif',full.names=T,recursive=T)

ee.alb <- stack(ee.alb)

boreal.laea <- spTransform(boreal.WGS,CRS(projection(ee.alb)))
## thermal Band 10 from 2014 ##
ee.2014.b10 <- list.files(path='/Users/mloranty/Google Drive/Documents/Research/kolyma_landsat_albedo/chersky_earth_engine_landsat/',
                          pattern=glob2rx('2014*.B10.tif'),full.names=T,recursive=T)
ee.2014.b10 <- stack(ee.2014.b10)
boreal.utm <- spTransform(boreal.WGS,CRS(projection(ee.2014.b10)))
## thermal Band 11 from 2014 ##
ee.2014.b11 <- list.files(path='/Users/mloranty/Google Drive/Documents/Research/kolyma_landsat_albedo/chersky_earth_engine_landsat/',
                          pattern=glob2rx('2014*.B11.tif'),full.names=T,recursive=T)
ee.2014.b11 <- stack(ee.2014.b11)


## extract the point data ##

stand.alb <- extract(ee.alb,boreal.laea,method='simple')
stand.b10 <- extract(ee.2014.b10,boreal.laea,method='simple')
stand.b11 <- extract(ee.2014.b11,boreal.laea,method='simple')

boreal.all <- cbind(boreal.dat,stand.alb,stand.b10,stand.b11)

plot(boreal.all$X20140423.SW.albedo.laea,boreal.all$Canopy.cover.....Densiometry)
plot(boreal.all$X20140821.SW.albedo.laea,boreal.all$Canopy.cover.....Densiometry)

plot(boreal.all$X20140423.B11,boreal.all$Canopy.cover.....Densiometry)
plot(boreal.all$X20140602.B11,boreal.all$Canopy.cover.....Densiometry)
plot(boreal.all$X20140821.B11,boreal.all$Canopy.cover.....Densiometry)

plot(boreal.all$X20140423.SW.albedo.laea,boreal.all$Larch.Biomass..g.C.m.2.)
plot(boreal.all$X20140821.SW.albedo.laea,boreal.all$Larch.Biomass..g.C.m.2.)

plot(boreal.all$X20140423.B11,boreal.all$Larch.Biomass..g.C.m.2.)
plot(boreal.all$X20140602.B11,boreal.all$Larch.Biomass..g.C.m.2.)
plot(boreal.all$X20140821.B11,boreal.all$Larch.Biomass..g.C.m.2.)

## get just density graident & dg +y4 data ##
boreal.dg <- boreal.all[which(boreal.all$Type=='DG'),]
boreal.ch <- boreal.all[which(boreal.all$Type!='TS'),]

## bin data by biomass
boreal.ch$biomass.bin[which(boreal.ch$Larch.Biomass..g.C.m.2.<750)] <- 'low'
boreal.ch$biomass.bin[which(boreal.ch$Larch.Biomass..g.C.m.2.>=750 & boreal.ch$Larch.Biomass..g.C.m.2.<1500)] <- 'med'
boreal.ch$biomass.bin[which(boreal.ch$Larch.Biomass..g.C.m.2.>=1500)] <- 'high'
recs <- c(7,9,11,12,14,16,18,20:44)

boreal.ch.bin <- aggregate(boreal.ch[,recs],by=list(boreal.ch$biomass.bin),FUN='mean',na.rm=T)
boreal.ch.bin.SE <- aggregate(boreal.ch[,recs],by=list(boreal.ch$biomass.bin),FUN=function(x){return=sd(x,na.rm=T)/sqrt(length(na.omit(x)))})
## MAKE SOME STANDARD RIGAMAROLE ALBEDO PLOTS ##

doy <- as.numeric(round(julian(strptime(substr(colnames(boreal.ch.bin[9:33]),6,9),format='%m%d'),origin="2014-12-31 EST"),digits=0))
y1 <- boreal.ch.bin[,9:14]-boreal.ch.bin.SE[,9:14]
y2 <- boreal.ch.bin[,9:14]+boreal.ch.bin.SE[,9:14]

##ALBEDO TRAJECTORY PLOT##
pdf(file='figures/biomass_albedo_trajectory.pdf',6,6)
par(cex.lab=1.5,cex.axis=1.5)

plot(doy[1:6],boreal.ch.bin[2,9:14],type='o',col='brown',pch=16,lwd=3,cex=1.5,
     xlab='Julian Day',ylab='Albedo',ylim=c(.1,.65),xlim=c(70,250),lty='dashed',main='Seasonal Trajectory')
lines(doy[1:6],boreal.ch.bin[3,9:14],type='o',col='tan',pch=16,lwd=3,cex=1.5,lty='dashed')
lines(doy[1:6],boreal.ch.bin[1,9:14],type='o',col='darkgreen',pch=16,lwd=3,cex=1.5,lty='dashed')

arrows(doy[1:6],as.numeric(y1[2,]),doy[1:6],as.numeric(y2[2,]),col='brown',lwd=1.5,length=0.1,angle=90,code=3)
arrows(doy[1:6],as.numeric(y1[3,]),doy[1:6],as.numeric(y2[3,]),col='tan',lwd=1.5,length=0.1,angle=90,code=3)
arrows(doy[1:6],as.numeric(y1[1,]),doy[1:6],as.numeric(y2[1,]),col='darkgreen',lwd=1.5,length=0.1,angle=90,code=3)
## 2015 data ##
# lines(doy[7:13],boreal.ch.bin[2,15:21],type='o',col='brown',pch=16,lwd=3,cex=1.5,lty='dashed')
# lines(doy[7:13],boreal.ch.bin[3,15:21],type='o',col='tan',pch=16,lwd=3,cex=1.5,lty='dashed')
# lines(doy[7:13],boreal.ch.bin[1,15:21],type='o',col='darkgreen',pch=16,lwd=3,cex=1.5,lty='dashed')

legend('topright',c('Low Biomass','Medium Biomass','High Biomass'),fill=c('brown','tan','darkgreen'),bty='n',inset=0.1,cex=1.25)
dev.off()







##LST PLOT## - not too exciting, maybe look @ residual from local air temp? 
plot(doy[20:25],boreal.ch.bin[2,28:33],type='o',col='brown',pch=16,lwd=3,cex=1.5,
     xlab='Julian Day',ylab='Albedo',lty='dashed')
lines(doy[20:25],boreal.ch.bin[3,28:33],type='o',col='tan',pch=16,lwd=3,cex=1.5,lty='dashed')
lines(doy[20:25],boreal.ch.bin[1,28:33],type='o',col='darkgreen',pch=16,lwd=3,cex=1.5,lty='dashed')

#############################################
### Spring Albedo for individual stands ###
plot.col <- vector(mode='character',length=nrow(boreal.ch))
plot.col[which(boreal.ch$Type=='DG')] <- 'blue'
plot.col[which(boreal.ch$Type!='DG')] <- 'red'

spr.alb <- lm(log(boreal.ch$X20140423.SW.albedo.laea)~boreal.ch$Larch.Biomass..g.C.m.2.)
exp(predict(spr.alb,list(agb))
    
pdf(file='figures/biomass_spring_albedo_.pdf',6,6)
par(cex.lab=1.5,cex.axis=1.5)
plot(boreal.ch$Larch.Biomass..g.C.m.2.,boreal.ch$X20140423.SW.albedo.laea,col=plot.col,pch=16,
     xlab=expression(paste('Larch Aboveground Biomass (g C m'^-2,")",sep="")),ylab='Albedo',
     cex=1.5,ylim=c(.25,.7),xlim=c(0,2500),main='DOY 114')
abline(v=750,lty='dashed')
abline(v=1500,lty='dashed')
arrows(boreal.ch$Larch.Biomass..g.C.m.2.-boreal.ch$SE.2,boreal.ch$X20140423.SW.albedo.laea,
       boreal.ch$Larch.Biomass..g.C.m.2.+boreal.ch$SE.2,boreal.ch$X20140423.SW.albedo.laea,
       col=plot.col,lwd=1.5,length=0.1,angle=90,code=3)

lines(sort(boreal.ch$Larch.Biomass..g.C.m.2.),sort(exp(predict(spr.alb)),decreasing=T),lwd=3,lty='dashed')
text(c(375,1125,1875),0.275,labels=c('Low','Medium','High'),cex=1.25)
text(2000,0.575,expression(R^2 == 0.67))
legend('topright',c('Burn Scar','Y4 Watershed'),pch=16,col=c('blue','red'),bty='n',cex=1.25)
dev.off()

### SUMMER ALBEDO PLOT ###
plot(boreal.ch$Larch.Biomass..g.C.m.2.,boreal.ch$X20140821.SW.albedo.laea,col=boreal.ch$Type,pch=16)

smr.alb <- lm(log(boreal.ch$X20140821.SW.albedo.laea)~boreal.ch$Larch.Biomass..g.C.m.2.)
exp(predict(smr.alb,list(agb))
    
    pdf(file='figures/biomass_summer_albedo_.pdf',6,6)
    par(cex.lab=1.5,cex.axis=1.5)
    plot(boreal.ch$Larch.Biomass..g.C.m.2.,boreal.ch$X20140821.SW.albedo.laea,col=plot.col,pch=16,
         xlab=expression(paste('Larch Aboveground Biomass (g C m'^-2,")",sep="")),ylab='Albedo',
         cex=1.5,ylim=c(.09,.16),xlim=c(0,2500),main='DOY 233')
    abline(v=750,lty='dashed')
    abline(v=1500,lty='dashed')
    arrows(boreal.ch$Larch.Biomass..g.C.m.2.-boreal.ch$SE.2,boreal.ch$X20140821.SW.albedo.laea,
           boreal.ch$Larch.Biomass..g.C.m.2.+boreal.ch$SE.2,boreal.ch$X20140821.SW.albedo.laea,
           col=plot.col,lwd=1.5,length=0.1,angle=90,code=3)
    
    lines(sort(boreal.ch$Larch.Biomass..g.C.m.2.),sort(exp(predict(smr.alb)),decreasing=T),lwd=3,lty='dashed')
    text(c(375,1125,1875),0.09,labels=c('Low','Medium','High'),cex=1.25)
    text(2000,0.13,expression(R^2 == 0.51))
    legend('topright',c('Burn Scar','Y4 Watershed'),pch=16,col=c('blue','red'),bty='n',cex=1.25)
    dev.off()
    
    

#plot(boreal.ch$Larch.Biomass..g.C.m.2.,boreal.ch$X20140423.B11,col=boreal.ch$Type,pch=16)
#plot(boreal.ch$Larch.Biomass..g.C.m.2.,boreal.ch$X20140517.B11,col=boreal.ch$Type,pch=16)

### Biomass ~ Age Plot ###
age.dat <- read.xlsx('/Users/mloranty/Google Drive/Documents/Research/manuscripts/Siberia_albedo/field_data/tree_age_biomass.xlsx',sheetIndex=1,header=T)

plot(age.dat$Age.med,age.dat$Larch.Biomass..g.C.m.2.)

age.dat$age.bin[which(age.dat$Age.mean<50)] <- 'young'
age.dat$bin[which(which(age.dat$Age.mean>=50) & which(age.dat$Age.mean<100))] <- 'mid'
age.dat$bin[which(age.dat$Age.mean>=100)] <- 'old'