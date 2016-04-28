
# RUN LYDIA EXAMPLE

# 1) why does the end of the track go to 8 N, -41 E? don't see anything in the
#    likelihood (or in distr near end) that would cause that. L indicates
#    a known pop-up position at T=181 but that isn't reflected after the filter/
#    smooth stages either. The distr array as a likelihood at T=181 that puts
#    the tag in the Caribbean.

library(locfit) # reqd for some stats
library(fields) # plotting only
library(raster) # REQD
library(imager) # reqd
library(ncdf) # reqd
library(dplyr) # reqd
library(rworldmap) # suggested for basemap of plots
library(spatial.tools) # reqd for syncing rasters, better way to do this?

# calculate light-based likelihood
setwd('C:/Users/ben/Google Drive/Camrin-WOA/hmmwoa_files/')
#setwd('~/Documents/WHOI/RData/WhiteSharks/2013/121325/')

source('C:\\Users\\ben\\Documents\\GitHub\\hmmwoa\\R\\findDateFormat.r')
source('C:\\Users\\ben\\Documents\\GitHub\\hmmwoa\\R\\extract.pdt.r')
source('C:\\Users\\ben\\Documents\\GitHub\\hmmwoa\\R\\extract.woa.r')
source('C:\\Users\\ben\\Documents\\GitHub\\hmmwoa\\R\\removePacific.r')
source('C:\\Users\\ben\\Documents\\GitHub\\hmmwoa\\R\\calc.pdt.r')
#source('C:\\Users\\ben\\Documents\\GitHub\\hmmwoa\\R\\plot.woa.r')
source('C:\\Users\\ben\\Documents\\GitHub\\hmmwoa\\R\\calc.locs.r')
source('C:\\Users\\ben\\Documents\\GitHub\\hmmwoa\\R\\misc_funs.r')
source('C:\\Users\\ben\\Documents\\GitHub\\hmmwoa\\R\\hmm2.r')

#---------------------------------------------------------------#
# add map data
#---------------------------------------------------------------#
data("countriesLow")

#---------------------------------------------------------------#
# read in tag data
#---------------------------------------------------------------#
ptt <- 121325

iniloc <- data.frame(matrix(c(3, 3, 2013, 30.3917, -81.3802, 
                              31, 8, 2013, 30.668, -79.972), nrow = 2, ncol = 5, byrow = T))
colnames(iniloc) = list('day','month','year','lat','lon')

pdt <- read.table(paste(ptt,'-PDTs.csv', sep=''), sep=',',header=T,blank.lines.skip=F, skip = 0)

pdt <- extract.pdt(pdt)
tag <- as.POSIXct(paste(iniloc[1,1], '/', iniloc[1,2], '/', iniloc[1,3], sep=''), format = '%d/%m/%Y')
pop <- as.POSIXct(paste(iniloc[2,1], '/', iniloc[2,2], '/', iniloc[2,3], sep=''), format = '%d/%m/%Y')
dts <- as.POSIXct(pdt$Date, format = findDateFormat(pdt$Date))
d1 <- as.POSIXct('1900-01-02') - as.POSIXct('1900-01-01')
didx <- dts >= (tag + d1) & dts <= (pop - d1)
pdt <- pdt[didx,]

lon = c(-90, -40)
lat = c(10, 55)

udates <- unique(as.Date(pdt$Date))
dateVec <- as.Date(seq(tag, pop, by = 'day'))

#---------------------------------------------------------------#
# LIGHT
#---------------------------------------------------------------#

# Light-based Longitude Likelihood (ellipse error is a work in progress)
# do light first so that g is setup for both

locs <- read.table(paste(ptt, '-Locations.csv', sep=''), sep=',', header = T, blank.lines.skip = F)
dts <- format(as.POSIXct(locs$Date, format = findDateFormat(locs$Date)), '%Y-%m-%d')
didx <- dts > (tag + d1) & dts < (pop - d1)
locs <- locs[didx,]

g <- setup.grid(locs, res = 'quarter') # make sure loading function from misc_funs.r
ngrid <- rev(dim(g$lon))
lon <- g$lon[1,]
lat <- g$lat[,1]

L.locs <- calc.locs2(locs, iniloc, g, raster = T, dateVec = dateVec, errEll=T)
L.locs.save <- L.locs

pdf('light ellipses_Lydia.pdf',width=10,height=8)
idx <- which(dateVec %in% as.Date(locs$Date))
for(i in idx){
  plot(L.locs[[i]], main=paste(dateVec[i]))
  points(spot)
}
dev.off()

# try quick plot to check, if raster = 'stack' or 'brick' above
plot(L.locs[[48]])
plot(countriesLow, add = T)

#plot WOA and light likelihood together
# plot(L.locs[[4]]*L.pdt[[4]]) This won't be right beacuse multiplying by zero..
# plot(countriesLow, add = T)

# sync resolutions of pdt to locs to match grid, g
#L.pdt <- spatial_sync_raster(L.pdt, L.locs)

#---------------------------------------------------------------#
# SST
#---------------------------------------------------------------#

tag.sst <- read.table(paste(ptt, '-SST.csv', sep=''), sep=',',header=T, blank.lines.skip=F)
dts <- as.POSIXct(tag.sst$Date, format = findDateFormat(tag.sst$Date))
didx <- dts >= (tag + d1) & dts <= (pop - d1)
tag.sst <- tag.sst[didx,]

if (sst){
  
  dts <- as.POSIXct(tag.sst$Date, format = findDateFormat(tag.sst$Date))
  udates <- unique(as.Date(dts))
  
  sst.dir <- paste('~/Documents/WHOI/RData/SST/OI/', ptt, '/',sep = '')
  lims <- c(min(lon),max(lon),min(lat),max(lat))
  
  for(i in 1:length(udates)){
    time <- as.Date(udates[i])
    repeat{
      get.oi.sst(lims[1:2],lims[3:4],time,filename=paste(ptt,'_',time,'.nc',sep=''),download.file=TRUE,dir=sst.dir) # filenames based on dates from above
      #err <- try(open.ncdf(paste(ohc.dir,ptt,'_',time,'.nc',sep='')),silent=T)
      tryCatch({
        err <- try(open.ncdf(paste(ohc.dir,ptt,'_',time,'.nc',sep='')),silent=T)
      }, error=function(e){print(paste('ERROR: Download of data at ',time,' failed. Trying call to server again.',sep=''))})
      if(class(err) != 'try-error') break
    }
  }
  
  L.sst <- calc.sst(tag.sst, sst.dir = sst.dir, dateVec = dateVec, g=g)
  L.sst.save <- L.sst
}

#---------------------------------------------------------------#
# OHC / HYCOM
#---------------------------------------------------------------#

ohc = FALSE
if (ohc){
  ohc.dir <- paste('~/Documents/WHOI/RData/HYCOM/', ptt, '/',sep = '')
  
  for(i in 1:length(udates)){
    time <- as.Date(udates[i])
    repeat{
      get.hycom(lon,lat,time,type='a',filename=paste('_-',time,'.nc',sep=''),download.file=TRUE,dir=ohc.dir, vars = 'water_temp') # filenames based on dates from above
      #err <- try(open.ncdf(paste(ohc.dir,ptt,'_',time,'.nc',sep='')),silent=T)
      tryCatch({
        err <- try(open.ncdf(paste(ohc.dir,ptt,'_',time,'.nc',sep='')),silent=T)
      }, error=function(e){print(paste('ERROR: Download of data at ',time,' failed. Trying call to server again.',sep=''))})
      if(class(err) != 'try-error') break
    }
  }
  
  # calc.ohc
  L.ohc <- calc.ohc(pdt, ohc.dir = ohc.dir, g=g, dateVec=dateVec, isotherm='', raster = 'stack', downsample = F)
  L.ohc.save <- L.ohc
}


#---------------------------------------------------------------#
# PDT / WOA
#---------------------------------------------------------------#

# set limits of interest
# limits = c(lon, lat) # (min lon, max lon, min lat, max lat)
limits = c(min(lon)-3, max(lon)+3, min(lat)-3, max(lat)+3)

# woa.dir = '/Users/Cam/Documents/WHOI/RData/pdtMatch/WOA_25deg/global/woa13_25deg_global_meantemp.nc'
woa.dir = "C:/Users/ben/Documents/WOA/woa13_25deg_global.nc"

return.woa = extract.woa(woa.dir, limits, resolution = 'quarter')

dat = return.woa$dat 
lon = as.numeric(return.woa$lon); 
lat = as.numeric(return.woa$lat); 
depth = as.numeric(return.woa$depth)

# eliminate Pacific from woa data
dat = removePacific(dat, lat, lon)

# check woa data
graphics.off()
image.plot(lon,lat,dat[,,1,1])
#image.plot(dat$lon,dat$lat,sd[,,1,1])

# perform matching
# 'stack' makes the end of this routine much slower than 'brick' or 'array'
# but is only 10 extra seconds or so

L.pdt <- calc.pdt.int(pdt, dat = dat, lat = lat, lon = lon, g, depth = depth, raster = 'stack', dateVec = dateVec)
L.pdt.save <- L.pdt

# try quick plot to check, if raster = 'stack' or 'brick' above
plot(L.pdt[[10]])
plot(countriesLow, add = T)


#---------------------------------------------------------------#
# multiply daily likelihood matrices

# check sums of L components

# need an index of where likelihoods are zeros.. for each L component

L.locs = as.array(L.locs)
L.pdt = as.array(L.pdt)
L.sst = as.array(L.sst)
L.ohc = as.array(L.ohc)

L.locs[is.na(L.locs)] = 0 # turn NA to 0
L.pdt[is.na(L.pdt)] = 0
L.sst[is.na(L.sst)] = 0
L.ohc[is.na(L.ohc)] = 0

# are all cells in a given likelihood surface == 0?
nalocidx = apply(L.locs,3, sum, na.rm=T)!=0 # does sum of likelihood surface
napdtidx = apply(L.pdt,3, sum, na.rm=T)!=0
nasstidx = apply(L.sst,3, sum, na.rm=T)!=0
naohcidx = apply(L.ohc,3, sum, na.rm=T)!=0

# indicates which L layers, if any, are all zeros for each day
#naLidx = nalocidx+nasstidx # where both are zeros. These will be interpolted in the filter
naLidx = nalocidx+nasstidx+napdtidx
#dateIdx = naLidx==0 # may not need this but here for now..

Lmat = L.pdt*0
# where naLidx==0, both likelihoods are zero
#       naLidx==1, one has data
#       naLidx==2, both have data
idx1 = naLidx==1
idx2 = naLidx==2
idx3 = naLidx==3

Lmat[,,idx1] = L.locs[,,idx1] + L.sst[,,idx1] + L.pdt[,,idx1] # when only 1 has data
#Lmat[,,idx2] = L.sst[,,idx2]*L.locs[,,idx2] # when both have data
Lmat[,,idx3] = L.locs[,,idx3] * L.sst[,,idx3] * L.pdt[,,idx3] # when all have data

for(b in which(idx2)){
  if(nasstidx[b] & nalocidx[b]){
    Lmat[,,b] = L.sst[,,b] * L.locs[,,b]
  } else if(nasstidx[b] & napdtidx[b]){
    Lmat[,,b] = L.sst[,,b] * L.pdt[,,b]
  } else if(nalocidx[b] & napdtidx[b]){
    Lmat[,,b] = L.locs[,,b] * L.pdt[,,b]
  }
}

crs <- "+proj=longlat +datum=WGS84 +ellps=WGS84"
list.pdt <- list(x = lon, y = lat, z = L.pdt)
ex <- extent(list.pdt)

T <- dim(Lmat)[3]
for(i in 1:T){
  L.i <- raster(Lmat[,,i], xmn=ex[1], xmx=ex[2], ymn=ex[3], ymx=ex[4], crs)
  if(i==1) L <- L.i else L <- stack(L, L.i)
}
#ex <- extent(L)
# now need to re-format the array to match dims in sphmm (time, lat, lon)

## bg: you don't need to worry about the time dimension. That can be in either position
L <- aperm(as.array(flip(L, direction = 'y')), c(3,2,1))
# check that it worked ok
#lon <- seq(ex[1], ex[2], length=dim(L)[2])
lon <- g$lon[1,]
lat <- g$lat[,1]
#lat <- seq(ex[3], ex[4], length=dim(L)[3])
image.plot(lon, lat, L[2,,])
plot(countriesLow,add=T)

## ******
## now insert portions of imager script

# from Jonsen et al 2013 - Supplement:
# The input parameter par0 contains the six parameters of the spatial HMM. 
# par0[1:2] are log-transformed movement parameters (diffusivities) pertaining
# to the first behavioural state. par0[3:4] are log-transformed movement 
# parameters (diffusivities) pertaining to the second behavioural state. 
# par0[5:6] are logit-transformed transition probabilities for switching 
# between the two behavioural states. It is not recommended to manually 
# input different parameter values via the par0 input since these values
# are transformed. Instead when analysing new data one should find optimal
# parameters by setting do.fit=TRUE.

# try the NLL
nllf <- get.nll.fun(parvec=c(D1, D2, p), g, L)

#=========
## When using fixed parameters...
par0=c(8.908,10.27,3,1,0.707,0.866) # what units are these?
#par0 = c(40, 10, 10, 5, .707, .866)

D1 <- par0[1:2] # parameters for kernel 1. this is behavior mode transit
D2 <- par0[3:4] # parameters for kernel 2. resident behavior mode
p <- par0[5:6] # logit-transformed transition probabilities for switching between the two behavioural states

# Probably need to express kernel movement in terms of pixels per time step. The sparse matrix work likely renders this unnecessary, but going back to gausskern, it is. For example, if we have .25 degree and daily time step, what would the speed of the fish be when moving fast? 4 pixels/day?

K1 = as.cimg(gausskern(D1[1], D1[2], muadv = 0))
K2 = as.cimg(gausskern(D2[1], D2[2], muadv = 0))
P <- matrix(c(p[1],1-p[1],1-p[2],p[2]),2,2,byrow=TRUE)

# make all NA's very tiny for the convolution
# the previous steps may have taken care of this...
L[L==0] = 1e-15
L[is.na(L)] = 1e-15

f = hmm.filter2(g,L,K1,K2,P)
res = apply(f$phi[1,,,],2:3,sum, na.rm=T)
image.plot(lon, lat, res/max(res), zlim = c(.05,1))

# smooth
# way faster w/o plotting
s = hmm.smoother2(f, K1, K2, P, plot = F)

#=========
## resume here after NLL?

sres = apply(s[1,,,],2:3,sum, na.rm=T)
image.plot(lon, lat, sres/max(sres), zlim = c(.05,1))

sres = apply(s[2,,,],2:3,sum, na.rm=T)
image.plot(lon, lat, sres/max(sres), zlim = c(.05,1))

distr = s
meanlat <- apply(apply(distr,c(2,4),sum)*repmat(t(as.matrix(g$lat[,1])),T,1),1,sum)
meanlon <- apply(apply(distr,c(2,3),sum)*repmat(t(as.matrix(g$lon[1,])),T,1),1,sum)

locs_sst_pdt_par2 <- cbind(dates = dateVec, lon = meanlon, lat = meanlat)

graphics.off()
plot(meanlon, meanlat)
plot(countriesLow, add = T)

#sr = raster(sres/max(sres),xmn = min(lon), xmx = max(lon), ymn = min(lat), ymx = max(lat))
spot = read.csv('C:/Users/ben/Google Drive/Camrin-WOA/hmmwoa_files/121325-SPOT.csv')
#spot = read.csv('~/Documents/WHOI/RData/WhiteSharks/2013/121325/121325-SPOT.csv')
dts <- as.POSIXct(spot$Date, format=findDateFormat(spot$Date))
didx <- dts >= tag & dts <= pop
spot <- spot[didx,]

sres = apply(s,c(3,4), sum, na.rm=T)
image.plot(lon, lat, sres/max(sres), zlim = c(.01,1),xlim=c(-86,-47),ylim=c(20,45))
plot(locs_sst_ohc_par1[,2], locs_sst_ohc_par1[,3], col=2,type='l')
plot(countriesLow, add = T)
lines(spot$Longitude, spot$Latitude)

dist <- as.numeric(unlist(geodetic.distance(cbind(spot[(2:length(spot[,1])),c(8,7)]),cbind(spot[(1:length(spot[,1])-1),c(8,7)]))))
times <- as.numeric(dts[2:length(dts)] - dts[1:(length(dts)-1)])
spd <- dist * 1000 / times # m/s


##########
## END
##########

