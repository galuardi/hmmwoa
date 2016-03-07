
# RUN BLUE EXAMPLE


library(locfit)
library(fields)
library(raster)
library(imager)
library(ncdf)
library(plyr)
library(rworldmap)
library(spatial.tools)
library(magic)

# calculate light-based likelihood
setwd('C:/Users/ben/Google Drive/Camrin-WOA/hmmwoa_files/')
setwd('~/Documents/WHOI/Data/Blues/2015/141254/')

source('C:\\Users\\ben\\Documents\\GitHub\\hmmwoa\\R\\findDateFormat.r')
source('C:\\Users\\ben\\Documents\\GitHub\\hmmwoa\\R\\extract.pdt.r')
source('C:\\Users\\ben\\Documents\\GitHub\\hmmwoa\\R\\extract.woa.r')
source('C:\\Users\\ben\\Documents\\GitHub\\hmmwoa\\R\\removePacific.r')
source('C:\\Users\\ben\\Documents\\GitHub\\hmmwoa\\R\\calc.pdt.r')
source('C:\\Users\\ben\\Documents\\GitHub\\hmmwoa\\R\\plot.woa.r')
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
# OHC / HYCOM
#---------------------------------------------------------------#

## LET'S IGNORE HYCOM / OHC FOR NOW. CURRENTLY LYDIA'S TIMESPAN ISN'T
## AVAILABLE IN THE UNIFORM SPATIAL PROJECTION. THIS ISN'T A HUGE
## ISSUE AS I'VE MANAGED TO DEAL WITH THAT BUT I'D RATHER GET
## THIS TO YOU NOW AND JUST USE WOA. ONCE THE REST OF THE ROUTINE
## IS WORKING, OHC IS PIECE OF CAKE TO DROP IN.

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
  L.ohc <- calc.ohc(pdt, ohc.dir = ohc.dir, dateVec, isotherm='', raster = 'stack')
  
  #plot.ohc(lik = L.ohc, ohc.dir = ohcdir, pdt = pdt.data, 
  #         filename = paste(ptt,'_ohclik.pdf', sep = ''), write.dir = getwd())
}

#---------------------------------------------------------------#
# Light-based Longitude Likelihood (ellipse error is a work in progress)
# do light first so that g is setup for both
#---------------------------------------------------------------#

locs <- read.table(paste(ptt, '-Locations.csv', sep=''), sep=',', header = T, blank.lines.skip = F)
dts <- format(as.POSIXct(locs$Date, format = findDateFormat(locs$Date)), '%Y-%m-%d')
didx <- dts > (tag + d1) & dts < (pop - d1)
locs <- locs[didx,]

g <- setup.grid(locs, res = 'quarter') # make sure loading function from misc_funs.r
ngrid <- rev(dim(g$lon))
lon <- g$lon[1,]
lat <- g$lat[,1]

L.locs <- calc.locs(locs, iniloc, g, raster = T, dateVec = dateVec)
# try quick plot to check, if raster = 'stack' or 'brick' above
plot(L.locs[[4]])
plot(countriesLow, add = T)

#plot WOA and light likelihood together
# plot(L.locs[[4]]*L.pdt[[4]]) This won't be right beacuse multiplying by zero..
# plot(countriesLow, add = T)

# sync resolutions of pdt to locs to match grid, g
#L.pdt <- spatial_sync_raster(L.pdt, L.locs)


#---------------------------------------------------------------#
# PDT / WOA
#---------------------------------------------------------------#

# set limits of interest
# limits = c(lon, lat) # (min lon, max lon, min lat, max lat)
limits = c(min(lon)-3, max(lon)+3, min(lat)-3, max(lat)+3)

# woa.dir = '/Users/Cam/Documents/WHOI/RData/pdtMatch/WOA_25deg/global/woa13_25deg_global_meantemp.nc'
woa.dir = "C:/Users/ben/Documents/WOA/woa13_25deg_global.nc"
#sd.dir = '/Users/Cam/Documents/WHOI/RData/pdtMatch/WOA_25deg/global/woa13_25deg_global_sd.nc'

return.woa = extract.woa(woa.dir, limits, resolution = 'quarter')
#return.sd = extract.woa(sd.dir, limits, resolution = 'quarter')
#dat = return.woa$dat; 
#lon = as.numeric(return.woa$lon); 
#lat = as.numeric(return.woa$lat); 
#depth = as.numeric(return.woa$depth)

dat = return.woa$dat 
lon = as.numeric(return.woa$lon); 
lat = as.numeric(return.woa$lat); 
depth = as.numeric(return.woa$depth)
#sd = return.sd$dat

# eliminate Pacific from woa data
dat$dat = removePacific(dat$dat, dat$lat, dat$lon)

# check woa data
graphics.off()
image.plot(dat$lon,dat$lat,dat$dat[,,1,1])
#image.plot(dat$lon,dat$lat,sd[,,1,1])

# perform matching
# 'stack' makes the end of this routine much slower than 'brick' or 'array'
# but is only 10 extra seconds or so

### something going wrong in the integration around day 34.. maybe not enough depths?? 
### Also pretty slow... looking into parallelization
#pdt.sub <- pdt[c(1:max(which(as.Date(pdt$Date) %in% dateVec[49]))),]
#dateVec.sub <- dateVec[1:49]
#dat1 <- dat$dat
#pdt.sub <- pdt[1:50,]
#dateVec.sub <- dateVec[1:11]
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
L.locs[is.na(L.locs)] = 0 # turn NA to 0
L.pdt[is.na(L.pdt)] = 0

# are all cells in a given likelihood surface == 0?
nalocidx = apply(L.locs,3, sum, na.rm=T)!=0 # does sum of likelihood surface
napdtidx = apply(L.pdt,3, sum, na.rm=T)!=0

# indicates which L layers, if any, are all zeros for each day
naLidx = nalocidx+napdtidx # where both are zeros. These will be interpolted in the filter
#dateIdx = naLidx==0 # may not need this but here for now..

Lmat = L.pdt*0
# where naLidx==0, both likelihoods are zero
#       naLidx==1, one has data
#       naLidx==2, both have data
idx1 = naLidx==1
idx2 = naLidx==2

Lmat[,,idx1] = L.pdt[,,idx1]+L.locs[,,idx1] # when only 1 has data
Lmat[,,idx2] = L.pdt[,,idx2]*L.locs[,,idx2] # when both have data

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
library(imager)

par0=c(8.908,10.27,3,1,0.707,0.866) # what units are these?
D1 <- par0[1:2] # parameters for kernel 1. this is behavior mode transit
D2 <- par0[3:4] # parameters for kernel 2. resident behavior mode
p <- par0[5:6] # not sure what these parameters are.. look like the diagonal of a 2x2 transition matrix.  


# Probably need to express kernel movement in terms of pixels per time step. The sparse matrix work likely renders this unnecessary, but going back to gausskern, it is. For example, if we have .25 degree and daily time step, what would the speed of the fish be when moving fast? 4 pixels/day?

K1 = as.cimg(gausskern(D1[1], D1[2], muadv = 0))
K2 = as.cimg(gausskern(D2[1], D2[2], muadv = 0))
P <- matrix(c(p[1],1-p[1],1-p[2],p[2]),2,2,byrow=TRUE)

# make all NA's very tiny for the convolution
# the previous steps may have taken care of this...
L[L==0] = 1e-15
L[is.na(L)] = 1e-15

# add a 'skip' index for missing days in the L.. 

# filter - moved function to sphmmfuns_hmm
f = hmm.filter2(g,L,K1,K2,P)
res = apply(f$phi[1,,,],2:3,sum, na.rm=T)
image.plot(lon, lat, res/max(res), zlim = c(.05,1))

# smooth
# way faster w/o plotting
s = hmm.smoother2(f, K1, K2, P, plot = F)

sres = apply(s[1,,,],2:3,sum, na.rm=T)
image.plot(lon, lat, sres/max(sres), zlim = c(.05,1))

sres = apply(s[2,,,],2:3,sum, na.rm=T)
image.plot(lon, lat, sres/max(sres), zlim = c(.05,1))

# calculate track
# don't use this
#calc.track(s, g)  # dimensions flipped...

# switch the dimensions in the calc.track.r... gives a weird output.. ON FIN LAND!
# this is either 1) right and we have to deal with the L and K elements or 2) the dimensions need adjusting..
distr = s
meanlat <- apply(apply(distr,c(2,4),sum)*repmat(t(as.matrix(g$lat[,1])),T,1),1,sum)
meanlon <- apply(apply(distr,c(2,3),sum)*repmat(t(as.matrix(g$lon[1,])),T,1),1,sum)

graphics.off()
plot(meanlon, meanlat)
plot(countriesLow, add = T)

## CDB: not sure what this raster is here but something about its orientation
##      is screwed up. worth having a look at.
#sr = raster(sres/max(sres),xmn = min(lon), xmx = max(lon), ymn = min(lat), ymx = max(lat))
spot = read.csv('C:/Users/ben/Google Drive/Camrin-WOA/hmmwoa_files/121325-SPOT.csv')
#spot = read.csv('~/Documents/WHOI/RData/WhiteSharks/2013/121325/121325-SPOT.csv')
dts <- as.POSIXct(spot$Date, format=findDateFormat(spot$Date))
didx <- dts >= tag & dts <= pop
spot <- spot[didx,]

sres = apply(s,c(3,4), sum, na.rm=T)
image.plot(lon, lat, sres/max(sres), zlim = c(.01,1),xlim=c(-83,-77),ylim=c(27,34))
lines(meanlon, meanlat, pch=19, col=2)
plot(countriesLow, add = T)
lines(spot.sub$Longitude, spot.sub$Latitude, typ='o', pch=19)

#plot(sr)
#plot(countriesLow, add = T)
#lines(meanlon, meanlat, pch=19, col=2)






##########
## END
##########


# try sphmm
## Number of time steps
T <- dim(L)[1]

## Fixed parameter values
par0=c(8.908,10.27,1.152,0.0472,0.707,0.866)
D1 <- par0[1:2]
D2 <- par0[3:4]
p <- par0[5:6]

if(do.fit){
  guess <- c(log(10),log(10),log(0.5),log(0.5),log(0.95/0.05),log(0.95/0.05))
  fit <- nlm(neg.log.lik.fun,guess,g,L,dt)
  D1 <- exp(fit$estimate[1:2])
  D2 <- exp(fit$estimate[3:4])
  p <- 1/(1+exp(-fit$estimate[5:6]))
}

## Setup transition matrices
dt <- 1
G1 <- make.kern(D1,g)
K1 <- uniformization(G1,dt)
#CDB: [1] Error in diag(A) : no method for coercing this S4 class to a vector

G2 <- make.kern(D2,g)
K2 <- uniformization(G2,dt)
P <- matrix(c(p[1],1-p[1],1-p[2],p[2]),2,2,byrow=TRUE)

## Run smoother and filter
f <- hmm.filter(g,L,K1,K2,P)
s <- hmm.smoother(f,K1,K2,P)
sphmm <- calc.track(s,g)
sphmm$date <- lsst$date
sphmm$p.resid <- apply(s,c(1,2),sum)[2,]

# can see known positions from her SPOT tag
tr <- read.table(paste(ptt,'-SPOT.csv', sep=''), sep=',',header=T,blank.lines.skip=F, skip = 0)
tr <- tr[,c(4, 7, 8)]
tr$b <- 1 # set an arbitrary state at this point
colnames(tr) <- list('date','lat','lon','b')
dts <- as.POSIXct(tr$date, format = findDateFormat(tr$date))
didx <- dts >= tag & dts <= pop
tr <- tr[didx,]

plot.results(save.plot = F)


