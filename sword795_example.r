# RUN SWORDFISH WITH DIFF TIME CONSTRAINTS/DATA FORMATTING
library(hmmwoa)

# SETWD
setwd('~/Documents/WHOI/Data/Swords/2013/106795/') 

#----------------------------------------------------------------------------------#
# ADD MAP DATA
library(rworldmap)
data("countriesLow")

#----------------------------------------------------------------------------------#
# READ IN TAG DATA
ptt <- 106795

# TAGGING LOCATION
iniloc <- data.frame(matrix(c(27, 9, 2013, 46.47683333, -45.5640, 
                              6, 11, 2013, 31.027, -39.684), nrow = 2, ncol = 5, byrow = T))
colnames(iniloc) = list('day','month','year','lat','lon')

# READ IN PDT DATA FROM WC FILES
pdt <- read.table(paste(ptt,'-PDTs.csv', sep=''), sep=',',header=T,blank.lines.skip=F, skip = 2)
pdt <- extract.pdt(pdt)
tag <- as.POSIXct(paste(iniloc[1,1], '/', iniloc[1,2], '/', iniloc[1,3], sep=''), format = '%d/%m/%Y')
pop <- as.POSIXct(paste(iniloc[2,1], '/', iniloc[2,2], '/', iniloc[2,3], sep=''), format = '%d/%m/%Y')
dts <- as.POSIXct(pdt$Date, format = findDateFormat(pdt$Date))
d1 <- as.POSIXct('1900-01-02') - as.POSIXct('1900-01-01')
didx <- dts >= (tag + d1) & dts <= (pop - d1)
pdt <- pdt[didx,]

# VECTOR OF DATES FROM DATA. THIS IS USED IN MANY FUNCTIONS 
udates <- unique(as.Date(pdt$Date))
dateVec <- as.Date(seq(tag, pop, by = 'day'))

#----------------------------------------------------------------------------------#
# LIGHT LIKELIHOOD
# Light-based Longitude Likelihood
#----------------------------------------------------------------------------------#
# READ IN LIGHT DATA FROM WC FILES
locs <- read.table(paste(ptt, '-Locations-GPE2.csv', sep=''), sep=',', header = T, blank.lines.skip = F)
dts <- format(as.POSIXct(locs$Date, format = findDateFormat(locs$Date)), '%Y-%m-%d')
didx <- dts > (tag + d1) & dts < (pop - d1)
locs <- locs[didx,]

# SPATIAL LIMITS
sp.lim <- list(lonmin = -82, lonmax = -25, latmin = 15, latmax = 50)

if (exists('sp.lim')){
  locs.grid <- setup.locs.grid(sp.lim)
} else{
  locs.grid <- setup.locs.grid(locs)
  sp.lim <- list(lonmin = min(locs.grid$lon[1,]), lonmax = max(locs.grid$lon[1,]),
                 latmin = min(locs.grid$lat[,1]), latmax = max(locs.grid$lat[,1]))
}

# GET THE LIKELIHOOD ELLIPSES
t <- Sys.time()
L.locs <- calc.locs(locs, iniloc, locs.grid, dateVec = dateVec, errEll = T)
Sys.time() - t # around 20 seconds with user-defined limits, up to 5 mins or so with locs limits

# something here for quick example plot to check L.locs?

#----------------------------------------------------------------------------------#
# SST LIKELIHOOD
#----------------------------------------------------------------------------------#

# READ IN TAG SST FROM WC FILES
tag.sst <- read.table(paste(ptt, '-SST.csv', sep=''), sep=',',header=T, blank.lines.skip=F)
dts <- as.POSIXct(tag.sst$Date, format = findDateFormat(tag.sst$Date))
didx <- dts >= (tag + d1) & dts <= (pop - d1)
tag.sst <- tag.sst[didx,]
if (length(tag.sst[,1]) <= 1){
  stop('Something wrong with reading and formatting of tags SST data. Check date format.')
}

# IF USING SST
{

dts <- as.POSIXct(tag.sst$Date, format = findDateFormat(tag.sst$Date))
udates <- unique(as.Date(dts))

sst.dir <- paste('~/Documents/WHOI/RData/SST/OI/', ptt, '/',sep = '')

for(i in 1:length(udates)){
  time <- as.Date(udates[i])
  repeat{
    get.oi.sst(sp.lim, time, filename = paste(ptt, '_', time, '.nc', sep = ''), download.file = TRUE, dir = sst.dir) # filenames based on dates from above
    tryCatch({
      err <- try(ncdf::open.ncdf(paste(sst.dir, ptt, '_', time, '.nc', sep = '')), silent = T)
    }, error=function(e){print(paste('ERROR: Download of data at ', time, ' failed. Trying call to server again.', sep = ''))})
    if(class(err) != 'try-error') break
  }
}

t <- Sys.time()
L.sst <- calc.sst(tag.sst, sst.dir = sst.dir, dateVec = dateVec)
Sys.time() - t
# focal for SD calc takes about .5 sec each t step, the dnorm is <.3 sec
# total 2 min for blue shk 259

}

#----------------------------------------------------------------------------------#
# OHC / HYCOM LIKELIHOOD(S)
#----------------------------------------------------------------------------------#

# IF USING OHC HYCOM
{

udates <- unique(as.Date(pdt$Date))
ohc.dir <- paste('~/Documents/WHOI/RData/HYCOM/', ptt, '/',sep = '')

for(i in 1:length(udates)){
  time <- as.Date(udates[i])
  repeat{
    get.hycom(sp.lim, time, type='a', filename = paste(ptt, '_', time, '.nc', sep = ''),
              download.file = TRUE, dir = ohc.dir, vars = 'water_temp') 
    tryCatch({
      err <- try(ncdf::open.ncdf(paste(ohc.dir,ptt,'_',time,'.nc',sep='')),silent=T)
    }, error=function(e){print(paste('ERROR: Download of data at ',time,' failed. Trying call to server again.',sep=''))})
    if(class(err) != 'try-error') break
  }
}

# calc.ohc
t <- Sys.time()
L.ohc <- calc.ohc(pdt, ohc.dir = ohc.dir, dateVec = dateVec, isotherm = '')
Sys.time() - t
# focal takes <8 secs and dnorm 2-7 secs for each t step (day)

}

#----------------------------------------------------------------------------------#
# PDT / WOA LIKELIHOOD
#----------------------------------------------------------------------------------#

# WOA DIRECTORY: LOCATION OF .NC FILE FOR EXTRACTION
# woa.dir = '/Users/Cam/Documents/WHOI/RData/pdtMatch/WOA_25deg/global/woa13_25deg_global_meantemp.nc'
woa.dir = "C:/Users/ben/Documents/WOA/woa13_25deg_global.nc"

# GET THE WOA SUBSET
return.woa = extract.woa(woa.dir, sp.lim, resolution = 'quarter')
woa = return.woa$dat 
lon = as.numeric(return.woa$lon); 
lat = as.numeric(return.woa$lat); 
depth = as.numeric(return.woa$depth)

# ELIMINATE PACIFIC FROM WOA DATA
woa = removePacific(woa, lat, lon)

t <- Sys.time()
L.pdt <- calc.pdt.int(pdt, dat = woa, lat = lat, lon = lon, depth = depth, dateVec = dateVec)
Sys.time() - t

#----------------------------------------------------------------------------------#
# SETUP A COMMON GRID
#----------------------------------------------------------------------------------#

L.rasters <- list(L.pdt = L.pdt, L.ohc = L.ohc, L.locs = L.locs, L.sst = L.sst)
t <- Sys.time()
L.res <- resample.grid(L.rasters, L.rasters$L.ohc)
# total of ~5 mins when resampling to ohc, faster when more coarse is desired
Sys.time() - t
L.mle.res <- L.res$L.mle.res
g <- L.res$g; lon <- g$lon[1,]; lat <- g$lat[,1]
g.mle <- L.res$g.mle

# save.image...
#base::save.image('blue259_example.RData')

#----------------------------------------------------------------------------------#
# MULTIPLY DAILY LIKELIHOOD MATRICES
# check sums of L components
#----------------------------------------------------------------------------------#

# INDEX WHERE LIKELIHOODS ARE ZEROS.. FOR EACH L COMPONENT
L.locs = raster::as.array(L.res[[1]]$L.locs)
L.pdt = raster::as.array(L.res[[1]]$L.pdt)
L.sst = raster::as.array(L.res[[1]]$L.sst)
L.ohc = raster::as.array(L.res[[1]]$L.ohc)

L.locs[is.na(L.locs)] = 0 # turn NA to 0
L.pdt[is.na(L.pdt)] = 0
L.sst[is.na(L.sst)] = 0
L.ohc[is.na(L.ohc)] = 0

# are all cells in a given likelihood surface == 0?
nalocidx = apply(L.locs, 3, sum, na.rm=T) != 0 # does sum of likelihood surface
napdtidx = apply(L.pdt, 3, sum, na.rm=T) != 0
nasstidx = apply(L.sst, 3, sum, na.rm=T) != 0
naohcidx = apply(L.ohc, 3, sum, na.rm=T) != 0

# HERE, WE CHOOSE WHICH L's TO USE
# INDICATES WHICH L LAYERS, IF ANY, ARE ALL ZEROS FOR EACH DAY

# WHERE BOTH ARE ZEROS. THESE WILL BE INTERPOLTED IN THE FILTER
naLidx = nalocidx + nasstidx + naohcidx

# MAKE AN ARRAY OF ZEROS
Lmat = L.pdt * 0
# where naLidx==0, both likelihoods are zero
#       naLidx==1, one has data
#       naLidx==2, both have data
idx1 = naLidx == 1
idx2 = naLidx == 2
idx3 = naLidx == 3

Lmat[,,idx1] = L.locs[,,idx1] + L.sst[,,idx1] + L.ohc[,,idx1] # when only 1 has data
#Lmat[,,idx2] = L.sst[,,idx2] * L.locs[,,idx2] # when both have data
Lmat[,,idx3] = L.locs[,,idx3] * L.sst[,,idx3] * L.ohc[,,idx3] # when all have data

# USE THE INDICES TO POPULATE L
for(b in which(idx2)){
  if(nasstidx[b] & nalocidx[b]){
    Lmat[,,b] = L.sst[,,b] * L.locs[,,b]
  } else if(nasstidx[b] & naohcidx[b]){
    Lmat[,,b] = L.sst[,,b] * L.ohc[,,b]
  } else if(nalocidx[b] & naohcidx[b]){
    Lmat[,,b] = L.locs[,,b] * L.ohc[,,b]
  }
}

# DEFINE A PROJECTION
crs <- "+proj=longlat +datum=WGS84 +ellps=WGS84"

# MAKE A LIST OF LIKELIHOOD
list.pdt <- list(x = lon, y = lat, z = Lmat)
ex <- raster::extent(list.pdt)

# MAKE A RASTER OUT OF IT
T <- dim(Lmat)[3]
for(i in 1:T){
  L.i <- raster::raster(Lmat[,,i], xmn=ex[1], xmx=ex[2], ymn=ex[3], ymx=ex[4], crs)
  if(i==1) L <- L.i else L <- stack(L, L.i)
}

# CREATE A MORE COARSE RASTER FOR PARAMETER ESTIMATION LATER
L.mle <- raster::resample(L, L.mle.res)

# MAKE BOTH RASTERS (COARSE AND FINE RES L's) INTO AN ARRAY
L <- aperm(raster::as.array(raster::flip(L, direction = 'y')), c(3, 2, 1))
L.mle <- aperm(raster::as.array(raster::flip(L.mle, direction = 'y')), c(3, 2, 1))

# CHECK THAT IT WORKED OK
# lon <- g$lon[1,]
# lat <- g$lat[,1]
# fields::image.plot(lon, rev(lat), L[1,,])
# par(mfrow=c(2,1))
# fields::image.plot(lon, lat, L2[1,,])
# image.plot(L1[1,,])
# plot(countriesLow,add=T)

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

#----------------------------------------------------------------------------------#
# MAKE ALL NA'S VERY TINY FOR THE CONVOLUTION
# the previous steps may have taken care of this...
#----------------------------------------------------------------------------------#
L[L == 0] = 1e-15
L[is.na(L)] = 1e-15
L.mle[L.mle == 0] = 1e-15
L.mle[is.na(L.mle)] = 1e-15

#----------------------------------------------------------------------------------#
# TRY THE MLE. SOME OTHER TIME.
{
t <- Sys.time()
par0=c(8.908,10.27,3,1,0.707,0.866) # from Pedersen 2011
fit <- nlm(get.nll.fun, par0, g.mle, L.mle)
Sys.time() - t

## **THESE OUTPUT PARAMETERS ARE PIXEL-BASED. DON'T FORGET TO CONVERT FOR USE
##  WITH THE HIGHER RESOLUTION LIKELIHOOD RESULTS STORED IN L 
D1 <- exp(fit$estimate[1:2]) # parameters for kernel 1. this is behavior mode transit
D2 <- exp(fit$estimate[3:4]) # parameters for kernel 2. resident behavior mode 
p <- 1/(1+exp(-fit$estimate[5:6])) # logit-transformed
#transition probabilities for switching between the two behavioural states 
#Probably need to express kernel movement in terms of pixels per time step.
#The sparse matrix work likely renders this unnecessary, but going back to 
#gausskern, it is. For example, if we have .25 degree and daily time step,
#what would the speed of the fish be when moving fast? 4 pixels/day?

}
#----------------------------------------------------------------------------------#
# OR... JUST DEFINE THE PARAMETERS
D1 <- par0[1:2] # parameters for kernel 1. this is behavior mode transit
D2 <- par0[3:4] # parameters for kernel 2. resident behavior mode
p <- par0[5:6]

#----------------------------------------------------------------------------------#
# GENERATE MOVEMENT KERNELS. D VALUES ARE MEAN AND SD PIXELS
K1 = as.cimg(gausskern(D1[1], D1[2], muadv = 0))
K2 = as.cimg(gausskern(D2[1], D2[2], muadv = 0))
P <- matrix(c(p[1],1-p[1],1-p[2],p[2]),2,2,byrow=TRUE)

#----------------------------------------------------------------------------------#
# RUN THE FILTER STEP
f = hmm.filter(g,L,K1,K2,P)
res = apply(f$phi[1,,,],2:3,sum, na.rm=T)
fields::image.plot(lon, lat, res/max(res), zlim = c(.05,1))

#----------------------------------------------------------------------------------#
# RUN THE SMOOTHING STEP
s = hmm.smoother(f, K1, K2, P, plot = F)

# PLOT IT IF YOU WANT TO SEE LIMITS (CI)
sres = apply(s[1,,,], 2:3, sum, na.rm=T)
fields::image.plot(lon, lat, sres/max(sres), zlim = c(.05,1))

#----------------------------------------------------------------------------------#
# GET THE MOST PROBABLE TRACK
#----------------------------------------------------------------------------------#
distr = s
meanlat <- apply(apply(distr, c(2, 4), sum) * repmat(t(as.matrix(g$lat[,1])), T, 1), 1, sum)
meanlon <- apply(apply(distr, c(2, 3), sum) * repmat(t(as.matrix(g$lon[1,])), T, 1), 1, sum)

# ADD THE DATES AND STORE THIS VERSION OF ESTIMATED TRACK
locs_sst_ohc_par1 <- cbind(dates = dateVec, lon = meanlon, lat = meanlat)

# PLOT IT!
# graphics.off()
# plot(meanlon, meanlat)
# plot(countriesLow, add = T)

#----------------------------------------------------------------------------------#
# COMPARE TO SPOT DATA
#----------------------------------------------------------------------------------#
# READ IN SPOT DATA
# spot = read.csv('C:/Users/ben/Google Drive/Camrin-WOA/hmmwoa_files/121325-SPOT.csv')
spot = read.table('~/Documents/WHOI/RData/sharkSiteData/AllArgosData.csv', sep=',', header = T)
spot <- spot[which(spot$ptt == 141261),]
dts <- as.POSIXct(spot$date, format=findDateFormat(spot$date))
didx <- dts >= tag & dts <= pop
spot <- spot[didx,]

# READ IN GPE3 DATA
gpe <- read.table('~/Documents/WHOI/Data/Blues/2015/141259/141259-6-GPE3.csv',
                  sep=',',header=T, skip = 5)


# PLOT IT
# sres = apply(s,c(3,4), sum, na.rm=T)
# image.plot(lon, lat, sres/max(sres), zlim = c(.01,1),xlim=c(-86,-47),ylim=c(20,45))
plot(meanlon, meanlat, col=2,type='l', xlim=c(-80,-35),ylim=c(20,46))
plot(countriesLow, add = T)
lines(spot$lon, spot$lat)
lines(gpe$Most.Likely.Longitude, gpe$Most.Likely.Latitude, col='orange')
lines(meanlon, meanlat, col=2)

# dist <- as.numeric(unlist(geodetic.distance(cbind(spot[(2:length(spot[,1])),c(8,7)]),cbind(spot[(1:length(spot[,1])-1),c(8,7)]))))
# times <- as.numeric(dts[2:length(dts)] - dts[1:(length(dts)-1)])
# spd <- dist * 1000 / times # m/s


#=======================================================================================#
## END
#=======================================================================================#


