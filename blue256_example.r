# RUN BLUE8 VIA HMMWOA
library(hmmwoa)

# SETWD
setwd('~/Documents/WHOI/Data/Blues/2015/141256/') 
#load('~/Documents/WHOI/RData/Blues/2015/141256/example256_16May.RData')

# READ IN TAG DATA
ptt <- '141256'

# TAG/POPUP DATES AND LOCATIONS (dd, mm, YYYY, lat, lon)
iniloc <- data.frame(matrix(c(13, 10, 2015, 41.575, -69.423, 
                              24, 2, 2016, 26.6798, -69.0147), nrow = 2, ncol = 5, byrow = T))
colnames(iniloc) = list('day','month','year','lat','lon')
tag <- as.POSIXct(paste(iniloc[1,1], '/', iniloc[1,2], '/', iniloc[1,3], sep=''), format = '%d/%m/%Y', tz='UTC')
pop <- as.POSIXct(paste(iniloc[2,1], '/', iniloc[2,2], '/', iniloc[2,3], sep=''), format = '%d/%m/%Y', tz='UTC')
# VECTOR OF DATES FROM DATA. THIS WILL BE THE TIME STEPS, T, IN THE LIKELIHOODS
dateVec <- as.Date(seq(tag, pop, by = 'day')) 

# READ IN DATA FROM WC FILES
myDir <- '~/Documents/WHOI/RCode/hmmwoa/inst/extdata/' # WHERE YOUR DATA LIVES, THIS IS THE EXAMPLE DATA
tag.sst <- read.wc(ptt, wd = myDir, type = 'sst', tag=tag, pop=pop); 
sst.udates <- tag.sst$udates; tag.sst <- tag.sst$data

pdt <- read.wc(ptt, wd = myDir, type = 'pdt', tag=tag, pop=pop); 
pdt.udates <- pdt$udates; pdt <- pdt$data

light <- read.wc(ptt, wd = myDir, type = 'light', tag=tag, pop=pop); 
light.udates <- light$udates; light <- light$data

#----------------------------------------------------------------------------------#
# LIGHT LIKELIHOOD
# Light-based Longitude Likelihood
#----------------------------------------------------------------------------------#

# SET SPATIAL LIMITS, IF DESIRED
sp.lim <- list(lonmin = -95, lonmax = -52, latmin = 10, latmax = 55)

if (exists('sp.lim')){
  locs.grid <- setup.locs.grid(sp.lim)
} else{
  locs.grid <- setup.locs.grid(locs)
  sp.lim <- list(lonmin = min(locs.grid$lon[1,]), lonmax = max(locs.grid$lon[1,]),
                 latmin = min(locs.grid$lat[,1]), latmax = max(locs.grid$lat[,1]))
}

# GET THE LIKELIHOOD ELLIPSES
L.locs <- calc.light(light, locs.grid = locs.grid, dateVec = dateVec)


####
spot <- read.table('121420-Locations.csv', sep=',', header = T)
spot$dtime <- as.POSIXct(spot$Date, format=findDateFormat(spot$Date), tz='UTC')
spot$yday <- yday(spot$dtime)
spot$daymins <- minute(spot$dtime)+hour(spot$dtime)*60
d1 <- as.POSIXct('1900-01-02') - as.POSIXct('1900-01-01')
didx <- spot$dtime >= (tag + d1) & spot$dtime <= (pop - d1)
spot <- spot[didx,]
spot$day <- as.Date(spot$dtime)

pdf('light try.pdf', height=8, width = 12)
for(t in 2:(length(dateVec)-1)){
  plot(L.locs[[t]])
  world(add=T)
  points(spot[which(spot$day %in% dateVec[t]),c(8,7)],pch=16)
}
dev.off()

#----------------------------------------------------------------------------------#
# SST LIKELIHOOD
#----------------------------------------------------------------------------------#

# IF USING SST:
sst.dir <- paste('~/Documents/WHOI/RData/SST/OI/', ptt, '/',sep = '')

# DOWNLOAD THE SST DATA
get.env(sst.udates, type = 'sst', spatLim = sp.lim, save.dir = sst.dir)

# GENERATE DAILY SST LIKELIHOODS
L.sst <- calc.sst(tag.sst, sst.dir = sst.dir, dateVec = dateVec)

#----------------------------------------------------------------------------------#
# OHC / HYCOM LIKELIHOOD(S)
#----------------------------------------------------------------------------------#

# IF USING OHC HYCOM
ohc.dir <- paste('~/Documents/WHOI/RData/HYCOM/', ptt, '/',sep = '')

# DOWNLOAD OHC(HYCOM) DATA
get.env(pdt.udates, type = 'ohc', spatLim = sp.lim, save.dir = ohc.dir)

# GENERATE DAILY OHC LIKELIHOODS
L.ohc <- calc.ohc(pdt, ohc.dir = ohc.dir, dateVec = dateVec, isotherm = '')

#----------------------------------------------------------------------------------#
# PDT / WOA LIKELIHOOD
#----------------------------------------------------------------------------------#
# FOR YOUR DATA, YOU WILL NEED THE DIR AND EXTRACT.WOA PORTIONS BELOW
# FOR THE EXAMPLE, JUST USE DATA(WOA.QUARTER)

# WOA DIRECTORY: LOCATION OF .NC FILE FOR EXTRACTION
#woa.dir <- 'your WOA file.nc'

# GET THE WOA SUBSET BASED ON SPATIAL LIMITS
#return.woa <- extract.woa(woa.dir, sp.lim, resolution = 'one')
#woa <- return.woa$dat 
#lon <- as.numeric(return.woa$lon); 
#lat <- as.numeric(return.woa$lat); 
#depth <- as.numeric(return.woa$depth)

# ELIMINATE PACIFIC FROM WOA DATA IF YOU'RE WORKING IN THE W ATLANTIC
#woa <- removePacific(woa, lat, lon)

# OR JUST USE EXAMPLE DATA FOR NOW
data(woa.quarter)
woa <- woa.quarter$watertemp 
lon <- as.numeric(woa.quarter$lon); 
lat <- as.numeric(woa.quarter$lat); 
depth <- as.numeric(woa.quarter$depth)

L.pdt <- calc.pdt.int(pdt, dat = woa, lat = lat, lon = lon, depth = depth, dateVec = dateVec)

#----------------------------------------------------------------------------------#
# SETUP A COMMON GRID
#----------------------------------------------------------------------------------#

L.rasters <- list(L.pdt = L.pdt, L.ohc = L.ohc, L.locs = L.locs$L.locs, L.sst = L.sst)
L.res <- resample.grid(L.rasters, L.rasters$L.ohc)
# total of ~5 mins when resampling to ohc, faster when more coarse is desired

L.mle.res <- L.res$L.mle.res
g <- L.res$g; lon <- g$lon[1,]; lat <- g$lat[,1]
g.mle <- L.res$g.mle

#----------------------------------------------------------------------------------#
# BUILD ARRAY OF KNOWN LOCATIONS, IF ANY
#----------------------------------------------------------------------------------#

# here we create L.locs with dim(L.ohc) but containing known locations, if any
# known positions may most often be Fastloc/GPS, Argos, acoustic, or sightings data
# that occur between the tag and pop-up times

#----------------------------------------------------------------------------------#
# COMBINE LIKELIHOOD MATRICES
#----------------------------------------------------------------------------------#

L <- make.L(L1 = L.res[[1]]$L.ohc , L2 = L.res[[1]]$L.sst, L3 = L.res[[1]]$L.locs, L.mle.res = L.mle.res)
L.mle <- L$L.mle; L <- L$L

#----------------------------------------------------------------------------------#
# TRY THE MLE.

t <- Sys.time()
par0=c(8.908,10.27,3,1,0.707,0.866) # from Pedersen 2011
fit <- nlm(get.nll.fun, par0, g.mle, L.mle)
Sys.time() - t

## **THESE OUTPUT PARAMETERS ARE PIXEL-BASED. DON'T FORGET TO CONVERT FOR USE
##  WITH THE HIGHER RESOLUTION LIKELIHOOD RESULTS STORED IN L 
D1 <- exp(fit$estimate[1:2]) # parameters for kernel 1. this is behavior mode transit. log-transformed movement parameters (diffusivities) pertaining
# to the first behavioural state.

D2 <- exp(fit$estimate[3:4]) # parameters for kernel 2. resident behavior mode. log-transformed movement 
# parameters (diffusivities) pertaining to the second behavioural state.

p <- 1/(1+exp(-fit$estimate[5:6])) # logit-transformed
#transition probabilities for switching between the two behavioural states 
#Probably need to express kernel movement in terms of pixels per time step.
#The sparse matrix work likely renders this unnecessary, but going back to 
#gausskern, it is. For example, if we have .25 degree and daily time step,
#what would the speed of the fish be when moving fast? 4 pixels/day?

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
mpt <- cbind(dates = dateVec, lon = meanlon, lat = meanlat)

# PLOT IT!
data("countriesLow") # ADD MAP DATA
graphics.off()
plot(meanlon, meanlat)
plot(countriesLow, add = T)


#=======================================================================================#
## END
#=======================================================================================#


