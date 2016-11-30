# RUN BLUE SHARK EXAMPLE, 141256, USING HMMWOA
library(hmmwoa)

# SETWD
setwd('~/Documents/WHOI/Data/Makos/2015/141257/') 

# READ IN TAG DATA
ptt <- '141257'
# spot is 141267

# TAG/POPUP DATES AND LOCATIONS (dd, mm, YYYY, lat, lon)
iniloc <- data.frame(matrix(c(15, 10, 2015, 41.637, -69.706, 
                              12, 4, 2016, 37.75090,	-71.64880), nrow = 2, ncol = 5, byrow = T))
colnames(iniloc) = list('day','month','year','lat','lon')
tag <- as.POSIXct(paste(iniloc[1,1], '/', iniloc[1,2], '/', iniloc[1,3], sep=''), format = '%d/%m/%Y', tz='UTC')
pop <- as.POSIXct(paste(iniloc[2,1], '/', iniloc[2,2], '/', iniloc[2,3], sep=''), format = '%d/%m/%Y', tz='UTC')
# VECTOR OF DATES FROM DATA. THIS WILL BE THE TIME STEPS, T, IN THE LIKELIHOODS
dateVec <- as.Date(seq(tag, pop, by = 'day')) 

# READ IN DATA FROM WC FILES
#myDir <- '~/Documents/WHOI/RCode/hmmwoa/inst/extdata/' # WHERE YOUR DATA LIVES, THIS IS THE EXAMPLE DATA
myDir <- getwd()
tag.sst <- read.wc(ptt, wd = myDir, type = 'sst', tag=tag, pop=pop); 
sst.udates <- tag.sst$udates; tag.sst <- tag.sst$data

pdt <- read.wc(ptt, wd = myDir, type = 'pdt', tag=tag, pop=pop); 
pdt.udates <- pdt$udates; pdt <- pdt$data

light <- read.wc(ptt, wd = myDir, type = 'light', tag=tag, pop=pop); 
light.udates <- light$udates; light <- light$data

#----------------------------------------------------------------------------------#
# FURTHER PREPARATION
# Set spatial limits and download env data
#----------------------------------------------------------------------------------#

# SET SPATIAL LIMITS, IF DESIRED
sp.lim <- list(lonmin = -85, lonmax = -50, latmin = 20, latmax = 50)

#**
## THE LOCS PART BELOW NO LONGER RUNS BC WE DONT USE 'LOCS' FILE
#**
if (exists('sp.lim')){
  locs.grid <- setup.locs.grid(sp.lim)
} else{
  locs.grid <- setup.locs.grid(locs)
  sp.lim <- list(lonmin = min(locs.grid$lon[1,]), lonmax = max(locs.grid$lon[1,]),
                 latmin = min(locs.grid$lat[,1]), latmax = max(locs.grid$lat[,1]))
}

# IF USING SST, DOWNLOAD THE SST DATA:
sst.dir <- paste('~/Documents/WHOI/RData/SST/OI/', ptt, '/',sep = '')
get.env(sst.udates[1], type = 'sst', spatLim = sp.lim, save.dir = sst.dir)

# IF USING OHC, DOWNLOAD HYCOM DATA
hycom.dir <- paste('~/Documents/WHOI/RData/HYCOM/', ptt, '/',sep = '')
get.env(pdt.udates, type = 'ohc', spatLim = sp.lim, save.dir = hycom.dir)

#----------------------------------------------------------------------------------#
# CALC LIKELIHOODS
#----------------------------------------------------------------------------------#

# LIGHT LIKELIHOOD
L.light <- calc.light(light, locs.grid = locs.grid, dateVec = dateVec)

#-------
# GENERATE DAILY SST LIKELIHOODS
L.sst <- calc.sst(tag.sst, sst.dir = sst.dir, dateVec = dateVec)

#-------
# GENERATE DAILY OHC LIKELIHOODS
L.ohc <- calc.ohc(pdt, ohc.dir = hycom.dir, dateVec = dateVec, isotherm = '')

#-------
# GENERATE DAILY PROFILE LIKELIHOODS
#L.prof <- calc.profile(pdt, dateVec = dateVec, envType = 'hycom', hycom.dir = hycom.dir)
#L.prof.woa <- calc.profile(pdt, dat = woa, lat = lat, lon = lon, dateVec = dateVec, envType = 'woa')

#-------
# GENERATE PROFILE/WOA LIKELIHOODS
# for your data, you will need the dir and extract.woa portions below
# for the example, just use data(woa.quarter) which has already been cropped for the example

# WOA DIRECTORY: LOCATION OF .NC FILE FOR EXTRACTION
#woa.dir <- 'your WOA file.nc'
data(woa.quarter)

# GET THE WOA SUBSET BASED ON SPATIAL LIMITS
#return.woa <- extract.woa(woa.dir, sp.lim, resolution = 'one')
woa <- woa.quarter$watertemp
lon <- woa.quarter$lon 
lat <- woa.quarter$lat
depth <- woa.quarter$depth

# ELIMINATE PACIFIC FROM WOA DATA IF YOU'RE WORKING IN THE W ATLANTIC
#woa <- removePacific(woa, lat, lon)

L.pdt <- calc.pdt.int(pdt, dat = woa, lat = lat, lon = lon, depth = depth, dateVec = dateVec)

#----------------------------------------------------------------------------------#
# SETUP A COMMON GRID
#----------------------------------------------------------------------------------#

#L.rasters <- list(L.ohc = L.ohc, L.sst = L.sst, L.pdt = L.prof, L.light = L.light)
L.rasters <- list(L.sst = L.sst, L.light = L.light, L.ohc = L.ohc)
L.res <- resample.grid(L.rasters, L.rasters$L.sst)
# total of ~5 mins when resampling to ohc, faster when more coarse is desired

L.mle.res <- L.res$L.mle.res
g <- L.res$g; lon <- g$lon[1,]; lat <- g$lat[,1]
g.mle <- L.res$g.mle

#----------------------------------------------------------------------------------#
# LOAD AND FORMAT DATAFRAME OF KNOWN LOCATIONS, IF ANY
#----------------------------------------------------------------------------------#

#colnames(known.locs) <- list('date','lat','lon')
#   where 'date' is from as.Date(known.locs$date)

#----------------------------------------------------------------------------------#
# COMBINE LIKELIHOOD MATRICES
#----------------------------------------------------------------------------------#

L <- make.L(L1 = L.res[[1]]$L.sst,
            L2 = L.res[[1]]$L.ohc,
            L.mle.res = L.mle.res, dateVec = dateVec,
            locs.grid = locs.grid, iniloc = iniloc)
L.mle <- L$L.mle; L <- L$L

Lp <- make.L(L1 = L.res[[1]]$L.pdt, L2 = L.res[[1]]$L.sst,
            L3 = L.res[[1]]$L.light, 
            L.mle.res = L.mle.res, dateVec = dateVec,
            locs.grid = locs.grid, iniloc = iniloc)
Lp.mle <- Lp$L.mle; Lp <- Lp$L

#----------------------------------------------------------------------------------#
# TRY THE MLE.

#t <- Sys.time()
#par0=c(9,10,3,1,0.707,0.866) # from Pedersen 2011
par0=c(8.908,10.27,1.152,0.0472,0.707,0.866)

#fit <- nlm(get.nll.fun, par0, g.mle, L.mle)
#Sys.time() - t

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
K1 = gausskern(D1[1], D1[2], muadv = 0)
K2 = gausskern(D2[1], D2[2], muadv = 0)
P <- matrix(c(p[1],1-p[1],1-p[2],p[2]),2,2,byrow=TRUE)

#----------------------------------------------------------------------------------#
# RUN THE FILTER STEP
f = hmm.filter(g, L, K1, K2, P)
res = apply(f$phi[1,,,],2:3,sum, na.rm=T)
#fields::image.plot(lon, lat, res/max(res), zlim = c(.05,1))
par(mfrow=c(2,2))
image.plot(lon,lat,f$phi[1,10,,])
image.plot(lon,lat,f$phi[2,10,,])
image.plot(lon,lat,f$phi[1,100,,])
image.plot(lon,lat,f$phi[2,100,,])

#----------------------------------------------------------------------------------#
# RUN THE SMOOTHING STEP
s = hmm.smoother(f, K1, K2, P, plot = F)
par(mfrow=c(2,2))
image.plot(lon,lat,s[1,1,,])
world(add=T,fill=T)
points(iniloc[1,c(5:4)])
image.plot(lon,lat,s[2,1,,])
world(add=T,fill=T)
points(iniloc[1,c(5:4)])
image.plot(lon,lat,f$phi[1,1,,])
world(add=T,fill=T)
#points(iniloc[2,c(5:4)])
image.plot(lon,lat,f$phi[2,1,,])
world(add=T,fill=T)
#points(iniloc[2,c(5:4)])

## set some options first
#oopt = ani.options(interval = 0.2, nmax = 3)
## use a loop to create images one by one
pdf('mako257_smooth.pdf',width=12,height=8)
for (i in T:1) {
  image.plot(lon,lat,s[1,i,,])
  world(add=T,fill=T)
  title(paste(i))
  #ani.pause()  ## pause for a while ('interval')
}
dev.off()
## restore the options
#ani.options(oopt)


# PLOT IT IF YOU WANT TO SEE LIMITS (CI)
#sres = apply(s[1,,,], 2:3, sum, na.rm=T)
#fields::image.plot(lon, lat, sres/max(sres), zlim = c(.05,1))

#----------------------------------------------------------------------------------#
# GET THE MOST PROBABLE TRACK
#----------------------------------------------------------------------------------#
distr = s
meanlat <- apply(apply(distr, c(2, 4), sum) * repmat(t(as.matrix(g$lat[,1])), T, 1), 1, sum)
meanlon <- apply(apply(distr, c(2, 3), sum) * repmat(t(as.matrix(g$lon[1,])), T, 1), 1, sum)

plot(meanlon,meanlat,type='l')
world(add=T, fill=T, col='grey')

# ADD THE DATES AND STORE THIS VERSION OF ESTIMATED TRACK
mpt <- cbind(dates = dateVec, lon = meanlon, lat = meanlat)

# PLOT IT!
#data("countriesLow") # ADD MAP DATA
#graphics.off()
spot <- read.table('~/Documents/WHOI/RData/sharkSiteData/AllArgosData.csv', sep=',', header = T)
spot$date <- as.POSIXct(spot$date, format=findDateFormat(spot$date))
spot <- spot[which(spot$ptt == 141267 & spot$date <= as.POSIXct('2016-04-12')),]
plot(spot$lon, spot$lat, ylim=c(35,45),xlim=c(-80,-57), type='l', col='blue')
world(add=T, fill=T)
lines(meanlon, meanlat)

## SOMETHING WEIRD HAPPENING AFTER L IS CREATED. FINAL TRACK DOESNT SHOW TAGGING LOCATION
## AND MAYBE NOT POP UP EITHER? IS 0,0 POSITION COMING FROM T-1 IN FILTER STEP?

#=======================================================================================#
## END
#=======================================================================================#
lon.g <- g.mle$lon[1,]; lat.g <- g.mle$lat[,1]
for (t in 1:181){
  image.plot(lon, lat, s[1,t,,], xlab='', ylab='')
  title(paste(t))
  Sys.sleep(0.5)
}

idx <- c(83,122,124,132,141,145:148,156,157,174)
pdf('check L_v2.pdf', height=14, width=10)
par(mfrow=c(3,1))

for (t in idx){
  #plot(L.res[[1]]$L.light[[t]], xlab='', ylab='')
  #world(add=T)
  #title(paste('light ', t))
  
  plot(L.res[[1]]$L.ohc[[t]], xlab='', ylab='')
  world(add=T)
  title(paste('ohc ', t))
  
  plot(L.res[[1]]$L.sst[[t]], xlab='', ylab='')
  world(add=T)
  title(paste('sst ', t))
  
  image.plot(lon, lat, L[t,,], xlab='', ylab='')
  world(add=T)
  title(paste('L ', t))
  
}

dev.off()

# need to normalize all L rasters before they go into make.L()
# figure out how these all red layers are created and fix it
# how do we get negatives?! (t=147)


