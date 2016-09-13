# RUN BLUE 256 VIA HMMWOA
library(hmmwoa)

# SETWD
setwd('~/Documents/WHOI/Data/Blues/2015/141256/') 

#----------------------------------------------------------------------------------#
# ADD MAP DATA
library(rworldmap)
data("countriesLow")

#----------------------------------------------------------------------------------#
# READ IN TAG DATA
ptt <- 141256

# TAG/POPUP LOCATION
iniloc <- data.frame(matrix(c(13, 10, 2015, 41.575, -69.423, 
                              24, 2, 2016, 26.6798, -69.0147), nrow = 2, ncol = 5, byrow = T))
colnames(iniloc) = list('day','month','year','lat','lon')

# READ IN DATA FROM WC FILES
a <- read.wc(ptt, iniloc); attach(a)

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
L.locs <- calc.locs(locs, gps = NULL, iniloc, locs.grid, dateVec = dateVec, errEll = T)

#----------------------------------------------------------------------------------#
# SST LIKELIHOOD
#----------------------------------------------------------------------------------#

# IF USING SST:
{
sst.dir <- paste('~/Documents/WHOI/RData/SST/OI/', ptt, '/',sep = '')

# DOWNLOAD THE SST DATA
get.env(sst.udates[1:3], type = 'sst', spatLim = sp.lim, save.dir = sst.dir)

# GENERATE DAILY SST LIKELIHOODS
L.sst <- calc.sst(tag.sst, sst.dir = sst.dir, dateVec = dateVec)

}

#----------------------------------------------------------------------------------#
# OHC / HYCOM LIKELIHOOD(S)
#----------------------------------------------------------------------------------#

# IF USING OHC HYCOM
{

ohc.dir <- paste('~/Documents/WHOI/RData/HYCOM/', ptt, '/',sep = '')

# DOWNLOAD OHC(HYCOM) DATA
get.env(pdt.udates, type = 'ohc', spatLim = sp.lim, save.dir = ohc.dir)

# GENERATE DAILY OHC LIKELIHOODS
L.ohc <- calc.ohc(pdt, ohc.dir = ohc.dir, dateVec = dateVec, isotherm = '')

}

#----------------------------------------------------------------------------------#
# PDT / WOA LIKELIHOOD
#----------------------------------------------------------------------------------#

# WOA DIRECTORY: LOCATION OF .NC FILE FOR EXTRACTION
 woa.dir <- '/Users/Cam/Documents/WHOI/RData/pdtMatch/WOA_25deg/global/woa13_25deg_global_meantemp.nc'
#woa.dir <- "C:/Users/ben/Documents/WOA/woa13_25deg_global.nc"

# GET THE WOA SUBSET BASED ON SPATIAL LIMITS
return.woa <- extract.woa(woa.dir, sp.lim, resolution = 'quarter')
woa <- return.woa$dat 
lon <- as.numeric(return.woa$lon); 
lat <- as.numeric(return.woa$lat); 
depth <- as.numeric(return.woa$depth)

# ELIMINATE PACIFIC FROM WOA DATA IF YOU'RE WORKING IN THE W ATLANTIC
#woa <- removePacific(woa, lat, lon)

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
# COMBINE LIKELIHOOD MATRICES
#----------------------------------------------------------------------------------#
t <- Sys.time()
L <- make.L(L1 = L.res[[1]]$L.ohc , L2 = L.res[[1]]$L.sst, L3 = L.res[[1]]$L.locs, L.mle.res = L.mle.res)
Sys.time() - t
L.mle <- L$L.mle; L <- L$L

#----------------------------------------------------------------------------------#
# TRY THE MLE. SOME OTHER TIME.

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
locs_sst_ohc_par1 <- cbind(dates = dateVec, lon = meanlon, lat = meanlat)

# PLOT IT!
# graphics.off()
 plot(meanlon, meanlat)
 plot(countriesLow, add = T)

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


