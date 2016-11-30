# RUN BLUE 259 VIA HMMWOA
library(hmmwoa)

# SETWD
setwd('~/Documents/WHOI/Data/Blues/2015/141259/') 

#----------------------------------------------------------------------------------#
# ADD MAP DATA
library(rworldmap)
data("countriesLow")

#----------------------------------------------------------------------------------#
# READ IN TAG DATA
ptt <- 141259

# TAGGING LOCATION
iniloc <- data.frame(matrix(c(13, 10, 2015, 41.3, -69.27, 
                              10, 4, 2016, 40.251, -36.061), nrow = 2, ncol = 5, byrow = T))
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

# SPATIAL LIMITS
sp.lim <- list(lonmin = -82, lonmax = -25, latmin = 15, latmax = 50)

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
#get.env(sst.udates[1], type = 'sst', spatLim = sp.lim, save.dir = sst.dir)

# IF USING OHC, DOWNLOAD HYCOM DATA
hycom.dir <- paste('~/Documents/WHOI/RData/HYCOM/', ptt, '/',sep = '')
#get.env(pdt.udates, type = 'ohc', spatLim = sp.lim, save.dir = hycom.dir)

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
#L.ohc <- calc.ohc(pdt, ohc.dir = hycom.dir, dateVec = dateVec, isotherm = '')

#-------

#----------------------------------------------------------------------------------#
# SETUP A COMMON GRID
#----------------------------------------------------------------------------------#

#L.rasters <- list(L.ohc = L.ohc, L.sst = L.sst, L.pdt = L.prof, L.light = L.light)
L.rasters <- list(L.sst = L.sst, L.light = L.light)
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
            L2 = L.res[[1]]$L.light,
            L.mle.res = L.mle.res, dateVec = dateVec,
            locs.grid = locs.grid, iniloc = iniloc)
L.mle <- L$L.mle; L <- L$L

#----------------------------------------------------------------------------------#
# TRY THE MLE.

# NOT RIGHT NOW

#----------------------------------------------------------------------------------#
# OR... JUST DEFINE THE PARAMETERS
par0=c(8.908,10.27,1.152,0.0472,0.707,0.866)
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
# PLOT IT IF YOU WANT TO SEE LIMITS (CI)
#res = apply(f$phi[1,,,],2:3,sum, na.rm=T)
#fields::image.plot(lon, lat, res/max(res), zlim = c(.05,1))

#----------------------------------------------------------------------------------#
# RUN THE SMOOTHING STEP
s = hmm.smoother(f, K1, K2, P, plot = F)

# PLOT IT IF YOU WANT TO SEE LIMITS (CI)
#sres = apply(s[1,,,], 2:3, sum, na.rm=T)
#fields::image.plot(lon, lat, sres/max(sres), zlim = c(.05,1))

#----------------------------------------------------------------------------------#
# GET THE MOST PROBABLE TRACK
#----------------------------------------------------------------------------------#
distr = s
T <- dim(distr)[2]
meanlat <- apply(apply(distr, c(2, 4), sum) * repmat(t(as.matrix(g$lat[,1])), T, 1), 1, sum)
meanlon <- apply(apply(distr, c(2, 3), sum) * repmat(t(as.matrix(g$lon[1,])), T, 1), 1, sum)

plot(meanlon,meanlat,type='l')
world(add=T, fill=T, col='grey')

#=======================================================================================#
## END
#=======================================================================================#


