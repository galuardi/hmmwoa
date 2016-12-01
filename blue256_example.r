# RUN BLUE SHARK EXAMPLE, 141256, USING HMMWOA
library(hmmwoa)

# SETWD
setwd('~/Documents/WHOI/Data/Blues/2015/141256/') 

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
L.light <- calc.light(light, locs.grid = locs.grid, dateVec = dateVec)

locs <- read.table('141256-Locations-GPE2.csv',sep=',',header=T, blank.lines.skip = F)
L.light.ell <- calc.locs(locs, gps = NULL, iniloc, locs.grid, dateVec, errEll = TRUE, gpeOnly = TRUE)

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
#get.env(pdt.udates, type = 'ohc', spatLim = sp.lim, save.dir = ohc.dir)

# GENERATE DAILY OHC LIKELIHOODS
L.ohc <- calc.ohc(pdt, ohc.dir = ohc.dir, dateVec = dateVec, isotherm = '')

#----------------------------------------------------------------------------------#
# PDT / WOA LIKELIHOOD
#----------------------------------------------------------------------------------#
# FOR YOUR DATA, YOU WILL NEED THE DIR AND EXTRACT.WOA PORTIONS BELOW
# FOR THE EXAMPLE, JUST USE DATA(WOA.QUARTER)

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

# OR JUST USE EXAMPLE DATA FOR NOW
#data(woa.quarter)
#woa <- woa.quarter$watertemp 
#lon <- as.numeric(woa.quarter$lon); 
#lat <- as.numeric(woa.quarter$lat); 
#depth <- as.numeric(woa.quarter$depth)

L.pdt <- calc.pdt.int(pdt, dat = woa, lat = lat, lon = lon, depth = depth, dateVec = dateVec)

#----------------------------------------------------------------------------------#
# SETUP A COMMON GRID
#----------------------------------------------------------------------------------#

L.rasters <- list(L.ohc = L.ohc, L.sst = L.sst, L.light = L.light.ell$L.locs)
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

L <- make.L(L1 = L.res[[1]]$L.ohc, L2 = L.res[[1]]$L.sst, 
            L3 = L.res[[1]]$L.light,
            L.mle.res = L.mle.res, dateVec = dateVec,
            locs.grid = locs.grid, iniloc = iniloc)
L <- make.L(L1 = L.res[[1]]$L.sst, 
            L2 = L.res[[1]]$L.light,
            L.mle.res = L.mle.res, dateVec = dateVec,
            locs.grid = locs.grid, iniloc = iniloc)
L.mle <- L$L.mle; L <- L$L

#----------------------------------------------------------------------------------#
# TRY THE MLE.

t <- Sys.time()
#par0=c(9, 10, 1.152, .0472, 0.707, 0.866) # from Pedersen 2011
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
K1 = gausskern(D1[1], D1[2], muadv = 0)
#K2 = as.cimg(gausskern(D2[1], D2[2], muadv = 0))
K2 = gausskern(D2[1], D2[2], muadv = 0)
P <- matrix(c(p[1],1-p[1],1-p[2],p[2]),2,2,byrow=TRUE)

#----------------------------------------------------------------------------------#
# RUN THE FILTER STEP
f = hmm.filter(g, L, K1, K2, P)
#res = apply(f$phi[1,,,],2:3,sum, na.rm=T)
#fields::image.plot(lon, lat, res/max(res), zlim = c(.05,1))

#----------------------------------------------------------------------------------#
# RUN THE SMOOTHING STEP
s = hmm.smoother(f, K1, K2, P)

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

# ADD THE DATES AND STORE THIS VERSION OF ESTIMATED TRACK
#mpt <- cbind(dates = dateVec, lon = meanlon, lat = meanlat)

# PLOT IT!
#data("countriesLow") # ADD MAP DATA
#graphics.off()
plot(meanlon, meanlat,type='l')
world(add=T)

#plot(countriesLow, add = T)


#=======================================================================================#
## END
#=======================================================================================#



pdf('try L_noknown.pdf', height=8,width=12)
for (i in 1:length(dateVec)){
  image.plot(lon,lat,L[i,,])
  world(add=T)
}
dev.off()


layout(matrix(c(1,1,1,2,3,4), 2, 3, byrow=TRUE))


#L1 = L.res[[1]]$L.ohc , L2 = L.res[[1]]$L.sst, L3 = L.res[[1]]$L.light

# check T = 48:50,54:56
tidx <- c(123,128,129,131,133:135)
pdf('try L_piecewise_1.pdf')#,height=20,width=12)
for(i in tidx){
  #par(mfrow=c(2,1))
  image.plot(lon,lat,L[i,,])
  world(add=T)
  
  par(mfrow=c(1,3))
  plot(L.res[[1]]$L.ohc[[i]])
  world(add=T)
  
  plot(L.res[[1]]$L.sst[[i]])
  world(add=T)
  
  plot(L.res[[1]]$L.light[[i]])
  world(add=T)
  
}
}
dev.off()


r <- raster::resample(L.sst[[123]], L.ohc[[123]], NAflag=NA)

# look at 131
i=131
par(mfrow=c(3,1))
plot(L.res[[1]]$L.ohc[[i]])
world(add=T)

plot(L.res[[1]]$L.sst[[i]])
world(add=T)

plot(L.res[[1]]$L.light[[i]])
world(add=T)


naL1idx = cellStats(L.res[[1]]$L.ohc, sum, na.rm=T) != 0
naL2idx = cellStats(L.res[[1]]$L.sst, sum, na.rm=T) != 0
naL3idx = cellStats(L.res[[1]]$L.light, sum, na.rm=T) != 0

