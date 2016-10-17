# The below code demonstrates how to use the Continuous Time Correlated Random Walk Model
# of Johnson et al. 2008. Ecology as implemented using the 'crawl' package for various animal movement modeling purposes including: 
# 1. Estimating missing locations with uncertainty- could be used for missing data imputation
# 2. Simulating complete movement paths at finer time scales than data were collected
# 3. Displaying uncertainty in movement paths
# 4. Estimating availability for RSFs

# Clear workspace

#rm(list=ls())
require(sp)
library(crawl)
library(ctmm)

# Load libraries

library(mvtnorm)
library(MASS)
; require(rgdal); require(manipulate);require(survival)

# For use with the crawl package the data frame needs coordinates and a column indicating the relative time
# This relative time column starts at 1 (the first location) and increases by the number of hours between subsequent locations

source('/Users/Cam/Documents/WHOI/RCode/pdtMatch/findDateFormat.r')
setwd('/Users/Cam/Documents/WHOI/Data/Blues/2015/141256/')

# TAG/POPUP DATES AND LOCATIONS (dd, mm, YYYY, lat, lon)
iniloc <- data.frame(matrix(c(13, 10, 2015, 41.575, -69.423, 
                              24, 2, 2016, 26.6798, -69.0147), nrow = 2, ncol = 5, byrow = T))

spot <- read.table('141268-Locations.csv',sep=',',header=T)
spot$t <- as.numeric(as.POSIXct(spot$Date,format=findDateFormat(spot$Date)))


idx = c(which(!duplicated(spot$t)),nrow(spot))  # identifies where each unique date occurs
#values occur and returns vector of those positions
uid = numeric(nrow(spot)) # makes a variable for unique identification
for(i in 1:(length(idx)-1)){
  uid[idx[i]:idx[i+1]] = i #fills uid with unique id #s
}
# was planning to filter exact times according to location class
# but seems that all duplicates are same B class so just thin data
#for(b in 1:nrow(lyd)){
#  lyd[which(uid==b),]
#}

spot <- spot[which(!duplicated(spot$t)),]
#coordinates(spot) <- ~Longitude+Latitude
#proj4string(spot) <- "+proj=longlat +datum=WGS84"
#spot <- spTransform(spot,CRS('+init=epsg:32618'))
#spot$y <- spot@coords[,2]
#spot$x <- spot@coords[,1]
#lyd$t <- as.numeric(lyd$date)
#v1 <- variogram(spot)
#plot.variogram(v1)
#variogram.fit(v1)
#variogram.fit.model
#beta <- variogram.fit.model$tau[1]
#log.beta <- log(beta)
#log.sigma <- (1/2)*log(variogram.fit.model$sigma/variogram.fit.model$tau[1])

spotDate <- as.POSIXct(spot$Date,format='%H:%M:%S %d-%b-%Y')
spot.rt <- as.numeric(spotDate)-as.numeric(spotDate[1])
spot.rt[1] <- 1

spot <- data.frame(relative.time=spot.rt,lon=spot$Longitude,lat=spot$Latitude,
                   date=spotDate, lc=spot$Quality)

# Plot to look at data
plot(spot$lon,spot$lat)

argosClasses <- c("3", "2", "1", "0", "A", "B")
ArgosMultFactors <- data.frame(lc=argosClasses,
                               errX=log(c(1, 1.5, 4, 14, 5.21, 20.78)),
                               errY=log(c(1, 1.5, 4, 14, 11.08, 31.03)))
spotNew <- merge(spot, ArgosMultFactors,
                by=c("lc"), all.x=TRUE)
spotNew <- spotNew[order(spotNew$relative.time),]
spotNew <- spotNew[1:(nrow(spotNew)-2),]

# Check out the parameters 
displayPar(mov.model=~1, err.model=list(x=~errX, y=~errY), drift.model=FALSE,
           data=spotNew, fixPar=c(NA, 1, NA, 1, NA, NA))

#initial.parm <- list(a1.x=c(spotNew$lon[1], 0), a1.y=c(spotNew$lat[1], 0),P1.x=diag(c(.1, .1)),P1.y=diag(c(.1, .1)))
initial.parm <- list(a=c(spotNew$lon[1], 0, spotNew$lat[1], 0), P=diag(c(.1, .1,.1, .1)))

displayPar(mov.model=~1, err.model=list(x=~errX, y=~errY), drift=FALSE,
data=spotNew, coord=c("lon", "lat"), #polar.coord=TRUE,
Time.name="relative.time", initial.state=initial.parm, 
fixPar=c(NA, 1, NA, 1, NA, NA))

fit <- crwMLE(mov.model=~1, err.model=list(x=~errX, y=~errY), drift=FALSE,
              data=spotNew, coord=c("lon", "lat"), #polar.coord=TRUE,
              Time.name="relative.time", initial.state=initial.parm, 
              fixPar=c(NA, 1, NA, 1, NA, NA), #theta=c(-.5,-.5,log.sigma,log.beta),
              control=list(maxit=2000,trace=1, REPORT=10),
              initialSANN=list(maxit=300, trace=1, REPORT=1))

# Now make predictions for hourly locations
# First create a vector of the times (relative to the initial time) at which you want predictions

predTime <- seq(ceiling(min(spotNew$relative.time)), floor(max(spotNew$relative.time)), 60^2*24)
predTime.5 <- seq(ceiling(min(spotNew$relative.time)), floor(max(spotNew$relative.time)), 60^2*24*5)


# Use the crwPredict function to predict locations for every hour
# takes fit argument from model we just fit above
# getuseavail: quantify resource selection

predObj <- crwPredict(object.crwFit=fit, predTime=predTime, speedEst=TRUE, flat=FALSE, getUseAvail=TRUE)
predObj.5 <- crwPredict(object.crwFit=fit, predTime=predTime.5, speedEst=TRUE, flat=FALSE, getUseAvail=TRUE)

# Plot predicted path
# This shows the most probable predicted path (will be straight lines between known locations)
# when you have large gaps track becomes brownian (straight), but when fine-scale
# you see some autocorr in velocity (rounded track)

crwPredictPlot(predObj, "map")

# Predict missing locations- in this case all GPS locations are present, so predicting "missing" hourly locations

# Extract the mean and variance from the predicted objects

mean.X <- predObj$alpha.hat.x
mean.Y <- predObj$alpha.hat.y
var.X <- predObj$V.hat.x
var.Y <- predObj$V.hat.y

# mu.x = position; nu.x = velocity
# var.X also var in position[,1], velocity[,2]

# Plot known locations and fill in with mean predicted locations

crwPredictPlot(predObj, "map")
points(mean.X$mu.x, mean.Y$mu.y, pch=16, col="red")
pred.a=NULL
# sample x,y from normal distribution to look at available locations
for(i in 1:nrow(predObj$useAvail$use)){
  ax <- rnorm(1000, predObj$useAvail$avail[i,2], sqrt(predObj$useAvail$avail[i,4]))
  ay <- rnorm(1000, predObj$useAvail$avail[i,3], sqrt(predObj$useAvail$avail[i,5]))
  pred.a <- rbind(pred.a, cbind(ax,ay,rep(i,1000)))
}


for(i in 2:nrow(spotNew)){
  plot(spotNew$lon, spotNew$lat, pch=16)
  #k <- kde2d(pred.a[pred.a[,3]==(i-1),1],pred.a[pred.a[,3]==(i-1),2], n=100, h=500)
  #image(k, add=T, col=cols)
  points(pred.a[pred.a[,3]==(i-1),1:2])
  #points(RG1.data$X, RG1.data$Y, pch=16)
  #points(RG1.data$X[i-1], RG1.data$Y[i-1], pch=16, col="green")
  #points(RG1.data$X[i], RG1.data$Y[i], pch=16, col="red")
  #	legend(x=727200, y=4420700, pch=c(16,16,1), col=c("red","green","black"), legend=c("Use at time t","Use at time t-1", "Availability"))
  browser()
  # type break to stop browser function
}

pred.a.sp <- data.frame(pred.a)
use <- data.frame(predObj$useAvail$use)
# need data frame containing both used and available
# binary col indicated used (1) or avail (0)
# and id of tstep the row belongs to
for (i in 1:max(pred.a[,3])){
  if(i==1){
    ssf <- cbind(use[i,],used=1,group=i)
  } else{
    ssf <- rbind(ssf,cbind(use[i,],used=1,group=i))
  }
  ssf <- rbind(ssf,cbind(relative.time=use[i,1],meanUse.x=pred.a.sp[which(pred.a.sp[,3]==i),1],
                         meanUse.y=pred.a.sp[which(pred.a.sp[,3]==i),2],
                         varUse.x=NA,varUse.y=NA,used=0,group=i))
}

ssf.coords <- ssf
coordinates(ssf.coords) <- ~meanUse.x+meanUse.y
proj4string(ssf.coords) <- "+proj=longlat +datum=WGS84"

load('lydia_crawl_extract.RData')

envmat <- extract(r,ssf.coords,method='simple')
ssf$sst <- envmat

ssf1 <- ssf[which(!is.na(ssf$sst)),]

m_1 <- clogit(used ~ poly(sst,2) + strata(group), data = ssf1)
m_1 <- clogit(used ~ sst + strata(group), data = ssf1)
summary(m_1)

# also tried thinning data to 5 days instead of 1 but made very little difference
# despite probably getting rid of much of the autocorr in the velocity according
# to the variogram

#' Model selection
#' ==================
#' Model selection is a vast topic. Here just use stepwise backward selection
#' based on AIC
#' comparison based on AIC makes a lot of sense if you think you have some 
#' idea of biologically significant model rather than just throwing
#' a bunch of variables at the model
m_2 <- update(m_1, .~. - slope)
AIC(m_1)
AIC(m_2)
#poly(sst, 2)


coordinates(pred.a.sp) <- ~ax+ay
proj4string(pred.a.sp) <- "+proj=longlat +datum=WGS84"

coordinates(spot) <- ~lon+lat
proj4string(spot) <- "+proj=longlat +datum=WGS84"


source('/Users/Cam/Documents/WHOI/RCode/pdtMatch/extract.woa.r')
nc.dir = '/Users/Cam/Documents/WHOI/RData/pdtMatch/WOA_25deg/global/'
limits = c(-100,-10,-40,60) # (min long, max long, min lat, max lat)
return.woa = extract.woa.crawl(nc.dir, limits, resolution = 'quarter', pdt)
dat = return.woa$dat; lon = return.woa$lon; lat = return.woa$lat; depth = return.woa$depth
image.plot(dat1[,,1]) # plot surface climatology
data <- list(x=lon,y=lat,z=dat1[,,1])
r <- raster(data,crs="+proj=longlat +datum=WGS84")

#sst <- raster(dat1[,,1],xmn=limits[1],xmx=limits[2],ymn=limits[3],ymx=limits[4])



extract.woa.crawl <- 
function(nc.dir,bbox,resolution){
  # Simple function to extract the desired temperature data from a global
  # dataset derived from monthly, 1/4 deg gridded climatology data 
  # contained in the 2013 World Ocean Atlas
  
  # Inputs (both required):
  #   NC.DIR is the directory to load the global nc file from; make sure it's
  #       the only .nc file in the given directory
  #   BBOX is a bounding box of form (long min, long max, lat min, lat max)
  #   RESOLUTION indicates whether oceanographic data is gridded at 'quarter'
  #       or 'one' degree resolution
  #   PDT is output dataframe from convert.pdt
  
  # Output:
  #   RETURNWOA is a list containing:
  #     DAT is an array of temperature data with dimensions (long, lat, depth, time)
  #       depth contains 57 standard depth levels by default and levels are defined
  #       in variable 'depth' contained here. time dimension covers the months of
  #       tag deployment as gathered from querying month data in variable 'pdt'
  #     LON/LAT are vectors of lon/lat bounds
  
  #--------------------------------------------------------------------------
  # - Written by Camrin Braun in R, December 20, 2014. cbraun@whoi.edu
  #--------------------------------------------------------------------------
  
  require(ncdf)
  
  # load global nc
  ncfiles=dir(nc.dir, pattern='.nc')
  nc = open.ncdf(paste(nc.dir,ncfiles,sep=''))
  
  # retrieve var bounds from global nc
  lon = get.var.ncdf(nc, 'Longitude')
  lat = get.var.ncdf(nc, 'Latitude')
  depth = get.var.ncdf(nc, 'Depth')
  
  # set bounds for extracting data
  xmin = which.min((bbox[1]-lon)^2); xmax = which.min((bbox[2]-lon)^2)
  ymin = which.min((bbox[3]-lat)^2); ymax = which.min((bbox[4]-lat)^2)
  
  if(resolution == 'quarter'){
    xlen = 4*(bbox[2]-bbox[1]) # for quarter degree
    ylen = 4*(bbox[4]-bbox[3]) 
  } else if(resolution == 'one'){
    xlen = bbox[2]-bbox[1] # for one degree
    ylen = bbox[4]-bbox[3]
  } else{
    stop('Resolution of input oceanographic data not defined.')
  }  
  
  # define time bounds using tag data
  #tmin = min(pdt$month); tmax = max(pdt$month); tlen = tmax - tmin + 1
  
  #if (tlen <= 1){
  #  tlen = 2
  #}
  
  dat = get.var.ncdf(nc, 'temp', start = c(xmin,ymin,1,tmin), count = c(xlen+1,ylen+1,57,tlen))
  
  returnWOA = list(dat = dat, lon = lon[xmin:xmax], lat = lat[ymin:ymax], depth = depth)
  
}
#fit.err <- crwMLE(mov.model=~1, drift.model=FALSE, err.model=list(x=~1, y=~1), fixPar=c(NA,NA,NA,NA), data=spotd[1:(nrow(spotd)-2),], coord=c("lon", "lat"), 
#                  polar.coord=TRUE, Time.name="relative.time", theta=c(-.5,-.5,log.sigma,log.beta), initial.state=initial.parm, 
#                  control=list(maxit=2000,trace=1, REPORT=10),initialSANN=list(maxit=1000, trace=1, REPORT=10), need.hess=1)


# sigma variance in velocity
# theta autocorr in velocity

displayPar(mov.model=~1, drift.model=FALSE, fixPar=c( NA, NA), data=RG1.data)

# Fit
# error model
# can fit time varying covariates
# allows error estimation with constant error across all locations
# fixPar argument allows fixed (NA) or a constant
# initialSANN - simulated annealing to get initial parameters
# maximum likelihood on model parameters
initial.parm <- list(a1.x=c(spotd$lon[1], 0), a1.y=c(spotd$lat[1], 0),P1.x=diag(c(.1, .1)),P1.y=diag(c(.1, .1)))

fit.err <- crwMLE(mov.model=~1, drift.model=FALSE, err.model=list(x=~1, y=~1), fixPar=c(NA,NA,NA,NA), data=spotd[1:(nrow(spotd)-2),], coord=c("lon", "lat"), 
                  polar.coord=TRUE, Time.name="relative.time", theta=c(-.5,-.5,log.sigma,log.beta), initial.state=initial.parm, 
                  control=list(maxit=2000,trace=1, REPORT=10),initialSANN=list(maxit=1000, trace=1, REPORT=10), need.hess=1)

