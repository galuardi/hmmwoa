# more light likelihood issues

liksrss <- function(obs, srss, srsd){
  # midT = (maxT + minT) / 2
  # Tsd = (maxT - minT) / 4
  #d = obs-srss
  # widx = w >= minT & w <= maxT & !is.na(w)
  sdf = data.frame(sr.gr = as.vector(srss), srsd = as.vector(srsd))
  sdf$srsd[is.na(sdf$srsd)] = 0
  # wint = apply(wdf, 1, function(x) pracma::integral(dnorm, minT, maxT, mean = x[1], sd = x[2]))
  # wint = apply(wdf, 1, function(x) integrate(dnorm, x[1]-x[2], x[1]+x[2], mean = midT, sd = Tsd * 2)$value) 
  res = dnorm(obs, sdf$sr.gr, sdf$srsd)
  srssout = srss
  values(srssout) = res
  # w = w * 0
  # w[widx] = wint
  # w
  srssout
} 

## run this after loading up our package from GitHub:
# RUN BLUE EXAMPLE VIA HMMWOA
library(hmmwoa)

# SETWD
setwd('~/Documents/WHOI/Data/Blues/2015/141256/') 

# READ IN TAG DATA
ptt <- '141256'

# TAG/POPUP DATES AND LOCATIONS (dd, mm, YYYY, lat, lon)
iniloc <- data.frame(matrix(c(13, 10, 2015, 41.575, -69.423, 
                              24, 2, 2016, 26.6798, -69.0147), nrow = 2, ncol = 5, byrow = T))
colnames(iniloc) = list('day','month','year','lat','lon')
tag <- as.POSIXct(paste(iniloc[1,1], '/', iniloc[1,2], '/', iniloc[1,3], sep=''), format = '%d/%m/%Y')
pop <- as.POSIXct(paste(iniloc[2,1], '/', iniloc[2,2], '/', iniloc[2,3], sep=''), format = '%d/%m/%Y')
# VECTOR OF DATES FROM DATA. THIS WILL BE THE TIME STEPS, T, IN THE LIKELIHOODS
dateVec <- as.Date(seq(tag, pop, by = 'day')) 

# READ IN DATA FROM WC FILES
myDir <- '~/Documents/WHOI/RCode/hmmwoa/inst/extdata/' # WHERE YOUR DATA LIVES, THIS IS THE EXAMPLE DATA

light <- read.wc(ptt, wd = myDir, type = 'light', tag=tag, pop=pop); 
light.udates <- light$udates; light <- light$data

#----------------------------------------------------------------------------------#
# LIGHT LIKELIHOOD
# Light-based Longitude Likelihood
#----------------------------------------------------------------------------------#

# SET SPATIAL LIMITS, IF DESIRED
sp.lim <- list(lonmin = -95, lonmax = -52, latmin = 10, latmax = 55)
locs.grid <- setup.locs.grid(sp.lim)

# GET THE LIGHT LIKELIHOOD ELLIPSES
# we'd normally run a calc function here
# L.locs <- calc.locs(locs, gps = NULL, iniloc, locs.grid, dateVec = dateVec, errEll = T)

# INSTEAD HERE ARE THE GUTS OF A CALC FUNCTION FOR LIGHT
#==================
# build the SRSS grids
#==================

# need lat/lon vectors. come from locs.grid??
lon <- locs.grid$lon[1,]
lat <- locs.grid$lat[,1]

# expand.grid and SpatialPoints establishes a grid
xy = as.matrix(expand.grid(lon,lat))
xy = SpatialPoints(xy, proj4string=CRS("+proj=longlat +datum=WGS84"))

# now do the building and rasterize
sr.grid = numeric(length = c(length(lon)*length(lat)*365))
dim(sr.grid) = c(length(lon),length(lat), 365)
ss.grid = sr.grid

fyear = seq(ISOdate(year(dateVec[1]), 1, 1, tz = 'UTC'), ISOdate(year(dateVec[1]), 12, 31, tz = 'UTC'), 'day')
#sr.grid[,,1:365] = sapply(1:365, function(i) matrix(sunriset(xy, fyear[i], direction="sunrise", POSIXct.out=TRUE)$day,length(lon),length(lat)))
#ss.grid[,,1:365] = sapply(1:365, function(i) matrix(sunriset(xy, fyear[i], direction="sunset", POSIXct.out=TRUE)$day,length(lon),length(lat)))
sr.grid[,,289] = matrix(sunriset(xy, fyear[289], direction="sunrise", POSIXct.out=TRUE)$day,length(lon),length(lat))
ss.grid[,,289] = matrix(sunriset(xy, fyear[289], direction="sunset", POSIXct.out=TRUE)$day,length(lon),length(lat))

crs <- "+proj=longlat +datum=WGS84 +ellps=WGS84"
list.ras <- list(x = lon, y = lat, z = sr.grid*24*60)
ex <- raster::extent(list.ras)
sr.ras <- raster::brick(list.ras$z, xmn=ex[1], xmx=ex[2], ymn=ex[3], ymx=ex[4], transpose=T, crs)
sr.ras <- raster::flip(sr.ras, direction = 'y')

list.ras <- list(x = lon, y = lat, z = ss.grid*24*60)
ss.ras <- raster::brick(list.ras$z, xmn=ex[1], xmx=ex[2], ymn=ex[3], ymx=ex[4], transpose=T, crs)
ss.ras <- raster::flip(ss.ras, direction = 'y')

# make some calculations on the tag data: yday, dtime, etc
light <- light[,c('Day','Time','Type')]
light$dtime <- dmy_hms(paste(light$Day, light$Time, sep = ' '))
light$yday <- yday(light$dtime)
light$daymins <- minute(light$dtime) + (hour(light$dtime) * 60)
lightDates <- as.Date(format(light$dtime, '%Y-%m-%d'))

# there will be a loop here:
t = 3
didx <- light$yday[t]
light.t <- light[which(lightDates %in% dateVec[t]),]
# get the SD for this day, T
f1 <- raster::focal(sr.ras[[didx]], w = matrix(1, nrow = 3, ncol = 3), fun = function(x) sd(x, na.rm = T))
# what's the SR?
sr <- light.t$daymins[which(light.t$Type == 'Dawn')]
# the SR likelihood
srlik <- liksrss(sr, srss = sr.ras[[didx]], srsd = f1)

# and sunset
f2 <- raster::focal(ss.ras[[didx]], w = matrix(1, nrow = 3, ncol = 3), fun = function(x) sd(x, na.rm = T))
ss <- light.t$daymins[which(light.t$Type == 'Dusk')]
sslik <- liksrss(ss, srss = ss.ras[[didx]], srsd = f2)

# and you'll see something is very wrong
plot(srlik * sslik) # the overlap is in Brazil
world(add=T)

plot(srlik)
world(add=T)
points(-68.5923,41.3027,pch=16, col='red') # known position for this day

# so what's happening here?
plot(sr.ras[[didx]])
world(add=T)
contour(sr.ras[[didx]],nlevels=30,add=T)
print(paste(sr,'is tag measured SR time. AND IS VERY WRONG. Seems off by one hour.'))

plot(ss.ras[[didx]])
world(add=T)
contour(ss.ras[[didx]],nlevels=30,add=T)
print(paste(ss,'is tag measured SS time. AND IS CORRECT.'))
