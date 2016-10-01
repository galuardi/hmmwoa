
library(maptools)
library(lubridate)
library(dplyr)
library(ggplot2)
library(raster)
library(fields)
library(reshape2)

#=========================================================================#
# Functions
#=========================================================================#

## matrix version of sunriset on a grid
get.srss.m<-function(lon,lat, date, direction = c('sunrise','sunset')){
  y1.v<-as.vector(lon)
  y2.v<-as.vector(lat)
  xy = as.matrix(expand.grid(lon,lat))
  srss.v = sapply(1:length(date), function(i) sunriset(xy, date[i], direction=direction, POSIXct.out=TRUE)$day)
  dim(srss.v) <- c(length(lon),length(lat),length(date))
  return(srss.v)
}


## Solar Error
sol.err = function(y0,a0,b0,j,s) {
 (y0^2)/(cos(2*pi*(j+(-1)^(s*b0))/365.25)+a0)
}

## Days since last solstice.. works with non posix dates
day.last.sol <- function(day0){
 dy = date.mdy(day0)$year
 s1 = mdy.date(12,21,dy-1)
 s2 = mdy.date(6,21,dy)
 s3 = mdy.date(12,21,dy)
 s4 = mdy.date(6,21,dy+1)
 svec = c(s1-(day0),s2-(day0),s3-(day0),s4-(day0))
 #min(svec[svec>0])
 svec[min(which(svec>0))]
}

daylen = function(x0,day = alldays[len]){
if(length(x0)>3){
  sr = sunriset(x0, day, direction='sunrise', POSIXct.out=TRUE)$day
  ss = sunriset(x0, day, direction='sunset', POSIXct.out=TRUE)$day
	}else{	 
	 #sr = sapply(1:nrow(x0), function(i) sunriset(matrix(as.numeric(x0[i,]),nrow=1), day[i], direction='sunrise', POSIXct.out=TRUE)$day)
	 #ss = sapply(1:nrow(x0), function(i) sunriset(matrix(as.numeric(x0[i,]),nrow=1), day[i], direction='sunset', POSIXct.out=TRUE)$day)
	 ss = sunriset(matrix(as.numeric(x0),nrow=1), day, direction='sunset', POSIXct.out=TRUE)$day
	 sr = sunriset(matrix(as.numeric(x0),nrow=1), day, direction='sunrise', POSIXct.out=TRUE)$day
	}
 ss-sr
}

# use dnorm to geenrate a likelihood
# use focal to get a sd field
# likint3 <- function(w, wsd, minT, maxT){
liksrss <- function(obs, srss, srsd){
  if(is.na(obs)){
    srssout = srss
    values(srssout) = 0
    
  } else{
    # midT = (maxT + minT) / 2
    # Tsd = (maxT - minT) / 4
    #d = obs-srss
    # widx = w >= minT & w <= maxT & !is.na(w)
    sdf = data.frame(sr = as.vector(srss), srsd = as.vector(srsd))
    sdf$srsd[is.na(sdf$srsd)] = 0
    # wint = apply(wdf, 1, function(x) pracma::integral(dnorm, minT, maxT, mean = x[1], sd = x[2]))
    # wint = apply(wdf, 1, function(x) integrate(dnorm, x[1]-x[2], x[1]+x[2], mean = midT, sd = Tsd * 2)$value) 
    res = dnorm(obs, sdf$sr, sdf$srsd)
    srssout = srss
    values(srssout) = res
    # w = w * 0
    # w[widx] = wint
    # w
  }
  srssout
} 

gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}


#=======================================================#
# Make the Likelihood array for SRSS
#=======================================================#
day0 = as.POSIXct(ymd('2011-1-01'))

#================================================#
# vector of interest 1 degree resolution
#================================================#
dres = 1

lon = seq(-80, -50, dres)
lat = seq(20, 50, dres)
days = day0 +86400*c(0:364)

#================================================#
# test the plotting
#================================================#
xy = as.matrix(expand.grid(lon,lat))
xy = SpatialPoints(xy, proj4string=CRS("+proj=longlat +datum=WGS84"))

#image(lon,lat,matrix(sunriset(xy, day0, direction="sunset", POSIXct.out=TRUE)$day,length(lon),length(lat)))
contour(lon,lat,matrix(sunriset(xy, day0, direction="sunrise", POSIXct.out=TRUE)$day*24*60, length(lon),length(lat)),add=T, nlevels = 50, col = 'lightblue')

contour(lon,lat,matrix(sunriset(xy, day0, direction="sunset", POSIXct.out=TRUE)$day*24*60, length(lon),length(lat)),add=T,nlevels = 50, col = 'salmon')

#================================================#
# Build the 3D grids
#================================================#
# make storage matrices 
sr.grid = numeric(length = c(length(lon)*length(lat)*365))
dim(sr.grid) = c(length(lon),length(lat), 365)
ss.grid = sr.grid

fyear = seq(ISOdate(2015,1,1), ISOdate(2015,12,31), 'day')
t1 = Sys.time()
sr.grid[,,1:365] = sapply(1:365, function(i) matrix(sunriset(xy, fyear[i], direction="sunrise", POSIXct.out=TRUE)$day,length(lon),length(lat)))
t2 = Sys.time()-t1
t2 # less than 2 mins with a 1/4 degree grid over the NWA

ss.grid[,,1:365] = sapply(1:365, function(i) matrix(sunriset(xy, fyear[i], direction="sunset", POSIXct.out=TRUE)$day,length(lon),length(lat)))

# now turn the grids into rasters
# ** this is where the problem was **
# raster needed t() and a flip
crs <- "+proj=longlat +datum=WGS84 +ellps=WGS84"
list.ras <- list(x = lon, y = lat, z = sr.grid*24*60)
ex <- raster::extent(list.ras)
sr.ras <- raster::brick(list.ras$z, xmn=ex[1], xmx=ex[2], ymn=ex[3], ymx=ex[4], transpose=T, crs)
sr.ras <- raster::flip(sr.ras, direction = 'y')

list.ras <- list(x = lon, y = lat, z = ss.grid*24*60)
ss.ras <- raster::brick(list.ras$z, xmn=ex[1], xmx=ex[2], ymn=ex[3], ymx=ex[4], transpose=T, crs)
ss.ras <- raster::flip(ss.ras, direction = 'y')


#r1 = raster::brick(sr.grid[,,didx]*24*60, xmn = min(lon), xmx = max(lon), ymn = min(lat), ymx = max(lat))
#f1 = raster::focal(sr.ras[[didx]], w = matrix(1, nrow = 3, ncol = 3), fun = function(x) sd(x, na.rm = T))

# sr = srss1$daymins[3]

#================================================================================#
# Use the example dataset to build likelihood
#================================================#

#load('C:/Users/benjamin.galuardi/Google Drive/Camrin-WOA/srss_example.RData')
load('~/Documents/WHOI/RData/Blues/2015/141256/srss_example.RData')
srss1 = light[,c('Day','Time','Type')]
srss1$dtime = dmy_hms(paste(srss1$Day, srss1$Time, sep = " "))#-3600*4
srss1$yday = yday(srss1$dtime)
srss1$daymins = minute(srss1$dtime)+hour(srss1$dtime)*60
### adjust for times less than 500 minutes. at GMT-4, we will never have sunrise times of less than 500 (~8:30 GMT or 4:30 EST). We would need to be pretty far west for this to be valid. ###
idx = srss1$daymins<500
srss1$daymins[idx] = srss1$daymins[idx]+500

ggplot(light, aes(x = dtime, y = daymins, colour = Type))+
  geom_point(size = 3)

# now add spot srss time calculations
spot$dtime <- as.POSIXct(spot$date, format='%Y-%m-%d %H:%M:%S', tz='UTC')
spot$yday <- yday(spot$dtime)
spot$daymins <- minute(spot$dtime)+hour(spot$dtime)*60
spot.xy = SpatialPoints(spot[,c(4,3)], proj4string=CRS("+proj=longlat +datum=WGS84"))
#for(i in 1:length(spot[,1])){
#  spot$sst[i] <- sunriset(spot.xy[i], spot$dtime[i], direction='sunset', POSIXct.out=T)$day*24*60
#}
srt <- sapply(1:length(spot[,1]), function(i) sunriset(spot.xy[i], spot$dtime[i], direction='sunrise'))
spot$srt <- srt*24*60
sst <- sapply(1:length(spot[,1]), function(i) sunriset(spot.xy[i], spot$dtime[i], direction='sunset'))
spot$sst <- sst*24*60

# and add extracted values for each spot location
spot$srt.ex <- NA; spot$sst.ex <- NA
for(ii in 1:length(spot[,1])){
  srt.ex <- raster::extract(sr.ras[[spot$yday[ii]]], spot.xy[ii])
  sst.ex <- raster::extract(ss.ras[[spot$yday[ii]]], spot.xy[ii])
  spot$srt.ex[ii] <- srt.ex
  spot$sst.ex[ii] <- sst.ex
}

#spot$srt <- spot$srt + 240
#spot$sst <- spot$sst + 240

# something weird is happening when this is melted and plotted because the original srss calcs are correct
# but when plotted up the line for srt is weird. can be ignored i think considering this is exploratory
# and the problem doesn't appear to be relevant to the likelihood calc
spot.m <- melt(spot,id.vars=c('ptt','date','lat','lon','class','dtime','yday','daymins'),measure.vars=c('srt','sst','srt.ex','sst.ex'))

cols <- gg_color_hue(6)
ggplot(light, aes(x = dtime, y = daymins, colour = Type))+
  geom_point(size = 3) +
  geom_line(data=spot.m, aes(x=dtime, y=value, colour=variable))+
  scale_color_manual(values=c(cols[1:5],'#000000',cols[6],'#000000'))
  
#=============
# now let's try the likelihood
#=============
# known date/time from the first spot position
tdate = as.POSIXct(ymd('2016-09-21', tz = 'UTC'))
didx = yday(tdate)
tloc <- matrix(as.numeric(c(-70.6731, 41.5265)),1,2) #woods hole

tloc <- matrix(as.numeric(c(-64.863, 34.617)),1,2) #Lydia 11 Jun 2013
tdate = as.POSIXct(ymd('2013-06-11', tz = 'UTC'))
didx = yday(tdate)

f1 = raster::focal(sr.ras[[didx]], w = matrix(1, nrow = 3, ncol = 3), fun = function(x) sd(x, na.rm = T))
sr = sunriset(tloc, tdate, direction="sunrise", POSIXct.out=TRUE)$day*60*24
print(paste(sr,' from sunriset')); print(paste(extract(sr.ras[[didx]],tloc),' from raster extract'))
srlik = liksrss(sr, srss = sr.ras[[didx]], srsd = f1)

# and sunset
f2 = raster::focal(ss.ras[[didx]], w = matrix(1, nrow = 3, ncol = 3), fun = function(x) sd(x, na.rm = T))
ss = sunriset(tloc, tdate, direction="sunset", POSIXct.out=TRUE)$day*60*24
print(paste(ss,' from sunriset')); print(paste(extract(ss.ras[[didx]],tloc),' from raster extract'))
sslik = liksrss(ss, srss = ss.ras[[didx]], srsd = f2)

pdf('test_lik.pdf')         
plot((srlik)*(sslik))         
contour(sr.ras[[didx]], add=T, nlevels = 20)
contour(ss.ras[[didx]], add=T, nlevels = 20)
world(add=T)
points(tloc, pch = 19, col = 2)
title('for today in Woods Hole')
dev.off()






#======================================================================================#
# plot(density(srss1$daymins[1]-sr.grid[,,90]*24*60))

mdat = melt(adply(sr.grid[,,1:90]*24*60, 3, function(x) srss1$daymins[1]-x), 'X1')
names(mdat) = c('day','diff','minutes')

ggplot(mdat, aes(x = day, y = minutes, colour = diff))+
  # geom_density
geom_line()


#================================================#
# try along the z axis 
#================================================#

get.srss.m(lon,lat,fyear[1], 'sunrise') 

k=6

#================================================#
#try for all days in a track
#================================================#

day0 <- day0#psat[[k]]$day0b
dayT <- day0+86400*365#psat[[k]]$dayTb
alldays = seq(day0,dayT,'day')
alldays = as.POSIXct(as.numeric(alldays), tz='GMT', origin="1970-01-01")

sr.grid = numeric(length = c(length(lon)*length(lat)*length(alldays)))
dim(sr.grid) = c(length(lon),length(lat), length(alldays))
ss.grid = sr.grid

t1 = Sys.time()
for(i in 1:length(alldays)){
  print(i)
  sr.grid[,,i] = get.srss.m(lon,lat,alldays[i], 'sunrise')*24*60#+4*60 # converts to minutes past midnight GMT
  ss.grid[,,i] = get.srss.m(lon,lat,alldays[i], 'sunset')*24*60#+4*60
}
t2=Sys.time()-t1


par(mfrow=c(2,5))
for( i in 1:10){
  contour(lon,lat,sr.grid[,,i],add=F, nlevels=50, col=rainbow(50))
  contour(lon,lat,ss.grid[,,i],add=T, nlevels=50, col=rainbow(50))
}

x11()
for( i in 1:11){
  if (i==1){ contour(lon,lat,sr.grid[,,i],add=F, nlevels=50, col=i)
  }else{ contour(lon,lat,sr.grid[,,i],add=T, nlevels=50, col=i)}
}
sapply(1:11, function(i) contour(lon,lat,ss.grid[,,i],col=rainbow(50)[i],add=T, nlevels=50))

#================================================#
# Day Length Grid
#================================================#

dl.grid = ss.grid-sr.grid
mid.grid = (ss.grid+sr.grid)/2
par(mfrow=c(1,2))
image.plot(lon,lat,dl.grid[,,1]/60); title('Day length in hours')
image.plot(lon,lat,mid.grid[,,1]/60); title('time of noon GMT')


temp = psat[[k]]$SRSS
len = nrow(temp)
x0 = rev(psat[[k]]$x0)
xT = rev(psat[[k]]$xT)

temp[1,2] = sunriset(matrix(as.numeric(x0),nrow=1), alldays[1], direction='sunrise', POSIXct.out=TRUE)$day*24*60+4*60
temp[1,3] = sunriset(matrix(as.numeric(x0),nrow=1), alldays[1], direction='sunset', POSIXct.out=TRUE)$day*24*60+4*60
temp[len,2] = sunriset(matrix(as.numeric(x0),nrow=1), alldays[len], direction='sunrise', POSIXct.out=TRUE)$day*24*60+4*60
temp[len,3] = sunriset(matrix(as.numeric(x0),nrow=1), alldays[len], direction='sunset', POSIXct.out=TRUE)$day*24*60+4*60

temp.dl = temp[,3]-temp[,2]
temp.md = (temp[,3]+temp[,2])/2

# ind = which.min((temp.dl[1]-dl.grid[,,1]/60)^2)
# txy = ind2sub(c(1:length(lat), 1:lenhe gth(lat)),ind);

image.plot(lon+360,lat,dl.grid[,,1]/60*mid.grid[,,1]/60);
points(x0[1]+360, x0[2], pch=19 ,col=2)
plot(map2$SP,add=T)

# how to define the SRSS likelihood???
#1 use solstice variance
#2 Diffusion in a day :kernel

j= 16				# single day at liberty

image.plot(lon+360,lat,(dl.grid[,,j])-(temp[j,3]-temp[j,2]),zlim=c(-5,5));  # this part should vary w/date
image.plot(lon+360,lat,(mid.grid[,,j])-(temp[j,3]+temp[j,2])/2,zlim=c(-5,5),add=T);  # this part should vary w/date


points(x0[1]+360, x0[2], pch=21 ,bg=3)



image(lon+360,lat,(dl.grid[,,j])-(temp[j,3]-temp[j,2]),zlim=c(-5,5), add=F);  # this part should vary w/date
image(lon+360,lat,(mid.grid[,,j])-(temp[j,3]+temp[j,2])/2,zlim=c(-5,5),add=T);  # this part should vary w/date
contour(bathy$lon+360,bathy$lat, t(bathy$data), levels=c(0), add=T, col=1)




day0 = temp[1,1]

dsls =  day.last.sol(temp[16+90,1])

plot(sol.err(1.5,.001,0,dsls,1))

# what is the relationship between day length and latitude?
# x0 is the position. 
# day is in POSIXct form
# output is length of day, in fraction of day

par(mfrow=c(3,1))
plot(alldays,daylen(c(-70,42),alldays))
plot(seq(40,50,length=100),daylen(cbind(rep(-72,100),seq(40,50,length=100)),alldays[1]))
plot(sunriset(cbind(rep(-72,100),seq(40,50,length=100)), alldays[1], direction='sunrise', POSIXct.out=TRUE)$day)

# put it together
len = length(alldays)
day0 = temp[1,1]
dsls = sapply(1:167, function(i) day.last.sol(temp[i,1]))

plot(sol.err(1.5,.001,0,dsls,1))

serr = abs(sol.err(1.5,.001,0,dsls,1))

x0=c(-69.373,41.457)


par(mfrow=c(6,1))

# predicted sunsets according to solstice error at a single point using observations of sunrise
plot(daylen(cbind(rep(x0[1],100),rnorm(100,x0[2],sqrt(serr[1:100])),alldays[1:100]))*24*60+temp[1:100,2])
lines(daylen(cbind(rep(x0[1],100),rep(x0[2],100),alldays[1:100]))*24*60+temp[1:100,2])

plot(alldays[1:100],daylen(cbind(rep(x0[1],100),rnorm(100,x0[2],sqrt(serr[1:100])),alldays[1:100]))*24*60+temp[1,2])
lines(alldays[1:100],daylen(cbind(rep(x0[1],100),rep(x0[2],100)),alldays[1:100])*24*60+temp[1,2])

# predicted sunsets according to solstice error at a single point using observations of sunrise
plot(temp[1:100,3]-daylen(cbind(rep(x0[1],100),rnorm(100,x0[2],sqrt(serr[1:100])),alldays[1:100]))*24*60, ylim=c(500,2000))
lines(temp[1:100,3]-daylen(cbind(rep(x0[1],100),rep(x0[2]),alldays[1:100]))*24*60,col=2)
points(daylen(cbind(rep(x0[1],100),rnorm(100,x0[2],sqrt(serr[1:100])),alldays[1:100]))*24*60+temp[1:100,2])
lines(daylen(cbind(rep(x0[1],100),rep(x0[2]),alldays[1:100]))*24*60+temp[1:100,2],col=2)



par(mfrow=c(2,2))
# plot latitude vs. sunrise and latitude vs. sunrise w/error
plot(alldays[1:100],rep(x0[2],100), ylab = 'Latitude')
lines(alldays[1:100],rnorm(100,x0[2],sqrt(serr[1:100])))
title('Latitude error w/solstice model')

# plot just the error
plot(alldays[1:100],rep(x0[2],100)-rnorm(100,x0[2],sqrt(serr[1:100])), ylab = 'Diff')
title('Latitude error w/solstice model')

# plot the day length with latitude error
dlerr = daylen(cbind(rep(x0[1],100),rnorm(100,x0[2],sqrt(serr[1:100])),alldays[1:100]))*24*60
plot(alldays[1:100],dlerr, ylab='minutes')
title('day length w/error')

# Plot 100 iterations
plot(dlerr,ylim=c(400,800))
for(i in 1:100){
  lines(daylen(cbind(rep(x0[1],100),rnorm(100,x0[2],sqrt(serr[1:100])),alldays[1:100]))*24*60)
}

# plot the daylength without error
dlne = daylen(x0,alldays[1:100])*24*60

plot(alldays[1:100],daylen(x0,alldays[1:100])*24*60, ylab='minutes')
title('daylength no error')

# sunrise no error
sr = sunriset(matrix(as.numeric(x0),nrow=1), alldays[1:100], direction='sunset', POSIXct.out=F)*24*60

plot(alldays[1:100],sr-dlne, ylab='minutes')
title('sunrise no error')

# sunrise no error
sr = sunriset(matrix(as.numeric(x0),nrow=1), alldays[1:100], direction='sunrise', POSIXct.out=F)*24*60

plot(alldays[1:100],sr+dlne, ylab='minutes')
title('sunset no error')


# sunrise with error
srerr = sunriset(cbind(rep(x0[1],100),rnorm(100,x0[2],sqrt(serr[1:100])),alldays[1:100]),alldays[1:100], direction='sunset', POSIXct.out=F)*24*60

plot(alldays[1:100],srerr, ylab='minutes')
lines(alldays[1:100],sr,typ='o',col=2)
title('sunrise with error')

# sunset with error
sserr = sunriset(cbind(rep(x0[1],100),rnorm(100,x0[2],sqrt(serr[1:100])),alldays[1:100]), alldays[1:100], direction='sunrise', POSIXct.out=F)*24*60

plot(alldays[1:100],sserr+dlne, ylab='minutes')
lines(alldays[1:100],sunriset(matrix(as.numeric(x0),nrow=1), alldays[1:100], direction='sunset', POSIXct.out=F)*24*60,typ='o',col=2)
title('sunset with error')

# day length with sunrise and sunset error at a single point over time
estsr = sunriset(cbind(rep(x0[1],100),rnorm(100,x0[2],sqrt(serr[1:100])),alldays[1:100]), alldays[1:100], direction='sunset', POSIXct.out=F)*24*60-daylen(matrix(as.numeric(x0),nrow=1),alldays[1:100])*24*60

estss = sunriset(cbind(rep(x0[1],100),rnorm(100,x0[2],sqrt(serr[1:100])),alldays[1:100]), alldays[1:100], direction='sunrise', POSIXct.out=F)*24*60+daylen(matrix(as.numeric(x0),nrow=1),alldays[1:100])*24*60

plot(alldays[1:100],estsr,ylim=c(500,1500), typ='o', bg=2, pch=21)
lines(alldays[1:100],estss, typ='o', bg=4, pch=21)
lines(alldays[1:100],ss, typ='l')
lines(alldays[1:100],sr, typ='l')
title('sunrise and sunset with error')


#======================== From Uffe========================================#
# t(lat, lon) = sunriset(x0, day, direction='sunrise', POSIXct.out=TRUE)$day

sr1 = sunriset(matrix(as.numeric(x0),nrow=1), alldays[1], direction='sunrise', POSIXct.out=TRUE)$day*24*60
Tobs = temp[1,2]  # observed SR from tag

tau=-100:100 # realted to crepuscular period

# so:


cdawn = crepuscule(matrix(as.numeric(x0),nrow=1), alldays, solarDep=12, direction='dawn', POSIXct.out=TRUE)$day*24*60
sr1 = sunriset(matrix(as.numeric(x0),nrow=1), alldays, direction='sunrise', POSIXct.out=TRUE)$day*24*60
noon = solarnoon(matrix(as.numeric(x0),nrow=1), alldays, POSIXct.out=TRUE)$day*24*60
ss1 = sunriset(matrix(as.numeric(x0),nrow=1), alldays, direction='sunset', POSIXct.out=TRUE)$day*24*60
cdusk = crepuscule(matrix(as.numeric(x0),nrow=1), alldays, solarDep=12, direction='dusk', POSIXct.out=TRUE)$day*24*60

plot(alldays,cdawn, ylim=c(500,1500), typ='l', ylab='min after midnight GMT')
lines(alldays,sr1,col=2)
lines(alldays,noon,col=3)
lines(alldays,ss1,col=4)
lines(alldays,cdusk)
tau = mean(sr1-cdawn)
like.sr = exp( (Tobs-sr1)/tau ) #if Tobs<t(lat,lon), else l(lat,lon) = 0.


