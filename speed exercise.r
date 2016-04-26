
## generate speeds for spot tagged fish
121418,121425,121420,132360,132352

spot = read.csv('~/Documents/WHOI/Data/WhiteSharks/2012/121418/121418-Locations.csv')
whites <- rbind(whites,spot)
dts <- as.POSIXct(spot$Date, format=findDateFormat(spot$Date))
dates <- c(dates,dts)
#spot <- blues[which(blues$ptt==141270),]
spot <- spot[which(spot$Type=='Argos'),]
dts <- as.POSIXct(spot$Date, format=findDateFormat(spot$Date))
#dts <- as.POSIXct(spot$date, format="%Y-%m-%d %H:%M:%S")

# remove duplicate date-times
spot$dts <- dts
spot <- spot[which(!duplicated(spot$dts)),] #modify this to take better loc class from duplicate times

# convert track to SpatialPointsDataFrame
coordinates(spot) <- ~Longitude + Latitude
#coordinates(spot) <- ~lon + lat
proj4string(spot)=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84")

# filter based on bathymetry - use raster::extract to sample
#   bathy at each point and filter >=0
ext <- raster::extract(bathy, spot) < 0
spot <- spot[ext,]

# then to trip object
tr <- trip(spot, c('dts','Ptt'))

# filter at maximum speed in km/hr
# e.g. 2 m/s = 7.2 km/hr
sf <- speedfilter2(tr, max.speed = speed)

#df1 <- list()
df1$spd132352 <- sf
str(df1)

all <- data.frame(NA, ncol=3, nrow=10^4)
for(i in 1:length(df1)){
  all[,2] <- df1[[1]]
}

# whites
121418,121425,121420,132360,132352
# blues
141261, 141264, 141268, 141270
# mako
141267

# blues/mako
blues <- read.table('/Users/Cam/Documents/WHOI/RData/sharkSiteData/AllArgosData.csv',sep=',',header=T)
# [1] 141261 141264 141267 141268 141270



df.spd <- as.data.frame(df1)
df.spd2 <- data.frame(t(df.spd))
df.spd <- melt(df.spd2)
df.spd <- df.spd[!is.na(df.spd$value),]
df.spd.filter <- df.spd[which(df.spd$value<=15),]
sm.density.compare(df.spd.filter$value,df.spd.filter$variable)

# overall speed density plots
# blues and mako
plot(density(df.spd.filter$value[which(df.spd.filter$variable == 'spd141261')]),col='red',
     main='compare blue speeds',xlab='Speed (km/hr)',xlim=c(0,15),ylim=c(0,.3))
lines(density(df.spd.filter$value[which(df.spd.filter$variable == 'spd141264')]),col='blue')
lines(density(df.spd.filter$value[which(df.spd.filter$variable == 'spd141267')]))
lines(density(df.spd.filter$value[which(df.spd.filter$variable == 'spd141270')]),col='green')
text(12,.25,'black is mako')

# whites
plot(density(df.spd.filter$value[which(df.spd.filter$variable == 'spd121418')]),col='red',
     main='compare whites speeds',xlab='Speed (km/hr)',xlim=c(0,15),ylim=c(0,.21))
lines(density(df.spd.filter$value[which(df.spd.filter$variable == 'spd121425')]),col='blue')
lines(density(df.spd.filter$value[which(df.spd.filter$variable == 'spd121420')]))
lines(density(df.spd.filter$value[which(df.spd.filter$variable == 'spd132360')]),col='green')
#lines(density(df.spd.filter$value[which(df.spd.filter$variable == 'spd132352')]),col='grey')
#text(12,.35,'grey is small coastal shark')


## now need to use BSAM to delineate behavioral modes then calc differential speeds based on behav mode
library(bsam)
blues.bsam <- blues[,c(1,2,5,4,3)]
blues.bsam <- blues[which(blues$ptt != 141267),]
colnames(blues.bsam) <- list('id','date','lc','lon','lat')
fit.corr.all = fitSSM(blues.bsam, model="hDCRWS", tstep=.5, adapt=3000, samples=1000, thin=10, chains=2)
fit.141261.dcrws <- fitSSM(blues.bsam[which(blues.bsam$id==141261),], model='DCRWS', tstep=.5, adapt=3000, samples=1000, thin=10, chains=2)
fit.141264.dcrws <- fitSSM(blues.bsam[which(blues.bsam$id==141264),], model='DCRWS', tstep=.5, adapt=3000, samples=1000, thin=10, chains=2)
fit.141268.dcrws <- fitSSM(blues.bsam[which(blues.bsam$id==141268),], model='DCRWS', tstep=.5, adapt=3000, samples=1000, thin=10, chains=2)
fit.141270.dcrws <- fitSSM(blues.bsam[which(blues.bsam$id==141270),], model='DCRWS', tstep=.5, adapt=3000, samples=1000, thin=10, chains=2)

whites1 <- whites[,c(2,4,6,8,7)]
whites1$Date <- factor(format(dates, format='%Y-%m-%d %H:%M:%S'))
colnames(whites1) <- colnames(blues.bsam)
whites.bsam <- whites1

dts <- as.POSIXct(whites.bsam$date, format=findDateFormat(whites.bsam$date))
whPtt <- c(121418,121425,121420,132360,132352)

for (i in whPtt){
  idx <- which(whites.bsam$id==i)
  spot1 <- spotFilter(whites.bsam[idx,], dts[idx], speed=14.4, bathy=bathy)
  
  if(i == whPtt[1]){
    spot <- spot1
  } else{
    spot <- rbind(spot,spot1)
  }
  
}

whites.bsam.filter <- spot[,c(1:5)]

fit.corr.all.wh = fitSSM(trim.whites, model="hDCRWS", tstep=.5, adapt=3000, samples=1000, thin=10, chains=2)
fit.121418.dcrws <- fitSSM(whites.bsam.filter[which(whites.bsam.filter$id==121418),], model='DCRWS', tstep=.5, adapt=3000, samples=1000, thin=10, chains=2)
fit.121425.dcrws <- fitSSM(whites.bsam.filter[which(whites.bsam.filter$id==121425),], model='DCRWS', tstep=.5, adapt=3000, samples=1000, thin=10, chains=2)
fit.121420.dcrws <- fitSSM(whites.bsam.filter[which(whites.bsam.filter$id==121420),], model='DCRWS', tstep=.5, adapt=3000, samples=1000, thin=10, chains=2)
fit.132360.dcrws <- fitSSM(whites.bsam.filter[which(whites.bsam.filter$id==132360),], model='DCRWS', tstep=.5, adapt=3000, samples=1000, thin=10, chains=2)
fit.132352.dcrws <- fitSSM(whites.bsam.filter[which(whites.bsam.filter$id==132352),], model='DCRWS', tstep=.5, adapt=3000, samples=1000, thin=10, chains=2)

trim.whites <- whites.bsam.filter[which(whites.bsam.filter$id != 121418 & 
                                          whites.bsam.filter$id != 132352),]

try1 <- whites.bsam.filter[which(whites.bsam.filter$id==132352),]
plot(try1$lon,try1$lat)
plot(countriesLow,add=T)


# essentially all one behav state when run this together
fit.corr.bl = fitSSM(blues.bsam, model="DCRW", tstep=.5, adapt=3000, samples=1000, thin=10, chains=2)

fit.141261 <- fitSSM(blues.bsam[which(blues.bsam$id==141261),], model='DCRW', tstep=.5, adapt=3000, samples=1000, thin=10, chains=2)

fit.sw.bl = fitSSM(blues.bsam, model="DCRWS", tstep=.5, adapt=3000, samples=1000, thin=10, chains=2)
fit.sw.bl = fitSSM(blues.bsam, model="DCRWS", tstep=.5, adapt=3000, samples=1000, thin=10, chains=2)




