library(fields); library(lubridate); library(GeoLight); library(reshape2)

matrix(c(as.POSIXct('2008-10-15 06:12:00',tz='utc'),as.POSIXct('2008-10-15 18:39:00'),1),1,3)

try[1,1] = as.POSIXct('2008-10-15 06:12:00', tz='UTC')
try[1,2] = as.POSIXct('2008-10-15 18:39:00', tz='UTC')
try[2,1] = as.POSIXct('2008-10-15 18:39:00', tz='UTC')
try[2,2] = as.POSIXct('2008-10-16 06:12:00', tz='UTC')


# try WOods Hole
# 0631 rise, 1837 EST set
try[1,1] = as.POSIXct('2016-07-23 09:29:00', tz='UTC')
try[1,2] = as.POSIXct('2016-07-24 00:08:00', tz='UTC')
try[2,1] = as.POSIXct('2016-07-24 00:08:00', tz='UTC')
try[2,2] = as.POSIXct('2016-07-24 09:31:00', tz='UTC')
coord(try)
#[1,] -70.49592 30.29507
#[2,] -70.74360 30.00940

tloc <- matrix(c(-70.6731,41.5265),1,2) # woods hole
tdte <- as.POSIXct('2016-07-23 12:00:00', tz='UTC')
sunriset(tloc,tdte,direction='sunrise',POSIXct.out=T)
sunriset(tloc,tdte,direction='sunset',POSIXct.out=T)
# spot on timing calculations

## RESULT:
# GeoLight has spot on Longitude and Lat is off about 10 deg



#====
## NOW TRY AN EXAMPLE TAG DATE THAT ISN'T IN THE EQUINOX
#====
geo.try <- light.geo[1:2,] # example tag day
geo.try[1,2] <- geo.try[2,1]
geo.try <- geo.try[1,]
geo.try[1,1] <- geo.try[1,1] + (60*60*4) # add 4 hours to get to UTC, maybe?
geo.try[1,2] <- geo.try[1,2] + (60*60*4)
coord(geo.try)
# [1,] -60.26903 45.57531 
# and its close...ish

geo.try <- light.geo[151:152,] # example tag day
geo.try[1,2] <- geo.try[2,1]
geo.try <- geo.try[1,]
geo.try[1,1] <- geo.try[1,1] + (60*60*4) # add 4 hours to get to UTC, maybe?
geo.try[1,2] <- geo.try[1,2] + (60*60*4)
coord(geo.try)
# [1,] -61.92447 38.40528
tloc <- matrix(c(-67.9475,27.9648),1,2) # spot location from this day, quality 2
tdte <- as.POSIXct('2015-12-31 12:00:00', tz='UTC')
srt <- sunriset(tloc,tdte,direction='sunrise',POSIXct.out=T);srt[1,2]
# [1] "2015-12-31 11:22:46 UTC"
geo.try[,1]
# "2015-12-31 10:55:00 UTC"
sst <- sunriset(tloc,tdte,direction='sunset',POSIXct.out=T);sst[1,2]
# [1] "2015-12-31 21:46:45 UTC"
geo.try[,2]
# [1] "2015-12-31 21:26:15 UTC"


#====
## LETS TRY THE WHOLE TAG DATASET
#===

# need to build a dataframe result containing
# colnames() <- list('tagsr','tagss','geolat','geolon','spotlat','spotlon','calcsr','calcss')

# first we clean the tag data
u.yday <- unique(light$yday)
res <- as.data.frame(matrix(NA, length(u.yday), 3))
for(i in 1:length(u.yday)){
  light.i <- light[which(light$yday == u.yday[i]),]

  # are there multiple (or none) of either dawn or dusk?
  srise <- ifelse(length(which(light.i$Type == 'Dawn')) == 0, NA, min(which(light.i$Type=='Dawn')))
  sset <- ifelse(length(which(light.i$Type == 'Dusk')) == 0, NA, min(which(light.i$Type=='Dusk')))
  
  if(!is.na(srise)){
    srise <- light.i[srise,]
  } else{
    srise <- as.data.frame(matrix(NA, 1, ncol(light.i)))
  }
  
  if(!is.na(sset)){
    sset <- light.i[sset,]
  } else{
    sset <- as.data.frame(matrix(NA, 1, ncol(light.i)))
  }
  
  # put the correct dawn/dusk in geotag df format
  if(any(is.na(srise) | is.na(sset))){
    
  } else{
    res[i,] <- c(srise$dtime, sset$dtime, 1) 
  }
  
}

res[,1] <- as.POSIXct(res[,1], origin = '1970-01-01', tz = 'UTC')
res[,2] <- as.POSIXct(res[,2], origin = '1970-01-01', tz = 'UTC')
colnames(res) <- list('tFirst','tSecond','type')
tag.crds <- coord(res)
res$yday <- yday(res[,1])
res$geoLat <- tag.crds[,2]
res$geoLon <- tag.crds[,1]

## MAYBE WE NEED TO CONVERT THESE TO UTC, CANT TELL IF THEY'RE ACTUALLY IN THAT NOW OR NOT

# also need to trim down the spot data; it's already speed filtered, just need to go to unique dates
spot$dtime <- as.POSIXct(spot$dtime,tz='UTC')
spotTrim <- as.data.frame(matrix(NA, length(u.yday), 5))
for(t in 1:length(u.yday)){
  if(length(which(spot$yday %in% u.yday[t])) > 0){
    spot.t <- spot[min(which(spot$yday %in% u.yday[t])),]
    #print(spot.t)
    spotTrim[t,c(1:3)] <- spot.t[1,c(6,3,4)]
    #print(spotTrim[c(t-1,t),])
  } else{}

}
spotTrim[,1] <- as.POSIXct(spotTrim[,1], origin = '1970-01-01', tz = 'UTC')
spotTrim <- spotTrim[which(!is.na(spotTrim[,2])),]

spot.sr <- sunriset(as.matrix(spotTrim[,c(3,2)]),spotTrim[,1], direction = 'sunrise', POSIXct.out = T)
spot.sst <- sunriset(as.matrix(spotTrim[,c(3,2)]),spotTrim[,1], direction = 'sunset', POSIXct.out = T)
spotTrim[,4] <- spot.sr$time
spotTrim[,5] <- spot.sst$time
colnames(spotTrim) <- list('spotTime','spotLat','spotLon','spotSRT','spotSST')
spotTrim$yday <- yday(spotTrim$spotTime)

# now we want to merge the spotTrim and res, resulting data frames
all <- merge(res, spotTrim, by='yday', all.x = T)
all.sel <- all[which(all$geoLon < -30 & all$geoLon > -100),]
pdf('geolight.pdf',height=12,width=6)
par(mfrow=c(1,2))
plot(all.sel$geoLon,all.sel$geoLat,pch=16, col='blue')
world(add=T)
lines(all.sel$spotLon,all.sel$spotLat,col='red')
plot(all.sel$geoLon,type='l',col='blue')
lines(all.sel$spotLon,col='red')
dev.off()

# need to figure out what's happening here; def due to structure of light data
# yday              tFirst             tSecond type      geoLat    geoLon            spotTime spotLat
# 1    1 2016-01-01 15:33:45 2016-01-01 17:12:30    1  72.6292511 -64.92995 2016-01-01 08:43:51 28.1413
# 2    2 2016-01-02 19:35:00 2016-01-02 15:37:30    1 -57.3094931  96.91189 2016-01-02 11:28:51 28.3071

