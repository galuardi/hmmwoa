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
setwd('~/Documents/WHOI/Data/WhiteSharks/2013/121325/') 

# READ IN TAG DATA
ptt <- '121325'

# TAG/POPUP DATES AND LOCATIONS (dd, mm, YYYY, lat, lon)
#iniloc <- data.frame(matrix(c(13, 10, 2015, 41.575, -69.423, 
#                              24, 2, 2016, 26.6798, -69.0147), nrow = 2, ncol = 5, byrow = T))

iniloc <- data.frame(matrix(c(3, 3, 2013, 30.3917, -81.3802, 
                              31, 8, 2013, 30.668, -79.972), nrow = 2, ncol = 5, byrow = T))

colnames(iniloc) = list('day','month','year','lat','lon')
tag <- as.POSIXct(paste(iniloc[1,1], '/', iniloc[1,2], '/', iniloc[1,3], sep=''), format = '%d/%m/%Y')
pop <- as.POSIXct(paste(iniloc[2,1], '/', iniloc[2,2], '/', iniloc[2,3], sep=''), format = '%d/%m/%Y')
# VECTOR OF DATES FROM DATA. THIS WILL BE THE TIME STEPS, T, IN THE LIKELIHOODS
dateVec <- as.Date(seq(tag, pop, by = 'day')) 

# READ IN DATA FROM WC FILES
#myDir <- '~/Documents/WHOI/RCode/hmmwoa/inst/extdata/' # WHERE YOUR DATA LIVES, THIS IS THE EXAMPLE DATA
#light.data <- read.wc(ptt, wd = myDir, type = 'light', tag=tag, pop=pop); 
light.data <- read.wc(121325, tag=tag, pop=pop,type='light')
light.udates <- light.data$udates; light.data <- light.data$data

#----------------------------------------------------------------------------------#
# LIGHT LIKELIHOOD
# Light-based Longitude Likelihood
#----------------------------------------------------------------------------------#

# SET SPATIAL LIMITS, IF DESIRED
sp.lim <- list(lonmin = -95, lonmax = -45, latmin = 10, latmax = 55)
locs.grid <- setup.locs.grid(sp.lim, res='one')

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
sr.grid[,,1:365] = sapply(1:365, function(i) matrix(sunriset(xy, fyear[i], direction="sunrise", POSIXct.out=TRUE)$day,length(lon),length(lat)))
ss.grid[,,1:365] = sapply(1:365, function(i) matrix(sunriset(xy, fyear[i], direction="sunset", POSIXct.out=TRUE)$day,length(lon),length(lat)))

crs <- "+proj=longlat +datum=WGS84 +ellps=WGS84"
list.ras <- list(x = lon, y = lat, z = sr.grid*24*60)
ex <- raster::extent(list.ras)
sr.ras <- raster::brick(list.ras$z, xmn=ex[1], xmx=ex[2], ymn=ex[3], ymx=ex[4], transpose=T, crs)
sr.ras <- raster::flip(sr.ras, direction = 'y')

list.ras <- list(x = lon, y = lat, z = ss.grid*24*60)
ss.ras <- raster::brick(list.ras$z, xmn=ex[1], xmx=ex[2], ymn=ex[3], ymx=ex[4], transpose=T, crs)
ss.ras <- raster::flip(ss.ras, direction = 'y')

# need to be able to cut SRSS times from tag that aren't within limits of the grid
min.sr <- sapply(1:365, function(i) cellStats(sr.ras[[i]],stat='min',na.rm=T))
max.sr <- sapply(1:365, function(i) cellStats(sr.ras[[i]],stat='max',na.rm=T))
min.ss <- sapply(1:365, function(i) cellStats(ss.ras[[i]],stat='min',na.rm=T))
max.ss <- sapply(1:365, function(i) cellStats(ss.ras[[i]],stat='max',na.rm=T))

# make some calculations on the tag data: yday, dtime, etc
light <- light.data[,c('Day','Time','Type')]
light$dtime <- dmy_hms(paste(light$Day, light$Time, sep = ' '))
light$yday <- yday(light$dtime)
light$daymins <- minute(light$dtime) + (hour(light$dtime) * 60)
lightDates <- as.Date(format(light$dtime, '%Y-%m-%d'))

for (t in 1:length(dateVec)){
  light.t <- light[which(lightDates %in% dateVec[t]),]
  
  if(length(light.t[,1]) <= 1){
    if(length(light.t[,1]) == 1){
      if(light.t$Type == 'Dawn'){
        light[which(lightDates %in% dateVec[t] & light$Type == 'Dawn'), 6] <- NA
      }
      if(light.t$Type == 'Dusk'){
        light[which(lightDates %in% dateVec[t] & light$Type == 'Dusk'), 6] <- NA
      }
    }

  } else{
    didx <- light.t$yday[1]
    # what's the SR?
    sr <- light.t$daymins[which(light.t$Type == 'Dawn')]
    ss <- light.t$daymins[which(light.t$Type == 'Dusk')]
    
    if(length(sr) == 0){
      
    } else{
      if(length(sr) > 1){
        sr <- sr[1]
      }
      
      if(sr < min.sr[didx] | sr > max.sr[didx]){
        sr <- NA
      }
      
      light[which(lightDates %in% dateVec[t] & light$Type == 'Dawn'), 6] <- sr
      
    }
    

    if(length(ss) == 0){
      
    } else{
      
      if(length(ss) > 1){
        ss <- ss[1]
      }
      
      if(ss < min.ss[didx] | ss > max.ss[didx]){
        ss <- NA
      }
      
      light[which(lightDates %in% dateVec[t] & light$Type == 'Dusk'), 6] <- ss
      
    }

  }
  
}
 


####


# now add spot srss time calculations
spot <- read.table('~/Documents/WHOI/Data/WhiteSharks/2013/121325/121420-Locations.csv', sep=',', header = T)
spot$dtime <- as.POSIXct(spot$Date, format=findDateFormat(spot$Date), tz='UTC')
spot$yday <- yday(spot$dtime)
spot$daymins <- minute(spot$dtime)+hour(spot$dtime)*60
d1 <- as.POSIXct('1900-01-02') - as.POSIXct('1900-01-01')
didx <- spot$dtime >= (tag + d1) & spot$dtime <= (pop - d1)
spot <- spot[didx,]

spot.xy = SpatialPoints(spot[,c(8,7)], proj4string=CRS("+proj=longlat +datum=WGS84"))
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
spot.m <- melt(spot,id.vars=c('Ptt','Date','Latitude','Longitude','Quality','dtime','yday','daymins'),measure.vars=c('srt','sst','srt.ex','sst.ex'))

cols <- gg_color_hue(6)
ggplot(light, aes(x = dtime, y = daymins, colour = Type))+
  geom_point(size = 3) +
  geom_line(data=spot.m, aes(x=dtime, y=value, colour=variable))+
  scale_color_manual(values=c(cols[1:5],'#000000',cols[6],'#000000'))#+
#geom_line(data=pred.d, aes(x=dtime, y=fit, colour='black'))+
#geom_line(data=pred.d, aes(x=dtime, y=fit+2*se.fit, colour='black'))+
#geom_line(data=pred.d, aes(x=dtime, y=fit-2*se.fit, colour='black'))

#=======
## what if we tried a regression on the SRSS times, then added some SD bounds and got rid of all the other SRSS times outside that window?
# make predictions based on the regression model earlier for the temperature at standard WOA depth levels for low and high temperature at that depth
light$idx <- seq(1,length(light[,1]),by=1)
fit.dawn <- locfit::locfit(light$daymins[light$Type=='Dawn'] ~ light$idx[light$Type=='Dawn'])
fit.dusk <- locfit::locfit(light$daymins[light$Type=='Dusk'] ~ light$idx[light$Type=='Dusk'])
n = length(light$yday)

#mod4 = gls(youth ~ cadult*minwage,
#           correlation = corAR1(form=~1), data = un)
#summary(mod4)


pred.dawn = predict(fit.dawn, newdata = light$idx, se = T, get.data = T)
pred.dusk = predict(fit.dusk, newdata = light$idx, se = T, get.data = T)

# data frame for next step
dawn.df = data.frame(low = pred.dawn$fit - 2 * pred.dawn$se.fit,
                     fit = pred.dawn$fit,
                     high = pred.dawn$fit + 2 * pred.dawn$se.fit,
                     dtime = light$dtime)

dusk.df = data.frame(low = pred.dusk$fit - 2 * pred.dusk$se.fit,
                     fit = pred.dusk$fit,
                     high = pred.dusk$fit + 2 * pred.dusk$se.fit,
                     dtime = light$dtime)

dawn.df$dtime <- as.POSIXct(dawn.df$dtime, origin = '1970-01-01', tz = 'UTC')
dusk.df$dtime <- as.POSIXct(dusk.df$dtime, origin = '1970-01-01', tz = 'UTC')

cols <- gg_color_hue(6)
ggplot(light, aes(x = dtime, y = daymins, colour = Type))+
  geom_point(size = 3) +
  geom_line(data=spot.m, aes(x=dtime, y=value, colour=variable))+
  scale_color_manual(values=c(cols[1:5],'#000000',cols[6],'#000000'))+
  geom_line(data=dawn.df, aes(x=dtime, y=fit, colour='black'))+
  geom_line(data=dawn.df, aes(x=dtime, y=low, colour='black'))+
  geom_line(data=dawn.df, aes(x=dtime, y=high, colour='black'))+
  geom_line(data=dusk.df, aes(x=dtime, y=fit, colour='black'))+
  geom_line(data=dusk.df, aes(x=dtime, y=low, colour='black'))+
  geom_line(data=dusk.df, aes(x=dtime, y=high, colour='black'))#+
#geom_line(data=srss.mat, aes(x=dtime, y=max.sr, colour='blue'))










#==================
## END



# load known positions
spot <- read.table('~/Documents/WHOI/RData/Blues/2015/141256/spot_locs.csv', sep=',', header = T)

# iterate through daily likelihood calculation
for (t in 1:length(dateVec)){
  
  if(t==1){
    pdf('256lik.pdf',height=12,width=8)
  }
  
  light.t <- light[which(lightDates %in% dateVec[t]),]
  
  if(length(light.t[,1]) <= 1){
    print(paste(dateVec[t],' no light'))
  } else{
    didx <- light.t$yday[1]
    # what's the SR?
    sr <- light.t$daymins[which(light.t$Type == 'Dawn')]
    
    if(min.sr[didx] < sr < max.sr[didx]){
      # get the SD for this day, T
      f1 <- raster::focal(sr.ras[[didx]], w = matrix(1, nrow = 3, ncol = 3), fun = function(x) sd(x, na.rm = T))
      # the SR likelihood
      srlik <- liksrss(sr, srss = sr.ras[[didx]], srsd = f1)
    } else{
      srlik <- matrix(NA, ncol=dim(sr.ras[[didx]])[2], nrow=dim(sr.ras[[didx]])[1])
    }

    ss <- light.t$daymins[which(light.t$Type == 'Dusk')]
    if(min.ss[didx] < ss < max.ss[didx]){
      # and sunset
      f2 <- raster::focal(ss.ras[[didx]], w = matrix(1, nrow = 3, ncol = 3), fun = function(x) sd(x, na.rm = T))
      sslik <- liksrss(ss, srss = ss.ras[[didx]], srsd = f2)
      
    } else{
      sslik <- matrix(NA, ncol=dim(ss.ras[[didx]])[2], nrow=dim(ss.ras[[didx]])[1])
    }

    par(mfrow=c(3,1))
    plot(ss.ras[[didx]])
    world(add=T)
    contour(ss.ras[[didx]],nlevels=30,add=T)
    points(spot[which(spot$yday == didx)[1],c(4,3)],pch=16, col='red') # known position for this day
    title(paste(dateVec[t],' tag SS is ',ss, sep=''))
    plot(sr.ras[[didx]])
    world(add=T)
    contour(sr.ras[[didx]],nlevels=30,add=T)
    points(spot[which(spot$yday == didx)[1],c(4,3)],pch=16, col='red') # known position for this day
    title(paste(dateVec[t],' tag SR is ',sr, sep=''))
    plot(srlik * sslik) # the overlap is in Brazil
    world(add=T)
    points(spot[which(spot$yday == didx)[1],c(4,3)],pch=16, col='red') # known position for this day
    
  }
  print(t)
}
dev.off()


# and you'll see something is very wrong
plot(srlik * sslik) # the overlap is in Brazil
world(add=T)

plot(srlik)
world(add=T)
points(-64.863,34.617,pch=16, col='red') # known position for this day
plot(sslik)
# so what's happening here?
plot(sr.ras[[didx]])
world(add=T)
contour(sr.ras[[didx]],nlevels=30,add=T)
print(paste(sr,'is tag measured SR time. AND IS VERY WRONG. Seems off by one hour.'))

