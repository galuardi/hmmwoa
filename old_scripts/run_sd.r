
# download directory
ohc.dir <- paste('~/Documents/WHOI/RData/HYCOM/sd/',sep = '')

# set our vars
lon = c(-90, -40)
lat = c(10, 55)
udates <- seq(as.Date('2016-01-01'),as.Date('2016-01-31'),by=1)

# download day's hycom data
for(i in 1:3){#length(udates)){
  time <- as.Date(udates[i])
  fname <- paste('glbu_08_91.1_temp_',time,'.nc',sep='')
  repeat{
    get.hycom(lon,lat,time,filename=fname,
              download.file=TRUE,dir=ohc.dir, vars = 'water_temp',
              type = 'a') # filenames based on dates from above
    #err <- try(open.ncdf(paste(ohc.dir,ptt,'_',time,'.nc',sep='')),silent=T)
    tryCatch({
      err <- try(open.ncdf(fname),silent=T)
    }, error=function(e){print(paste('ERROR: Download of data at ',time,' failed. Trying call to server again.',sep=''))})
    if(class(err) != 'try-error') break
  }
  
  # open day's hycom data
  setwd(ohc.dir)
  nc <- open.ncdf(fname)
  dat <- get.var.ncdf(nc, 'water_temp')
  
  if(i==1){
    jan <- array(NA, dim=c(dim(dat)[1:2],length(udates)))
    depth <- get.var.ncdf(nc, 'depth')
    lon.vec <- get.var.ncdf(nc, 'lon')
    lat.vec <- get.var.ncdf(nc, 'lat')
  }
  
  # isolate surface and store in array
  jan[,,i] <- dat[,,1]
   
}


japp <- apply(jan,1:2,mean)
image.plot(japp)
jvar <- apply(jan,1:2,function(x) var.hycom(x,mu=japp))

jan.err <- apply(jan, 1:2, function(x) (sd(x)^2)*3)
jan.sd <- apply(jan, 1:2, sd)

jet.colors <- colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan", "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))

brks <- c(seq(0,1.2,by=.05),max(jan.sd,na.rm=T))
image.plot(lon.vec-360,lat.vec,jan.sd,xlab='Lon',ylab='Lat',breaks=brks,col=jet.colors(length(brks)-1))
world(add=T,col='grey',fill=T)


http://ncss.hycom.org/thredds/ncss/GLBu0.08/expt_91.1?var=water_temp&north=50&west=-90&east=-40&south=10&disableProjSubset=on&horizStride=1&time_start=2016-01-13T00%3A00%3A00Z&time_end=2016-01-13T00%3A00%3A00Z&timeStride=1&vertCoord=&addLatLon=true&accept=netcdf
http://ncss.hycom.org/thredds/ncss/GLBu0.08/expt_91.1?var=water_temp&north=55.000000&west=-90.000000&east=-40.000000&south=10.000000&horizStride=1&time_start=2016-01-13T00%3A00%3A00Z&time_end=2016-01-13T00%3A00%3A00Z&timeStride=1&addLatLon=true&disableProjSubset=on&vertCoord=&accept=netcdf

