
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

