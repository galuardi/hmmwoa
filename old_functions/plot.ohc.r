
plot.ohc <- function(lik, filename, ohc.dir = ohc.dir, pdt = pdt, spot = '', write.dir = getwd()){
 
   # plot likelihood outputs from calc.ohc
  
  #' @param: likelihood is array of likelihood surfaces representing
  #' matches between daily tag-based ohc and hycom ohc maps
  #' @param filename is name of output pdf
  #' @param: write.dir is output directory to write pdf plots to. defaults
  #' to current working directory.
  #' @return: nothing in workspace. writes pdf of plots to disk at write.dir
  
  data(countriesLow)
  
  if(class(spot) == 'data.frame'){
    
    spotExists <- TRUE
    sdates <- as.POSIXct(spot$Date, format = findDateFormat(spot$Date))
    sdays <- as.Date(sdates)
    
  } else{spotExists <- FALSE}
  
  udates <- unique(pdt$Date)
  
  # calculate midpoint of tag-based min/max temps
  pdt$MidTemp <- (pdt$MaxTemp + pdt$MinTemp) / 2
  
  pdf(paste(write.dir, '/', filename, sep = ''), width = 11, height = 8)
  
  for(i in 1:length(udates)){

    time <- udates[i]
    
    nc <- open.ncdf(paste(ohc.dir, ptt, '_', as.Date(time), '.nc', sep=''))
    dat <- get.var.ncdf(nc, 'temperature')
    
    pdt.i <- pdt[which(pdt$Date == time),]
    
    if(i == 1){
      # get dat variables
      depth <- get.var.ncdf(nc, 'Depth')
      lon <- get.var.ncdf(nc, 'Longitude') - 360
      lat <- get.var.ncdf(nc, 'Latitude')

      }
    
    if(spotExists) par(mfrow = c(1, 2))
    
    image.plot(lon, lat, lik[,,i] / max(lik[,,i], na.rm = T))
    plot(countriesLow, add=T, col = 'grey')
    title(paste(time))
    
    if(spotExists){
      # get day's spot locations and find mean position
      spot.i <- spot[which(sdays == as.Date(time)), c(8,7)]
      mlon <- mean(spot.i[,1])
      mlat <- mean(spot.i[,2])
      points(mlon, mlat, col = 'white')
      
      # panel 2
      # sample hycom at that position
      idx <- c(which.min(abs(lon.idx-(mlon))), which.min(abs(lat.idx-(mlat))))
      dat.i <- dat[idx[1], idx[2],]
      rm(list=c('mlon', 'mlat'))
      
      plot(dat.i[1:19], depth[1:19], ylim = c(1000,0), type = 'l', xlab = 'temp', ylab = 'depth', lwd = 2, xlim = c(8,31))
      lines(pdt.i$MidTemp, pdt.i$Depth, col='blue', lwd = 2)
      legend(26,800,c('hycom','tag'),col=c('black','blue'),lty=c(1,1),lwd=c(3,3))
      
    }
    
    print(time)
    
  }
  
  dev.off()
    
}
