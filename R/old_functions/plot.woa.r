
plot.woa <- function(lik, return.woa, filename, pdt = pdt, spot = '', write.dir = getwd()){
  
  # plot likelihood outputs from calc.ohc
  
  #' @param: likelihood is array of likelihood surfaces representing
  #' matches between daily tag-based ohc and hycom ohc maps
  #' @param filename is name of output pdf
  #' @param: write.dir is output directory to write pdf plots to. defaults
  #' to current working directory.
  #' @return: nothing in workspace. writes pdf of plots to disk at write.dir
  
  data(countriesLow)
  
  dat = return.woa$dat; lon = return.woa$lon; lat = return.woa$lat; depth = return.woa$depth
  
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
    pdtMonth <- as.numeric(format(time,'%m'))
    
    pdt.i <- pdt[which(pdt$Date == time),]
    
    if(spotExists) par(mfrow = c(1, 2))
    
    image.plot(lon, lat, lik[,,i] / max(lik[,,i], na.rm = T))
    plot(countriesLow, add=T, col = 'grey')
    
    title(paste(time))
    
    if(spotExists){
      # get day's spot locations and find mean position
      spot.i <- spot[which(sdays == time), c(8,7)]
      mlon <- mean(spot.i[,1])
      mlat <- mean(spot.i[,2])
      points(mlon, mlat, col = 'red')
      
      # panel 2
      # sample woa at that position
      idx <- c(which.min(abs(lon-(mlon))), which.min(abs(lat-(mlat))))
      dat.i <- dat[idx[1], idx[2],, pdtMonth]
      rm(list=c('mlon', 'mlat'))
      
      plot(dat.i, depth, ylim = c(1000,0), type = 'l', xlab = 'temp', ylab = 'depth', lwd = 2, xlim = c(4.5,31))
      lines(pdt.i$MidTemp, pdt.i$Depth, col='blue', lwd = 2)
      plot(pdt.i$MidTemp, pdt.i$Depth, col='blue', lwd = 2, ylim = c(1000, 0), xlab = 'temp',
           type = 'l', ylab = 'depth', xlim = c(4.5, 31))
      
    } #else{
    #plot(pdt.i$MidTemp, pdt.i$Depth, col='blue', lwd = 2, ylim = c(1000, 0), xlab = 'temp',
         #type = 'l', ylab = 'depth', xlim = c(4.5, 31))
      
    #}
    
   # legend(26,800,c('woa','tag'),col=c('black','blue'),lty=c(1,1),lwd=c(3,3))
    print(time)
    
  }
  
  dev.off()
  
}
