calc.woa <- function(pdt, woa.dir, ptt, sdx){
  # compare tag data to ohc map and calculate likelihoods
  
  #' @param: pdt is variable containing tag-collected PDT data
  #' @param: woa.dir is local directory where woa ncdf is stored
  #' @return: likelihood is array of likelihood surfaces representing
  #' matches between daily tag-based pdt and woa maps
  
  # calculate midpoint of tag-based min/max temps
  pdt$MidTemp <- (pdt$MaxTemp + pdt$MinTemp) / 2
  
  udates <- unique(as.Date(pdt$Date))
  
  for(i in 1:length(udates)){
    # define time based on tag data
    time <- as.Date(udates[i])
    
    nc <- open.ncdf(paste(woa.dir, 'woa13_25deg_global.nc', sep=''))
    dat <- get.var.ncdf(nc, 'temperature')
    depth <- get.var.ncdf(nc, 'Depth')
    
    pdt.i <- pdt[which(as.Date(pdt$Date) == time),]
    
    # calculate daily isotherm based on tag data
    if(iso.def == FALSE) isotherm <- min(pdt.i$MidTemp, na.rm = T)
    
    dat[dat<isotherm] <- NA
    
    # Perform hycom integration
    dat <- dat - isotherm
    ohc <- cp * rho * apply(dat, 1:2, sum, na.rm = T) / 10000 
    
    # perform tag data integration
    tag <- approx(pdt.i$Depth, pdt.i$MidTemp, xout = depth)
    tag <- tag$y - isotherm
    tag.ohc <- cp * rho * sum(tag, na.rm = T) / 10000
    
    # compare hycom to that day's tag-based ohc
    #lik.dt <- matrix(dtnorm(ohc, tag.ohc, sdx, 0, 150), dim(ohc)[1], dim(ohc)[2])
    lik <- dnorm(ohc, tag.ohc, sdx) 
    print(paste(max(lik), time))
    
    # result should be array of likelihood surfaces
    if(i == 1){
      
      likelihood <- as.array(lik)
      
    } else{
      
      likelihood <- abind(likelihood, lik, along = 3)
      
    }
    
    print(paste(time, ' finished.', sep=''))
    
  }
  
  # return ohc likelihood surfaces as an array
  return(likelihood)
  
}
