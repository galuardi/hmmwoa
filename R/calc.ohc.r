calc.ohc <- function(pdt, isotherm = '', ohc.dir, ptt, sdx){
  # compare tag data to ohc map and calculate likelihoods
  
  #' @param: tagdata is variable containing tag-collected PDT data
  #' @param: time is vector of unique dates (daily) used to step
  #' through the integration / ohc calculations
  #' @param: isotherm is default '' in which isotherm is calculated
  #' on the fly based on daily shark data. Otherwise, numeric isotherm
  #' constraint can be specified.
  #' @param: ohc.dir is local directory where get.hycom downloads are
  #' stored.
  #' @return: likelihood is array of likelihood surfaces representing
  #' matches between daily tag-based ohc and hycom ohc maps
  
  # constants for OHC calc
  cp <- 3.993 # kJ/kg*C
  rho <- 1025 # kg/m3
  
  # calculate midpoint of tag-based min/max temps
  pdt$MidTemp <- (pdt$MaxTemp + pdt$MinTemp) / 2
  
  udates <- unique(pdt$Date)
  
  if(isotherm != '') iso.def <- TRUE else iso.def <- FALSE
  
  for(i in 1:length(udates)){
    # define time based on tag data
    time <- udates[i]
    
    nc <- open.ncdf(paste(ohc.dir, ptt, '_', as.Date(time), '.nc', sep=''))
    dat <- get.var.ncdf(nc, 'temperature')
    depth <- get.var.ncdf(nc, 'Depth')
    
    pdt.i <- pdt[which(pdt$Date == time),]
    
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
    lik <- (lik / max(lik, na.rm = T)) - .05
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
