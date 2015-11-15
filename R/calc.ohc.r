calc.ohc = function(pdt,isotherm='',ohc.dir,ptt,sdx,plot=F){
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
  
  require(ncdf); require(abind)
  
  # constants for OHC calc
  cp = 3.993 # kJ/kg*C
  rho = 1025 # kg/m3
  
  # calculate midpoint of tag-based min/max temps
  pdt$MidTemp <- (pdt$MaxTemp + pdt$MinTemp) / 2
  
  udates <- unique(as.Date(pdt$Date))
  
  if(isotherm!='') iso.def = TRUE else iso.def = FALSE
  
  for(i in 1:length(udates)){
    # define time based on tag data
    time <- as.Date(udates[i])
    
    nc = open.ncdf(paste(ohc.dir,ptt,'_',time,'.nc',sep=''))
    dat = get.var.ncdf(nc, 'temperature')
    
    pdt.i <- pdt[which(as.Date(pdt$Date)==time),]
    
    # calculate daily isotherm based on tag data
    if(iso.def==FALSE) isotherm = min(pdt.i$MidTemp,na.rm=T)
    
    dat[dat<isotherm] = NA
    
    # Perform hycom integration
    dat = dat - isotherm
    ohc = cp*rho*apply(dat, 1:2, sum, na.rm=T)/10000 
    
    # perform tag data integration
    tag = pdt.i$MidTemp - isotherm
    tag.ohc = cp*rho*sum(tag,na.rm=T)/10000
    
    # compare hycom to that day's tag-based ohc
    lik = dnorm(ohc, tag.ohc, sdx) # how to represent sd of tag-based ohc?
    print(paste(max(lik),time))
    # result should be array of likelihood surfaces
    if(i==1){
      likelihood = as.array(lik)
      
      # also get lon, lat, depth, variables in case we need them later
      #lon.length = get.var.ncdf(nc, 'X')
      #lat.length = get.var.ncdf(nc, 'Y')
      #lon = seq(lon[1], lon[2], length = length(lon.length))
      #lat = seq(lat[1], lat[2], length = length(lat.length))
      #depth = get.var.ncdf(nc, 'Depth')
      
      if(plot==T){
        pdf('lydia likelihood2.pdf')
        image.plot(lik,zlim=c(0,1))
        title(paste(time))
      }
      
    } else{
      likelihood = abind(likelihood,lik,along=3)
      if(plot==T){
        image.plot(lik,zlim=c(0,1))
        title(paste(time))
      }
    }
    
    print(paste(time,' finished.',sep=''))
    
  }
    
  
  # maybe add plot option in which each day is added to a pdf on file
  if(plot==T) dev.off()
  # return ohc likelihood surfaces as an array
  return(likelihood)
  
}

