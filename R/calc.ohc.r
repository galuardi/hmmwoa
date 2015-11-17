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
  
  for(i in vec){#1:length(udates)){
    # define time based on tag data
    time <- as.Date(udates[i])
    
    nc = open.ncdf(paste(ohc.dir,ptt,'_',time,'.nc',sep=''))
    dat = get.var.ncdf(nc, 'temperature')
    depth = get.var.ncdf(nc, 'Depth')
    # also get lon, lat, depth, variables in case we need them later
    lon.length = get.var.ncdf(nc, 'X')
    lat.length = get.var.ncdf(nc, 'Y')
    lon = seq(lon[1], lon[2], length = length(lon.length))
    lat = seq(lat[1], lat[2], length = length(lat.length))
    
    pdt.i <- pdt[which(as.Date(pdt$Date)==time),]
    
    # calculate daily isotherm based on tag data
    if(iso.def==FALSE) isotherm = min(pdt.i$MidTemp,na.rm=T)
    
    dat[dat<isotherm] = NA
    
    # Perform hycom integration
    dat = dat - isotherm
    ohc = cp*rho*apply(dat, 1:2, sum, na.rm=T)/10000 
    
    # perform tag data integration
    tag = approx(pdt.i$Depth,pdt.i$MidTemp,xout=depth)
    tag = tag$y - isotherm
    tag.ohc = cp*rho*sum(tag,na.rm=T)/10000
    
    # compare hycom to that day's tag-based ohc
    lik = dnorm(ohc, tag.ohc, sdx) 
    print(paste(max(lik),time))
    
    # result should be array of likelihood surfaces
    if(i==1){
      likelihood = as.array(lik)
      
      if(plot==T){
        pdf('lydia ohc mis.pdf')
        par(mfrow=c(1,2))
        image.plot(lon,lat,lik)
        if(length(which(sdays==time))>0){
          points(spot[which(sdays==time),c(8,7)],col='white')
        }
        title(paste(time))
        plot(dat.i[1:19],depth[1:19],ylim=c(1000,0),type='l',xlab='temp',ylab='depth',lwd=2,xlim=c(8,31))
        lines(pdt.i$MidTemp,pdt.i$Depth,col='blue',lwd=2)
        
      }
      
    } else{
      likelihood = abind(likelihood,lik,along=3)
      if(plot==T){
        par(mfrow=c(1,2))
        image.plot(lon,lat,lik)
        if(length(which(sdays==time))>0){
          points(spot[which(sdays==time),c(8,7)],col='white')
        }
        title(paste(time))
        plot(dat.i[1:19],depth[1:19],ylim=c(1000,0),type='l',xlab='temp',ylab='depth',lwd=2,xlim=c(8,31))
        lines(pdt.i$MidTemp,pdt.i$Depth,col='blue',lwd=2)
        
      }
    }
    
    print(paste(time,' finished.',sep=''))
    
  }
    
  
  if(plot==T) dev.off()
  
  # return ohc likelihood surfaces as an array
  return(likelihood)
  
}
