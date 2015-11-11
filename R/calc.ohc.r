calc.ohc = function(tagdata,lon,lat,isotherm,ohc.dir){
  # compare tag data to ohc map and calculate likelihoods
  
  #' @param: tagdata is variable containing tag data
  #' @param: lon is vector of length 2 containing (min,max) lon
  #' @param: lat is vector of length 2 containing (min,max) lat
  #' @param: isotherm is default '' in which isotherm is calculated
  #' on the fly based on daily shark data. Otherwise, numeric isotherm
  #' constraint can be specified.
  #' @param: ohc.dir is local directory where get.hycom downloads are
  #' stored.
  #' @return: likelihood is array of likelihood surfaces representing
  #' matches between daily tag-based ohc and hycom ohc maps
  
  # define time based on tag data
  
  for(i in days_of_tag_data){
    nc = open.ncdf(paste())
    lon.length = get.var.ncdf(nc, 'X')
    lat.length = get.var.ncdf(nc, 'Y')
    lon = seq(lon[1], lon[2], length = length(lon.length))
    lat = seq(lat[1], lat[2], length = length(lat.length))
    depth = get.var.ncdf(nc, 'Depth')
    dat = get.var.ncdf(nc, 'temperature')
    cp = 3.993 # kJ/kg*C
    rho = 1025 # kg/m3
    
    if(isotherm=''){
      # calculate daily isotherm based on tag data
      isotherm = min(,na.rm=T)
    }
    
    dat[dat<isotherm] = NA
    
    # Perform hycom integration
    dat = dat-isotherm
    ohc = cp*rho*apply(dat, 1:2, sum, na.rm=T)/10000 
    
    # perform tag data integration
    tag = tag-isotherm
    tag.ohc = cp*rho*apply(tag,?,sum,na.rm=T)/10000
    
    # compare hycom to that day's tag-based ohc
    lik = dnorm(ohc, tag.ohc, sdx) # how to represent sd of tag-based ohc?
    
    # result should be array of likelihood surfaces
    if(i==1){
      likelihood = as.array(lik)
    } else{
      likelihood[,,i] = lik
    }
  }
  
  # return ohc likelihood surfaces as an array
  return(likelihood)
  
  
}

