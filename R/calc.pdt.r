
calc.pdt <- function(pdt, dat, lat, lon){
  
  ##  This program matches depth temperature profiles collected by a WC PSAT
  ##  tag to climatological profiles from the World Ocean Atlas.
  
  #' @param pdt is -PDT data from WC psat tag summarizing depth/temperature
  #'        data over a programmed time interval
  #' @param dat is monthly global 1/4deg climatology data from WOA13
  #' @param lat is vector of latitudes from dat
  #' @param lon is vector of longitudes from dat
  #' @return lik is array of likelihoods for depth-temp profile
  #'        matching between tag data and WOA
  
  pdt$MidTemp <- (pdt$MaxTemp + pdt$MinTemp) / 2
  
  udates <- unique(pdt$Date)
  
  for(i in 1:length(udates)){
    # define time based on tag data
    time <- udates[i]
    pdt.i <- pdt[which(pdt$Date == time),]
    y <- pdt.i$Depth[!is.na(pdt.i$Depth)] #extracts depth from tag data for day i
    y[y<0] <- 0
    
    if (length(y) > 3){
      x <- pdt.i$MidTemp[!is.na(pdt.i$Depth)]  #extract temperature from tag data for day i
      pdtMonth <- as.numeric(format(pdt.i$Date, format='%m'))[1]
      
      dat.i = dat[,,,pdtMonth] #extract months climatology
      
      depIdx <- findInterval(y, depth) #locates climatology dep points nearest to tag's recorded depths
      datDep = depth[depIdx] 
      tag = approx(y, x, xout=datDep, rule=2) #interpolates temperatures in y to relevant WOA depths
      names(tag) = list('y', 'x')
      
      for (b in depIdx){
        # calculate likelihood at each depth for a given tag time point
        lik.b <- dnorm(dat[,, b, pdtMonth], tag$x[which(depIdx == b)], .5) 
        #lik.b <- (lik.b / max(lik.b, na.rm = T)) - .05
        if(min(which(depIdx == b)) == 1){
          lik.pdt <- as.array(lik.b)
        } else{
          lik.pdt <- abind(lik.pdt, lik.b, along = 3)
        }
      }
      
      lik.pdt <- apply(lik.pdt, 1:2, prod)
      lik.pdt <- (lik.pdt / max(lik.pdt, na.rm = T)) - .05

    } else{
      # what to do if y<3?
    }
    
    if(i == 1){
      lik <- lik.pdt
    } else{
      lik <- abind(lik, lik.pdt, along = 3)
    }
    
    print(time)
    
  }
  
  return(lik)
  
}
