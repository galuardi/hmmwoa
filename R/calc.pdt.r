
calc.pdt <- function(pdt, dat, lat, lon){
  
  ##  This program matches depth temperature profiles collected by a WC PSAT
  ##  tag to climatological profiles from World Ocean Database. Program 
  ##  calculates a MSE for each oceanographic profile compared to each dive,
  ##  and writes the grid cell number, tag day, and MSE to variable
  
  #==============================
  # Creates MSE variable in the workspace where data will be stored
  #==============================
  
  #cellidx = gridPts$cellNo
  
  # mse: Mean Square Error value calculated for each oceanographic profile
  # fill array with MSE values for each day and each grid cell
  #mse = array(dim = c(lon = length(lon), lat = length(lat), uid = max(pdt[,9]))) #creates MSE array
  
  #============================
  #For loops calculate MSE for each profile for each tag uid
  #============================
  
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
      
      depIdx <- findInterval(y, depth)
      datDep = depth[depIdx] #locates climatology dep points nearest to tag's recorded depths
      tag = approx(y,x,xout=datDep,rule=2) #interpolates temperatures in y to 1 m intervals in DepInt
      names(tag) = list('y','x')
      
      for (b in depIdx){
        lik.b <- dnorm(dat[,, b, pdtMonth], tag$x[which(depIdx == b)], .5) 
        if(min(which(depIdx == b)) == 1){
          lik.pdt <- as.array(lik.b)
        } else{
          lik.pdt <- abind(lik.pdt, lik.b, along = 3)
        }
      }
      
      lik.pdt <- apply(lik.pdt, 1:2, prod)

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
