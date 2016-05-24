
calc.pdt <- function(pdt, dat, lat, lon, g, depth, raster = 'stack', dateVec){
  
  ##  This program matches depth temperature profiles collected by a WC PSAT
  ##  tag to climatological profiles from the World Ocean Atlas.
  
  #' @param pdt is -PDT data from WC psat tag summarizing depth/temperature
  #'        data over a programmed time interval
  #' @param dat is monthly global 1/4deg climatology data from WOA13
  #' @param lat is vector of latitudes from dat
  #' @param lon is vector of longitudes from dat
  #' @param raster is character indicating whether likelihood 'array',
  #'        'stack' or 'brick' should be output
  #' @return lik is array of likelihoods for depth-temp profile
  #'        matching between tag data and WOA
  
  udates <- unique(pdt$Date)
  T <- length(udates)
   
  pdt$MidTemp <- (pdt$MaxTemp + pdt$MinTemp) / 2
  
  L.pdt <- array(0, dim = c(dim(dat)[1:2], length(dateVec)))
  
  for(i in 1:T){
    # define time based on tag data
    time <- udates[i]
    pdt.i <- pdt[which(pdt$Date == time),]
    y <- pdt.i$Depth[!is.na(pdt.i$Depth)] #extracts depth from tag data for day i
    y[y<0] <- 0
    
    if (length(y) >= 3){
      x <- pdt.i$MidTemp[!is.na(pdt.i$Depth)]  #extract temperature from tag data for day i
      pdtMonth <- as.numeric(format(as.Date(pdt.i$Date), format='%m'))[1]
      
      dat.i = dat[,,,pdtMonth] #extract months climatology
      
      # do the regression
      fit <- locfit(x ~ y)
      n = length(y)
      
      # find the standard depth levels that correspond to the depths the tag data is measured at
      # use the which.min
      depIdx = apply(as.data.frame(y), 1, FUN=function(z) which.min((z-depth)^2))
      woaDep <- depth[depIdx] 
      
      # make predictions based on the regression model earlier for the temperature at standard WOA depth levels for low and high temperature at that depth
      pred = predict(fit, newdata = woaDep, se = T, get.data = T)
      
      df <- data.frame(depth = woaDep, low=pred$fit-pred$se.fit*sqrt(n), 
                       high=pred$fit+pred$se.fit*sqrt(n))
      
      #depIdx <- findInterval(y, depth) #locates climatology dep points nearest to tag's recorded depths
      #datDep = depth[depIdx] 
      #tag = approx(y, x, xout=datDep, rule=2) #interpolates temperatures in y to relevant WOA depths
      #names(tag) = list('y', 'x')
      
      #sdx <- sd((pdt.i$MaxTemp - pdt.i$MinTemp), na.rm = T)
      sdx <- .7
      
      
      for (b in 1:length(depIdx)){
        # sequence
        dt = seq(df$low[b], df$high[b], length=10)
        
        # likelihood array for one depth
        lik0 = aaply(dt, 1, .fun = function(z) dnorm(x[b], dat[,,depIdx[b],pdtMonth], sdx))
        lik.b.int = apply(lik0, 2:3, sum)
        
        # calculate likelihood at each depth for a given tag time point
        #lik.b <- dnorm(dat[,, b, pdtMonth], tag$x[which(depIdx == b)], sdx) 
        #lik.b <- (lik.b / max(lik.b, na.rm = T)) - .05
        if(b == 1){
          lik.pdt <- as.array(lik.b.int)
        } else{
          lik.pdt <- abind(lik.pdt, lik.b.int, along = 3)
        }
      }
      
      lik.pdt <- apply(lik.pdt, 1:2, prod)

      idx <- which(dateVec == as.Date(time))
      L.pdt[,,idx] = lik.pdt
      
      print(time)
      
    } else{
      # what to do if y<3?
    }
    
  }

  
  crs <- "+proj=longlat +datum=WGS84 +ellps=WGS84"
  list.pdt <- list(x = lon, y = lat, z = L.pdt)
  ex <- extent(list.pdt)
  L.pdt <- brick(list.pdt$z, xmn=ex[1], xmx=ex[2], ymn=ex[3], ymx=ex[4], transpose=T, crs)
  L.pdt <- flip(L.pdt, direction = 'y')
 
  # make L.pdt match resolution/extent of g
  row <- dim(g$lon)[1]
  col <- dim(g$lon)[2]
  ex <- extent(c(min(g$lon[1,]), max(g$lon[1,]), min(g$lat[,1]), max(g$lat[,1])))
  crs <- "+proj=longlat +datum=WGS84 +ellps=WGS84"
  rasMatch <- raster(ex, nrows=row, ncols=col, crs = crs)
  L.pdt <- spatial_sync_raster(L.pdt, rasMatch)

  if(raster == 'brick'){
    s <- L.pdt
  } else if(raster == 'stack'){
    s <- stack(L.pdt)
  } else if(raster == 'array'){
    s <- raster::as.array(L.pdt, transpose = T)
  }
  
  print(class(s))
  return(s)
  
}