
calc.pdt.int <- function(pdt, dat = dat$dat, lat = dat$lat, lon = dat$lon, g, depth = dat$depth, raster = 'stack', dateVec){
  
  ##  This program matches depth temperature profiles collected by a WC PSAT
  ##  tag to climatological profiles from the World Ocean Atlas.
  
  #' @param pdt is -PDT data from WC psat tag summarizing depth/temperature
  #'        data over a programmed time interval
  #' @param dat is monthly global 1/4deg climatology data from WOA13
  #' @param lat is vector of latitudes from dat
  #' @param lon is vector of longitudes from dat
  #' @param g
  #' @param depth
  #' @param raster is character indicating whether likelihood 'array',
  #'        'stack' or 'brick' should be output
  #' @param dateVec
  #' 
  #' @return lik is array of likelihoods for depth-temp profile
  #'        matching between tag data and WOA
  
  require(locfit)
  
  udates <- unique(pdt$Date)
  T <- length(udates)
  
  pdt$MidTemp <- (pdt$MaxTemp + pdt$MinTemp) / 2
  
  L.pdt <- array(0, dim = c(dim(dat)[1:2], length(dateVec)))
  
for(i in 1:T){
  # should be 2:(T-1) 
  # CDB: no because we've already trimmed to tag+1, pop-1
  
  # define time based on tag data
  time <- udates[i]
  pdt.i <- pdt[which(pdt$Date == time),]
  
  #extracts depth from tag data for day i
  y <- pdt.i$Depth[!is.na(pdt.i$Depth)] 
  y[y<0] <- 0
  
  #extract temperature from tag data for day i
  x <- pdt.i$MidTemp[!is.na(pdt.i$Depth)]  

  # use the which.min
  depIdx = apply(as.data.frame(pdt.i$Depth), 1, FUN=function(x) which.min((x-depth)^2))
  woaDep <- depth[depIdx] 
    
  # make predictions based on the regression model earlier for the temperature at standard WOA depth levels for low and high temperature at that depth
  fit.low <- locfit(pdt.i$MinTemp ~ pdt.i$Depth)
  fit.high <- locfit(pdt.i$MaxTemp ~ pdt.i$Depth)
  n = length(depth[depIdx])

  pred.low = predict(fit.low, newdata = depth[depIdx], se = T, get.data = T)
  pred.high = predict(fit.high, newdata = depth[depIdx], se = T, get.data = T)
    
  # data frame for next step
  df = data.frame(low=pred.low$fit-pred.low$se.fit*sqrt(n)
                  , high=pred.high$fit+pred.high$se.fit*sqrt(n)
                  , depth = depth[depIdx])
  
  if(i == 1) pdtMonth <- as.numeric(format(as.Date(pdt.i$Date), format='%m'))[1]
  
  newMonth <- as.numeric(format(as.Date(pdt.i$Date), format='%m'))[1]
  
  if(i==1 | newMonth != pdtMonth) {
    pdtMonth <- as.numeric(format(as.Date(pdt.i$Date), format='%m'))[1]
    dat.i = dat[,,,pdtMonth] #extract months climatology
    #dat.i[is.na(dat.i)] = -9999

    # calculate sd using Le Bris neighbor method and focal()
    sd.i = array(NA,dim=dim(dat.i))
    for(ii in 1:57){
      r = flip(raster(t(dat.i[,,ii])),2)
      #plot(r, col = tim.colors(100))
      f1 = focal(r, w=matrix(1/9,nrow=3,ncol=3), fun=sd)
      f1 = t(as.matrix(flip(f1,2)))
      #plot(f1, add=T)
      sd.i[,,ii] = f1
    } 
  }
  
  # setup the likelihood array for each day. Will have length (dim[3]) = n depths
  lik.pdt = array(NA, dim=c(dim(dat)[1], dim(dat)[2], length(depIdx)))
        
  for (b in 1:length(depIdx)) {
    
    #calculate the likelihood for each depth level, b
    lik.pdt[,,b] = likint(dat.i[,,b], df[b,1], df[b,2], sd.i[,,b])
    
    image.plot(likint2(dat.i[,,b], sd.i[,,b], df[b,1], df[b,2]))
    
    likint2 <- function(woa, woasd, minT, maxT){
      wlist = array(1e-6, dim=c(dim(woa)[1], dim(woa)[2], 2))
      wlist[,,1] = woa
      wlist[,,2] = woasd
      wlist[is.na(wlist)] = 1e-6
      as.matrix(aaply(wlist, 1:2, .fun = function(x) integrate(dnorm, lower = minT, upper = maxT , mean = x[1], sd = x[2])$value))
    }
    
    #image.plot(dnorm(df[b,1],mean=dat.i[,,b],sd=.7))
    
    }
  
  # multiply likelihood across depth levels for each day
  #lik.pdt.naomit <- apply(lik.pdt, 1:2, function(x) prod(na.omit(x)))
  lik.pdt <- apply(lik.pdt, 1:2, prod, na.rm=T)
  
  # identify date index and add completed likelihood to L.pdt array    
  idx <- which(dateVec == as.Date(time))
  L.pdt[,,idx] = lik.pdt
  
  print(time)
  
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
