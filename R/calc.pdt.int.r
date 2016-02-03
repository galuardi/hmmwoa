
calc.pdt.int <- function(pdt, dat = dat$dat, lat = dat$lat, lon = dat$lon, g, depth = dat$depth, raster = 'stack', dateVec){
  
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
# should be 2:(T-1)  
  
  # define time based on tag data
  time <- udates[i]
  pdt.i <- pdt[which(pdt$Date == time),]
  y <- pdt.i$Depth[!is.na(pdt.i$Depth)] #extracts depth from tag data for day i
  y[y<0] <- 0
  
  # if (length(y) >= 3){
    x <- pdt.i$MidTemp[!is.na(pdt.i$Depth)]  #extract temperature from tag data for day i

    # use the which.min
    # depIdx = apply(as.data.frame(deps), 1, FUN=function(x) which.min((x-stdDepth)^2))
    depIdx = apply(as.data.frame(pdt.i$Depth), 1, FUN=function(x) which.min((x-depth)^2))
    
    woaDep <- depth[depIdx] 
    
    # make predictions based on the regression model earlier for the temperature at standard WOA depth levels for low and high temperature at that depth
    
    fit.low <- locfit(pdt.i$MinTemp ~ pdt.i$Depth)
    fit.high <- locfit(pdt.i$MaxTemp ~ pdt.i$Depth)
    n = length(depth)
#     pred = predict(fit, newdata = depth, se = T, get.data = T)
    pred.low = predict(fit.low, newdata = depth, se = T, get.data = T)
    pred.high = predict(fit.high, newdata = depth, se = T, get.data = T)
    
    # data frame for next step
    df = data.frame(low=pred.low$fit[depIdx]-pred.low$se.fit[depIdx]*sqrt(n)
                    , high=pred.high$fit[depIdx]+pred.high$se.fit[depIdx]*sqrt(n)
                    , row.names = depth[depIdx])
    
#     df1 <- data.frame(low=pred$fit[depIdx]-pred$se.fit[depIdx]*sqrt(n), 
#                       high=pred$fit[depIdx]+pred$se.fit[depIdx]*sqrt(n),
#                       row.names = depth[depIdx])
    
    pdtMonth <- as.numeric(format(as.Date(pdt.i$Date), format='%m'))[1]
    
    dat.i = dat[,,,pdtMonth] #extract months climatology
    dat.i[is.na(dat.i)] = -9999
    
    # sdx <- .7
    sdr = (df[,2]-df[,1])/4

# setup the likelihood array for each day. Will have length (dim[3]) = n depths
    lik.pdt = array(1e-15, dim=c(dim(dat)[1], dim(dat)[2], length(depIdx)))
    
        
  for (b in 1:length(depIdx)) {
      # sequence
      # NO..
      # dt = seq(df$low[b], df$high[b], length = 10)
      
      # likelihood array for one depth
#       lik0 = aaply(
#         dt, 1, .fun = function(z)
#           dnorm(x[b], dat[,,depIdx[b],pdtMonth], sdx)
#       )
      # lik.b.int = apply(lik0, 2:3, sum)
      
    lik.pdt[,,b] = likint(dat.i[,,b], df[b,1], df[b,2], sdr[b])
      
      # calculate likelihood at each depth for a given tag time point
      #lik.b <- dnorm(dat[,, b, pdtMonth], tag$x[which(depIdx == b)], sdx)
      #lik.b <- (lik.b / max(lik.b, na.rm = T)) - .05
#       if (b == 1) {
#         lik.pdt <- as.array(lik.b.int)
#       } else{
#         lik.pdt <- abind(lik.pdt, lik.b.int, along = 3)
#       }
    }
      
      lik.pdt <- apply(lik.pdt, 1:2, prod)
      
      idx <- which(dateVec == as.Date(time))
      L.pdt[,,idx] = lik.pdt
      
      print(time)
      
#     } else{
#       # what to do if y<3?
#     }
    
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
