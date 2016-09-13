#' Calculate Depth-temperature based Likelihood
#' 
#' \code{calc.pdt.int} calculates likelihood of animal position based on 
#' summarized depth-temperature profiles
#' 
#' Tag-based depth-temperature profile summaries are compared to climatological 
#' profiles from the World Ocean Atlas (WOA) and "matched" to generate position 
#' likelihoods. This essentially attempts to estimate animal position based on 
#' the water mass it is in, particularly if extensive diving performs thorough 
#' sampling of the environment. However, remember the in situ data is being 
#' compared to climatological means.
#' 
#' @param pdt is -PDT data from WC psat tag summarizing depth/temperature data
#'   over a programmed time interval
#' @param dat is monthly global 1/4deg climatology data from WOA13
#' @param lat is vector of latitudes from dat
#' @param lon is vector of longitudes from dat
#' @param depth is vector of depths from dat
#' @param dateVec is vector of dates from tag to pop-up in 1 day increments.
#' 
#' @return array or raster of likelihoods for depth-temp profile
#'        matching between tag data and WOA
#'   

calc.pdt.int <- function(pdt, dat = dat, lat = lat, lon = lon, depth = depth, dateVec){
  
  start.t <- Sys.time()
  
  udates <- unique(pdt$Date)
  T <- length(udates)
  
  pdt$MidTemp <- (pdt$MaxTemp + pdt$MinTemp) / 2
  
  L.pdt <- array(0, dim = c(dim(dat)[1:2], length(dateVec)))
  
  for(i in 1:T){
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
    fit.low <- locfit::locfit(pdt.i$MinTemp ~ pdt.i$Depth)
    fit.high <- locfit::locfit(pdt.i$MaxTemp ~ pdt.i$Depth)
    n = length(depth[depIdx])
    
    pred.low = predict(fit.low, newdata = depth[depIdx], se = T, get.data = T)
    pred.high = predict(fit.high, newdata = depth[depIdx], se = T, get.data = T)
    
    # data frame for next step
    df = data.frame(low = pred.low$fit - pred.low$se.fit * sqrt(n),
                    high = pred.high$fit + pred.high$se.fit * sqrt(n),
                    depth = depth[depIdx])
    
    if(i == 1) pdtMonth <- as.numeric(format(as.Date(pdt.i$Date), format = '%m'))[1]
    
    newMonth <- as.numeric(format(as.Date(pdt.i$Date), format = '%m'))[1]
    
    if(i == 1 | newMonth != pdtMonth) {
      # calculates sd but "if" statement ensures it is only calculated at
      # the beginning and when the date switches into a new month
      # because it's relatively computationally intensive
      
      pdtMonth <- as.numeric(format(as.Date(pdt.i$Date), format = '%m'))[1]
      dat.i = dat[,,,pdtMonth] #extract months climatology

      # calculate sd using Le Bris neighbor method and focal()
      sd.i = array(NA, dim = dim(dat.i))
      
      for(ii in 1:57){
        r = raster::flip(raster::raster(t(dat.i[,,ii])), 2)
        f1 = raster::focal(r, w = matrix(1, nrow = 3, ncol = 3), fun = function(x) sd(x, na.rm = T))
        f1 = t(raster::as.matrix(raster::flip(f1, 2)))
        sd.i[,,ii] = f1
      } 
    }
    
    # setup the likelihood array for each day. Will have length (dim[3]) = n depths
    lik.pdt = array(NA, dim = c(dim(dat)[1], dim(dat)[2], length(depIdx)))
    
    for (b in 1:length(depIdx)) {
      #calculate the likelihood for each depth level, b
      lik.pdt[,,b] = likint3(dat.i[,,depIdx[b]], sd.i[,,depIdx[b]], df[b,1], df[b,2])
      
    }
    
    # multiply likelihood across depth levels for each day
    lik.pdt <- apply(lik.pdt, 1:2, prod, na.rm = F)
    
    # identify date index and add completed likelihood to L.pdt array    
    idx <- which(dateVec == as.Date(time))
    L.pdt[,,idx] = lik.pdt
    
  }
  
  crs <- "+proj=longlat +datum=WGS84 +ellps=WGS84"
  list.pdt <- list(x = lon, y = lat, z = L.pdt)
  ex <- raster::extent(list.pdt)
  L.pdt <- raster::brick(list.pdt$z, xmn = ex[1], xmx = ex[2], ymn = ex[3], ymx = ex[4], transpose = T, crs)
  L.pdt <- raster::flip(L.pdt, direction = 'y')
  
  print(paste('PDT calculations took ', Sys.time() - start.t, '.'))
  
  return(L.pdt)
  
}


