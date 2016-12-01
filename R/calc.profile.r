#' Calculate Depth-temperature profile based likelihood
#' 
#' \code{calc.profile} calculates likelihood of animal position based on 
#' summarized depth-temperature profiles
#' 
#' Tag-based depth-temperature profile summaries are compared to climatological 
#' profiles from the World Ocean Atlas (WOA) or HYbrid Coordinate Ocean Model (HYCOM) 
#' and "matched" to generate position likelihoods. This essentially attempts to estimate animal position based on 
#' the water mass it is in, particularly if extensive diving performs thorough 
#' sampling of the environment. However, remember the in situ data is being 
#' compared to climatological means or the results of an oceanographic model.
#' 
#' @param pdt is -PDT data from WC psat tag summarizing depth/temperature data
#'   over a programmed time interval
#' @param dat is monthly global 1/4deg climatology data from WOA13. NULL if envType = 'hycom'.
#' @param lat is vector of latitudes from dat. NULL if envType = 'hycom'.
#' @param lon is vector of longitudes from dat. NULL if envType = 'hycom'.
#' @param dateVec is vector of dates from tag to pop-up in 1 day increments.
#' @param envType is character indicating whether to compare tag-based profile to World Ocean Atlas ('woa') or HYCOM ('hycom').
#' @param hycom.dir is path to stored HYCOM model outputs if envType = 'hycom'
#' 
#' @return raster of likelihoods for depth-temp profile
#'        matching between tag data and oceanographic profiles
#'   

calc.profile <- function(pdt, dat = NULL, lat = NULL, lon = NULL, dateVec, envType = 'woa', hycom.dir = NULL){
  
  options(warn=-1)
  start.t <- Sys.time()
  
   if(envType == 'woa'){
    if(is.null(dat) | is.null(lat) | is.null(lon)){
      stop('Error: dat, lat, lon must all be specified if envType == woa')
    }
    depth <- c(0, seq(2.5, 97.5, by=5), seq(112.5, 487.5, by=25), seq(525, 1475, by=50))
    
  } else if(envType == 'hycom'){
    if(is.null(hycom.dir)){
      stop('Error: hycom.dir must be specified if envType == hycom')
    }
    depth <- c(seq(0, 12, by=2), seq(15, 50, by=5), seq(60, 100, by=10), 125, 150, 200,
               250, 300, 350, seq(400, 1000, by=100), 1250, 1500, 2000, 2500, 3000)
    
  }
  
  udates <- unique(pdt$Date)
  T <- length(udates)
  
  pdt$MidTemp <- (pdt$MaxTemp + pdt$MinTemp) / 2
  
  print(paste('Starting iterations through time ', '...'))
  
  for(i in 1:T){
    
    # define time based on tag data
    time <- udates[i]
    pdt.i <- pdt[which(pdt$Date == time),]
    
    #extracts depth from tag data for day i
    y <- pdt.i$Depth[!is.na(pdt.i$Depth)] 
    y[y < 0] <- 0
    
    #extract temperature from tag data for day i
    x <- pdt.i$MidTemp[!is.na(pdt.i$Depth)]  
    
    # use the which.min
    depIdx = apply(as.data.frame(pdt.i$Depth), 1, FUN=function(x) which.min((x-depth)^2))
    
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
    
    if (envType == 'woa'){
      
      if(i == 1){
        pdtMonth <- as.numeric(format(as.Date(pdt.i$Date), format = '%m'))[1]
        
        L.prof <- array(0, dim = c(dim(dat)[1:2], length(dateVec)))
        
      } 
      
      newMonth <- as.numeric(format(as.Date(pdt.i$Date), format = '%m'))[1]
      
      if(i == 1 | newMonth != pdtMonth) {
        # calculates sd but "if" statement ensures it is only calculated at
        # the beginning and when the date switches into a new month
        # because it's relatively computationally intensive
        
        pdtMonth <- as.numeric(format(as.Date(pdt.i$Date), format = '%m'))[1]
        dat.i = dat[,,,pdtMonth] #extract months climatology
        
        # calculate sd using Le Bris neighbor method and focal()
        sd.i = array(NA, dim = c(dim(dat.i)[1:2], length(depth)))
        
        for(ii in 1:length(depth)){
          r = raster::flip(raster::raster(t(dat.i[,,ii])), 2)
          f1 = raster::focal(r, w = matrix(1, nrow = 3, ncol = 3), fun = function(x) sd(x, na.rm = T))
          f1 = t(raster::as.matrix(raster::flip(f1, 2)))
          sd.i[,,ii] = f1
        }
        
      }
      
    } else if (envType == 'hycom'){
      
      # load hycom data for this day
      # open day's hycom data
      t <- Sys.time()
      nc <- RNetCDF::open.nc(paste(hycom.dir, ptt,'_', as.Date(time), '.nc', sep=''))
      dat.i <- RNetCDF::var.get.nc(nc, 'water_temp') * RNetCDF::att.get.nc(nc, 'water_temp', attribute='scale_factor') + 
        RNetCDF::att.get.nc(nc, variable='water_temp', attribute='add_offset')
      
      if(i == 1){
        L.prof <- array(0, dim = c(dim(dat.i)[1:2], length(dateVec)))
        
        #depth <- RNetCDF::var.get.nc(nc, 'depth')
        lon <- RNetCDF::var.get.nc(nc, 'lon') - 360
        lat <- RNetCDF::var.get.nc(nc, 'lat')
      }
      
      # create a bathymetry mask
      mask <- dat.i[,,max(depIdx)]
      mask[is.na(mask)] <- NA
      mask[!is.na(mask)] <- 1
      for(bb in 1:length(depth)){
        dat.i[,,bb] <- dat.i[,,bb] * mask
      }
      
      # calculate sd using Le Bris neighbor method and focal()
      sd.i = array(NA, dim = c(dim(dat.i)[1:2], length(depIdx)))
      
      for(ii in 1:length(depIdx)){
        r = raster::flip(raster::raster(t(dat.i[,,depIdx[ii]])), 2)
        f1 = raster::focal(r, w = matrix(1, nrow = 9, ncol = 9), fun = function(x) sd(x, na.rm = T))
        f1 = t(raster::as.matrix(raster::flip(f1, 2)))
        sd.i[,,ii] = f1
      } 
      
    }
    
    print(paste('Calculating likelihood for ', as.Date(time), '...', sep=''))

    # setup the likelihood array for each day. Will have length (dim[3]) = n depths
    lik.pdt = array(NA, dim = c(dim(dat.i)[1], dim(dat.i)[2], length(depIdx)))
    
    if (envType == 'woa'){
      for (b in 1:length(depIdx)) {
        #calculate the likelihood for each depth level, b
        lik.pdt[,,b] = likint3(dat.i[,,depIdx[b]], sd.i[,,depIdx[b]], df[b,1], df[b,2])
        
      }
    } else{
      for (b in 1:length(depIdx)) {
        #calculate the likelihood for each depth level, b
        lik.pdt[,,b] = likint3(dat.i[,,depIdx[b]], sd.i[,,b], df[b,1], df[b,2])
        
      }
    }

    # multiply likelihood across depth levels for each day
    lik.pdt <- apply(lik.pdt, 1:2, prod, na.rm = F)
    
    # identify date index and add completed likelihood to L.pdt array    
    idx <- which(dateVec == as.Date(time))
    L.prof[,,idx] = (lik.pdt / max(lik.pdt, na.rm=T)) - 0.2
    
  }
  
  print(paste('Making final likelihood raster...'))
  
  crs <- "+proj=longlat +datum=WGS84 +ellps=WGS84"
  list.pdt <- list(x = lon, y = lat, z = L.prof)
  ex <- raster::extent(list.pdt)
  L.prof <- raster::brick(list.pdt$z, xmn = ex[1], xmx = ex[2], ymn = ex[3], ymx = ex[4], transpose = T, crs)
  L.prof <- raster::flip(L.prof, direction = 'y')
  
  options(warn = 2)
  return(L.prof)
  
}


