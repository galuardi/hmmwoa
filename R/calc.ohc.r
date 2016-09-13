#' Calculate Ocean Heat Content (OHC) probability surface
#' 
#' Compare tag data to OHC map and calculate likelihoods
#' 
#' @param pdt is variable containing tag-collected PDT data
#' @param isotherm default '' in which isotherm is calculated on the fly based 
#'   on daily tag data. Otherwise, numeric isotherm constraint can be specified 
#'   (e.g. 20 deg C).
#' @param ohc.dir local directory where get.hycom downloads are stored.
#' @param dateVec is vector of dates from tag to pop-up in 1 day increments.
#'   
#' @return likelihood is raster brick of likelihood surfaces representing 
#'   estimated position based on tag-based OHC compared to calculated OHC using 
#'   HYCOM

calc.ohc <- function(pdt, isotherm = '', ohc.dir, dateVec, bathy = TRUE){

  start.t <- Sys.time()
  
  # constants for OHC calc
  cp <- 3.993 # kJ/kg*C <- heat capacity of seawater
  rho <- 1025 # kg/m3 <- assumed density of seawater
  
  # calculate midpoint of tag-based min/max temps
  pdt$MidTemp <- (pdt$MaxTemp + pdt$MinTemp) / 2
  
  # get unique time points
  udates <- unique(pdt$Date)
  T <- length(udates)
  
  if(isotherm != '') iso.def <- TRUE else iso.def <- FALSE
  
  for(i in 1:T){
    
    time <- udates[i]
    pdt.i <- pdt[which(pdt$Date == time),]
    print(pdt.i)
    
    # open day's hycom data
    nc <- ncdf::open.ncdf(paste(ohc.dir, ptt,'_', as.Date(time), '.nc', sep=''))
    dat <- ncdf::get.var.ncdf(nc, 'water_temp')
    
    if(i == 1){
      depth <- ncdf::get.var.ncdf(nc, 'depth')
      lon <- ncdf::get.var.ncdf(nc, 'lon')
      lat <- ncdf::get.var.ncdf(nc, 'lat')
    }
    
    #extracts depth from tag data for day i
    y <- pdt.i$Depth[!is.na(pdt.i$Depth)] 
    y[y<0] <- 0
    
    #extract temperature from tag data for day i
    x <- pdt.i$MidTemp[!is.na(pdt.i$Depth)]  
    
    # use the which.min
    depIdx = unique(apply(as.data.frame(pdt.i$Depth), 1, FUN=function(x) which.min((x - depth) ^ 2)))
    hycomDep <- depth[depIdx]
    
    if(bathy){
      mask <- dat[,,max(depIdx)]
      mask[is.na(mask)] <- NA
      mask[!is.na(mask)] <- 1
      for(bb in 1:length(depth)){
        dat[,,bb] <- dat[,,bb] * mask
      }
    }

    # make predictions based on the regression model earlier for the temperature at standard WOA depth levels for low and high temperature at that depth
    fit.low <- locfit::locfit(pdt.i$MinTemp ~ pdt.i$Depth)
    fit.high <- locfit::locfit(pdt.i$MaxTemp ~ pdt.i$Depth)
    n = length(hycomDep)
    
    pred.low = predict(fit.low, newdata = hycomDep, se = T, get.data = T)
    pred.high = predict(fit.high, newdata = hycomDep, se = T, get.data = T)
    
    # data frame for next step
    df = data.frame(low = pred.low$fit - pred.low$se.fit * sqrt(n),
                    high = pred.high$fit + pred.high$se.fit * sqrt(n),
                    depth = hycomDep)
    print(df)
    
    # isotherm is minimum temperature recorded for that time point
    if(iso.def == FALSE) isotherm <- min(df$low, na.rm = T)
    
    # perform tag data integration at limits of model fits
    minT.ohc <- cp * rho * sum(df$low - isotherm, na.rm = T) / 10000
    maxT.ohc <- cp * rho * sum(df$high - isotherm, na.rm = T) / 10000
    print(minT.ohc)
    print(maxT.ohc)

    # Perform hycom integration
    dat[dat<isotherm] <- NA
    dat <- dat - isotherm
    ohc <- cp * rho * apply(dat[,,depIdx], 1:2, sum, na.rm = T) / 10000 
    ohc[ohc == 0] <- NA
    
    # calc sd of OHC
    # focal calc on mean temp and write to sd var
    t = Sys.time()
    r = raster::flip(raster::raster(t(ohc)), 2)
    sdx = raster::focal(r, w = matrix(1, nrow = 9, ncol = 9),
                        fun = function(x) sd(x, na.rm = T))
    sdx = t(raster::as.matrix(raster::flip(sdx, 2)))
    print(paste('finishing sd for ', time,'. Section took ', Sys.time() - t))
    
    # compare hycom to that day's tag-based ohc
    t = Sys.time()
    lik.ohc <- likint3(ohc, sdx, minT.ohc, maxT.ohc)
    print(paste('finishing lik.ohc for ', time,'. Section took ', Sys.time() - t))
    
    if(i == 1){
      # result will be array of likelihood surfaces
      L.ohc <- array(0, dim = c(dim(lik.ohc), length(dateVec)))
    }
    
    idx <- which(dateVec == as.Date(time))
    L.ohc[,,idx] = lik.ohc

    print(paste(time, ' finished.', sep=''))
    
  }

  raster.t <- Sys.time()

  crs <- "+proj=longlat +datum=WGS84 +ellps=WGS84"
  list.ohc <- list(x = lon-360, y = lat, z = L.ohc)
  ex <- raster::extent(list.ohc)
  L.ohc <- raster::brick(list.ohc$z, xmn=ex[1], xmx=ex[2], ymn=ex[3], ymx=ex[4], transpose=T, crs)
  L.ohc <- raster::flip(L.ohc, direction = 'y')

  print(paste('OHC calculations and raster manipulation took ', Sys.time() - start.t, ' and ', Sys.time() - raster.t, ', respesctively.'))
  
  # return ohc likelihood surfaces
  return(L.ohc)
  
}

