#' Calculate Position-based Likelihood
#' 
#' \code{calc.locs} calculates likelihood estimates for each day of animal tag 
#' data.
#' 
#' Light errors are parameterized using elliptical error values output in 
#' '-Locations.csv' (WC tags). GPS and Argos positions are also given a 
#' "likelihood" using this function but are currently both considered to be 
#' "known" positions without error.
#' 
#' @param locs is -Locations file output from DAP/Tag Portal for WC tags and 
#'   contains GPS, Argos, and GPE locations as applicable.
#' @param gps is -FastGPS file output from WC Tag Portal
#' @param iniloc is 2 x 5 dataframe containing day, month, year, lat, lon for 
#'   both tag and pop locations
#' @param dateVec is vector of dates from tag to pop-up in 1 day increments.
#' @param locs.grid is list output from \code{setup.locs.grid}
#' @param errEll is logical indicating whether error ellipses should be 
#'   generated for light-based likelihoods as given from output of WC-GPE. False
#'   if only longitude should be used. If False, standard deviation on light 
#'   measurements is currently fixed at 0.7 deg longitude following Musyl et al 
#'   2001. Default is FALSE and will use longitude only.
#' @param gpeOnly is logical. If TRUE, locs input is trimmed to GPE positions only. This is most applicable in scenarios with FastGPS data and you're adding this as a gps input (see "gps" argument).
#' @return L is an array of lon x lat likelihood surfaces (matrices) for each
#'   time point (3rd dimension)

calc.locs <- function(locs, gps = '', iniloc, locs.grid, dateVec, errEll = T, gpeOnly = T){

  if(gpeOnly == TRUE){
    locs <- locs[which(locs$Type == 'GPE'),]
    
    if(length(locs[,1]) < 1){
      stop('No GPE positions available in input locs file.')
    }
 
  }
  
  row <- dim(locs.grid$lon)[1]
  col <- dim(locs.grid$lon)[2]
  
  L.locs <- array(0, dim = c(col, row, length(dateVec)))
  
  # Initial location is known
  ilo <- which.min(abs(locs.grid$lon[1,] - iniloc$lon[1]))
  ila <- which.min(abs(locs.grid$lat[,1] - iniloc$lat[1]))
  L.locs[ilo, ila, 1] <- 1
  
  locDates <- as.Date(locs$Date, format = findDateFormat(locs$Date))
  if(gps != ''){
    gpsDates <- as.Date(gps$InitTime, format = findDateFormat(gps$InitTime))  
  }
  
  # set up a larger grid to base ellipse on and to shift that error, if necessary (GPE only)
  ngrid <- rev(dim(locs.grid$lon))
  lon1 <- seq(min(locs.grid$lon[1,]) - 10, max(locs.grid$lon[1,]) + 10, by = locs.grid$dlo)
  lat1 <- seq(min(locs.grid$lat[,1]) - 10, max(locs.grid$lat[,1]) + 10, by = locs.grid$dla)
  g1 <- meshgrid(lon1, lat1)
  locs$Offset[which(is.na(locs$Offset))] <- 0
  
  for(t in 2:(length(dateVec)-1)){
    if(gps != ''){
      
      if(gpeOnly == FALSE){
        warning('Specifying a gps input and having gpeOnly == FALSE may cause wacky things to happen. Recommended to switch gpeOnly to TRUE.')
      }
      
      if(any(dateVec[t] %in% gpsDates)){
        # if GPS exists then other forms of data for that time point are obsolete
        idx <- which(gpsDates == dateVec[t])
        glo <- which.min(abs(locs.grid$lon[1,] - gps$InitLon[idx]))
        gla <- which.min(abs(locs.grid$lat[,1] - gps$InitLat[idx]))
        L.locs[glo, gla, t] <- 1
        
      } else if(any(dateVec[t] %in% locDates)){
        # do the light thing
        if(errEll == FALSE){
          # create longitude likelihood based on GPE data
          # SD for light-based longitude from Musyl et al. (2001)
          slon.sd <- 35 / 111 # Converting from kms to degrees
          
          # use normally distributed error from position using fixed std dev
          L.light <- dnorm(t(locs.grid$lon), locs$Longitude[t], slon.sd)
          
          L.locs[,,which(dateVec == locDates[t])] <- L.light
          
        } else if(errEll == TRUE){
          # arithmetic converts from meters to degrees
          # add transformation due to projection?
          slon.sd <- locs$Error.Semi.minor.axis[t] / 1000 / 111 #semi minor axis
          L.light.lon <- dnorm(t(g1$X), locs$Longitude[t], slon.sd) # Longitude data
          slat.sd <- locs$Error.Semi.major.axis[t] / 1000 / 111 #semi major axis
          L.light.lat <- dnorm(t(g1$Y), locs$Latitude[t], slat.sd)
          
          #image.plot(g$lon[1,],g$lat[,1],L.light.lat*L.light.lon)
          
          L <- raster::flip(raster::raster(t(L.light.lat * L.light.lon), xmn = min(lon1), 
                                           xmx = max(lon1), ymn = min(lat1), ymx = max(lat1)), direction = 'y')
          
          # offset, assuming shift should be to the south
          shiftDist <- (-1 * (locs$Offset[t] / 1000 / 111))
          
          if(shiftDist >= -10){
            Ls <- raster::shift(L, y = shiftDist)
            L.ext <- raster::flip(raster::raster(locs.grid$lon, xmn = min(locs.grid$lon[1,]), 
                                                 xmx = max(locs.grid$lon[1,]), ymn = min(locs.grid$lat[,1]),
                                                 ymx = max(locs.grid$lat[,1])), direction = 'y')
            # create blank raster
            L.ext[L.ext <= 0] = 1
            
            # then crop our shifted raster
            Lsx <- raster::crop(Ls, L.ext)
            rr <- raster::resample(Lsx, L.ext)
            #image.plot(lon,lat,t(as.matrix(flip(rr,direction='y'))))
            L.locs[,,which(dateVec == locDates[t])] <- t(raster::as.matrix(raster::flip(rr, direction = 'y')))
            
          } else{
            # if supposed shift in error ellipse is >10 degrees, we revert to longitude only
            slon.sd <- 35/111 # Converting from kms to degrees
            
            L.light <- dnorm(t(locs.grid$lon), locs$Longitude[t], slon.sd)
            
            L.locs[,,which(dateVec == locDates[t])] <- L.light
            
          }
          
        }
        
      } else{} # do nothing
      
    } else{
      # this happens if no gps input is specified, the first two if statements here become obsolete if gpeOnly = TRUE
      
      if(locs$Type[t] == 'GPS'){
        # if GPS exists then other forms of data for that time point are obsolete
        glo <- which.min(abs(locs.grid$lon[1,] - locs$Longitude[t]))
        gla <- which.min(abs(locs.grid$lat[,1] - locs$Latitude[t]))
        L.locs[glo, gla, which(dateVec == locDates[t])] <- 1
        
      } else if(locs$Type[t] == 'Argos'){
        # if Argos exists, GPE positions are obsolete
        alo <- which.min(abs(locs.grid$lon[1,] - locs$Longitude[t]))
        ala <- which.min(abs(locs.grid$lat[,1] - locs$Latitude[t]))
        L.locs[alo, ala, which(dateVec == locDates[t])] <- 1
        
      } else if(locs$Type[t] == 'GPE'){
        if(errEll == FALSE){
          # create longitude likelihood based on GPE data
          # SD for light-based longitude from Musyl et al. (2001)
          slon.sd <- 35 / 111 # Converting from kms to degrees
          
          # use normally distributed error from position using fixed std dev
          L.light <- dnorm(t(locs.grid$lon), locs$Longitude[t], slon.sd)
          
          L.locs[,,which(dateVec == locDates[t])] <- L.light
          
        } else if(errEll == TRUE){
          # arithmetic converts from meters to degrees
          # add transformation due to projection?
          slon.sd <- locs$Error.Semi.minor.axis[t] / 1000 / 111 #semi minor axis
          L.light.lon <- dnorm(t(g1$X), locs$Longitude[t], slon.sd) # Longitude data
          slat.sd <- locs$Error.Semi.major.axis[t] / 1000 / 111 #semi major axis
          L.light.lat <- dnorm(t(g1$Y), locs$Latitude[t], slat.sd)
          
          #image.plot(g$lon[1,],g$lat[,1],L.light.lat*L.light.lon)
          
          L <- raster::flip(raster::raster(t(L.light.lat * L.light.lon), xmn = min(lon1), 
                                           xmx = max(lon1), ymn = min(lat1), ymx = max(lat1)), direction = 'y')
          
          # offset, assuming shift should be to the south
          shiftDist <- (-1 * (locs$Offset[t] / 1000 / 111))
          
          if(shiftDist >= -10){
            Ls <- raster::shift(L, y = shiftDist)
            L.ext <- raster::flip(raster::raster(locs.grid$lon, xmn = min(locs.grid$lon[1,]), 
                                                 xmx = max(locs.grid$lon[1,]), ymn = min(locs.grid$lat[,1]),
                                                 ymx = max(locs.grid$lat[,1])), direction = 'y')
            # create blank raster
            L.ext[L.ext <= 0] = 1
            
            # then crop our shifted raster
            Lsx <- raster::crop(Ls, L.ext)
            rr <- raster::resample(Lsx, L.ext)
            #image.plot(lon,lat,t(as.matrix(flip(rr,direction='y'))))
            L.locs[,,which(dateVec == locDates[t])] <- t(raster::as.matrix(raster::flip(rr, direction = 'y')))
            
          } else{
            # if supposed shift in error ellipse is >10 degrees, we revert to longitude only
            slon.sd <- 35/111 # Converting from kms to degrees
            
            L.light <- dnorm(t(locs.grid$lon), locs$Longitude[t], slon.sd)
            
            L.locs[,,which(dateVec == locDates[t])] <- L.light
            
          }
          
        }
        
      } else{
        stop('No data type for this days location.')
      }
    }
  }
  
  # End location is known
  elo <- which.min(abs(locs.grid$lon[1,] - iniloc$lon[2]))
  ela <- which.min(abs(locs.grid$lat[,1] - iniloc$lat[2]))
  L.locs[elo, ela, length(dateVec)] <- 1
  
  # this performs some transformations to the likelihood array to convert to useable raster
  crs <- "+proj=longlat +datum=WGS84 +ellps=WGS84"
  list.locs <- list(x = locs.grid$lon[1,], y = locs.grid$lat[,1], z = L.locs)
  ex <- raster::extent(list.locs)
  L.locs <- raster::brick(list.locs$z, xmn = ex[1], xmx = ex[2], ymn = ex[3], ymx = ex[4], transpose = T, crs)
  L.locs <- raster::flip(L.locs, direction = 'y')
  
  return(L.locs)
  
}

