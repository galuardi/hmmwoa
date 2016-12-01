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

calc.locs <- function(locs, gps = NULL, iniloc, locs.grid, dateVec, errEll = TRUE, gpeOnly = TRUE){
  
  print(paste('Starting likelihood calculation...'))
  
  # get rid of non-GPE locations if gpeOnly == TRUE and check length of resulting locs file
  if(gpeOnly == TRUE){
    locs <- locs[which(locs$Type == 'GPE'),]
    if(length(locs[,1]) < 1){
      stop('No GPE positions available in input locs file. Make sure you have calculated GPE positions and exported the resulting -Locations.csv file from WC DAP/GPE2 software.')
    }
  } else{
    stop('Error: Currently gpeOnly must be set to TRUE otherwise duplicate dates are handled improperly.')
  }
  
  # set some date vectors for our input data
  locDates <- as.Date(locs$Date, format = findDateFormat(locs$Date))
  if(!is.null(gps)){
    gpsDates <- as.Date(gps$InitTime, format = findDateFormat(gps$InitTime))  
  } else{
    gpsDates = NULL
  }
  
  if(any(duplicated(locDates))){
    # run a simplify function
    locList <- simplifyLocs(locs, locDates)
    locs <- locList$locs; locDates <- locList$locDates
  }
  
  if(any(duplicated(gpsDates))){
    # run a simplify function
    locList <- simplifyLocs(gps, gpsDates)
    gps <- locList$locs; gpsDates <- locList$locDates
    
  }
  
  # set up results array
  row <- dim(locs.grid$lon)[1]
  col <- dim(locs.grid$lon)[2]
  L.locs <- array(0, dim = c(col, row, length(dateVec)))
  
  # add tag/pop locations as known
  ilo <- which.min(abs(locs.grid$lon[1,] - iniloc$lon[1]))
  ila <- which.min(abs(locs.grid$lat[,1] - iniloc$lat[1]))
  L.locs[ilo, ila, 1] <- 1   # Initial location is known
  
  elo <- which.min(abs(locs.grid$lon[1,] - iniloc$lon[2]))
  ela <- which.min(abs(locs.grid$lat[,1] - iniloc$lat[2]))
  L.locs[elo, ela, length(dateVec)] <- 1  # End location is known

  for(t in 2:(length(dateVec)) - 1){
    
    print(paste('Starting iterations through deployment period...'))
    
    if(!is.null(gps) & dateVec[t] %in% gpsDates){
      
      # if GPS exists then other forms of data for that time point are obsolete
      idx <- which(gpsDates == dateVec[t]) # set index to identify position in gps file
      glo <- which.min(abs(locs.grid$lon[1,] - gps$InitLon[idx]))
      gla <- which.min(abs(locs.grid$lat[,1] - gps$InitLat[idx]))
      L.locs[glo, gla, t] <- 1
      
    } else if(!is.null(locs) & dateVec[t] %in% locDates){
      # set index to identify position in locs file
      idx <- which(locDates == dateVec[t])
      
      if(locs$Type[idx] == 'GPS'){ #locs includes GPS
        # if GPS exists then other forms of data for that time point are obsolete
        glo <- which.min(abs(locs.grid$lon[1,] - locs$Longitude[idx]))
        gla <- which.min(abs(locs.grid$lat[,1] - locs$Latitude[idx]))
        L.locs[glo, gla, t] <- 1
        
      } else if(locs$Type[idx] == 'Argos'){ #locs includes Argos
        # if Argos exists, GPE positions are obsolete
        alo <- which.min(abs(locs.grid$lon[1,] - locs$Longitude[idx]))
        ala <- which.min(abs(locs.grid$lat[,1] - locs$Latitude[idx]))
        L.locs[alo, ala, t] <- 1
        
      } else if(locs$Type[idx] == 'GPE'){ #locs includes GPE
        
        if(errEll == FALSE){
          # create longitude likelihood based on GPE data
          slon.sd <- locs$Error.Semi.minor.axis[idx] / 1000 / 111 #semi minor axis
          # use normally distributed error from position using fixed std dev
          L.light <- dnorm(t(locs.grid$lon), locs$Longitude[idx], slon.sd)
          
          L.locs[,,t] <- L.light
          
        } else if(errEll == TRUE){
          
          L.locs[,,t] <- calc.errEll(locs[idx,], locs.grid)
          
        }
        
      } else{
        
        stop('Error: Error ellipse unspecified.')
        
      }
      
    } else{
      
      # no data so we skip this day
      
    }
    
  }
  
  print(paste('Making final likelihood raster...'))
  
  # this performs some transformations to the likelihood array to convert to useable raster
  crs <- "+proj=longlat +datum=WGS84 +ellps=WGS84"
  list.locs <- list(x = locs.grid$lon[1,], y = locs.grid$lat[,1], z = L.locs)
  ex <- raster::extent(list.locs)
  L.locs <- raster::brick(list.locs$z, xmn = ex[1], xmx = ex[2], ymn = ex[3], ymx = ex[4], transpose = T, crs)
  L.locs <- raster::flip(L.locs, direction = 'y')
  
  return(list(L.locs = L.locs, gpsIdx = which(dateVec %in% gpsDates)))
  
}
