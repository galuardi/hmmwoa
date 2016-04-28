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
#'   contains GPS, Argos, and GPE locations as applicable
#' @param iniloc is 2 x 5 dataframe containing day, month, year, lat, lon for 
#'   both tag and pop locations
#' @param g is output from setup.grid and indicates extent and resolution of 
#'   grid used to calculate likelihoods
#' @param raster is logical indicating whether to return likelihood as a raster
#'   or an array. Default is TRUE and will return a raster.
#' @param dateVec is vector of dates from tag to pop-up in 1 day increments.
#' @param errEll is logical indicating whether error ellipses should be 
#'   generated for light-based likelihoods as given from output of WC-GPE. False
#'   if only longitude should be used. If False, standard deviation on light 
#'   measurements is currently fixed at 0.7 deg longitude following Musyl et al 
#'   2001. Default is FALSE and will use longitude only.
#' @return L is an array of lon x lat likelihood surfaces (matrices) for each
#'   time point (3rd dimension)

calc.locs <- function(locs, iniloc, g, raster = TRUE, dateVec, errEll = F){

  T <- length(locs$Longitude)
  row <- dim(g$lon)[1]
  col <- dim(g$lon)[2]
  
  L.locs <- array(0, dim = c(col, row, length(dateVec)))
  
  # Initial location is known
  ilo <- which.min(abs(g$lon[1,] - iniloc$lon[1]))
  ila <- which.min(abs(g$lat[,1] - iniloc$lat[1]))
  L.locs[ilo, ila, 1] <- 1
  
  locDates <- as.Date(locs$Date, format = findDateFormat(locs$Date))
  
  # set up a larger grid to base ellipse on and to shift that error, if necessary (GPE only)
  ngrid <- rev(dim(g$lon))
  lon1 <- seq(min(g$lon[1,]) - 10, max(g$lon[1,]) + 10, by = g$dlo)
  lat1 <- seq(min(g$lat[,1]) - 10, max(g$lat[,1]) + 10, by = g$dla)
  g1 <- meshgrid(lon1, lat1)
  locs$Offset[which(is.na(locs$Offset))] <- 0
  
  for(t in 1:T){
    if(locs$Type[t] == 'GPS'){
      # if GPS exists then other forms of data for that time point are obsolete
      glo <- which.min(abs(g$lon[1,] - locs$Longitude[t]))
      gla <- which.min(abs(g$lat[,1] - locs$Latitude[t]))
      L.locs[glo, gla, which(dateVec == locDates[t])] <- 1
      
    } else if(locs$Type[t] == 'Argos'){
      # if Argos exists, GPE positions are obsolete
      alo <- which.min(abs(g$lon[1,] - locs$Longitude[t]))
      ala <- which.min(abs(g$lat[,1] - locs$Latitude[t]))
      L.locs[alo, ala, which(dateVec == locDates[t])] <- 1
      
    } else if(locs$Type[t] == 'GPE'){
      if(errEll == FALSE){
        # create longitude likelihood based on GPE data
        # SD for light-based longitude from Musyl et al. (2001)
        slon.sd <- 35 / 111 # Converting from kms to degrees
        
        # use normally distributed error from position using fixed std dev
        L.light <- dnorm(t(g$lon), locs$Longitude[t], slon.sd)
        
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
          L.ext <- raster::flip(raster::raster(g$lon, xmn = min(lon), 
                               xmx = max(lon), ymn = min(lat), ymx = max(lat)), direction = 'y')
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
          
          L.light <- dnorm(t(g$lon), locs$Longitude[t], slon.sd)
          
          L.locs[,,which(dateVec == locDates[t])] <- L.light
          
        }
         
      }
      
    } else{
      stop('No data type for this days location.')
      }
    #image.plot(lon,lat,L.locs[,,which(dateVec == locDates[t])], main=paste(dateVec[which(dateVec==locDates[t])]))
  }
  
  # End location is known
  elo <- which.min(abs(g$lon[1,] - iniloc$lon[2]))
  ela <- which.min(abs(g$lat[,1] - iniloc$lat[2]))
  L.locs[elo, ela, length(dateVec)] <- 1
  
  if(raster){
    # this performs some transformations to the likelihood array to convert to useable raster
    crs <- "+proj=longlat +datum=WGS84 +ellps=WGS84"
    list.locs <- list(x = g$lon[1,], y = g$lat[,1], z = L.locs)
    ex <- raster::extent(list.locs)
    L.locs <- raster::brick(list.locs$z, xmn = ex[1], xmx = ex[2], ymn = ex[3], ymx = ex[4], transpose = T, crs)
    L.locs <- raster::flip(L.locs, direction = 'y')
    L.locs <- raster::stack(L.locs)
  }
  
  print(class(L.locs))
  return(L.locs)
  
}

