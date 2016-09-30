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

calc.locs2 <- function(light = NULL, locs = NULL, gps = NULL, iniloc, locs.grid, dateVec, res = 1){
  
  start.t <- Sys.time()
  
  #=====
  # BUILD A SPOT FOR THE RESULTS TO GO
  #=====
  # set up results array
  row <- dim(locs.grid$lon)[1]
  col <- dim(locs.grid$lon)[2]
  
  # need lat/lon vectors. come from locs.grid??
  lon <- locs.grid$lon[1,]
  lat <- locs.grid$lat[,1]
  crs <- "+proj=longlat +datum=WGS84 +ellps=WGS84"
  L.grid = numeric(length = c(length(lon)*length(lat)*length(dateVec)))
  dim(L.grid) = c(length(lon),length(lat), length(dateVec))
  list.ras <- list(x = lon, y = lat, z = L.grid)
  ex <- raster::extent(list.ras)
  L.locs <- raster::brick(list.ras$z, xmn=ex[1], xmx=ex[2], ymn=ex[3], ymx=ex[4], transpose=T, crs)
  L.locs <- raster::flip(L.locs, direction = 'y')
  
  #L.locs <- array(0, dim = c(col, row, length(dateVec)))
  
  if(!is.null(light)){
    #==================
    # build the SRSS grids
    #==================
    
    # expand.grid and SpatialPoints establishes a grid
    xy = as.matrix(expand.grid(lon,lat))
    xy = SpatialPoints(xy, proj4string=CRS("+proj=longlat +datum=WGS84"))
    
    # now do the building and rasterize
    sr.grid = numeric(length = c(length(lon)*length(lat)*365))
    dim(sr.grid) = c(length(lon),length(lat), 365)
    ss.grid = sr.grid
    
    fyear = seq(ISOdate(year(dateVec[1]), 1, 1, tz = 'UTC'), ISOdate(year(dateVec[1]), 12, 31, tz = 'UTC'), 'day')
    sr.grid[,,1:365] = sapply(1:365, function(i) matrix(sunriset(xy, fyear[i], direction = "sunrise", POSIXct.out = TRUE)$day,length(lon),length(lat)))
    ss.grid[,,1:365] = sapply(1:365, function(i) matrix(sunriset(xy, fyear[i], direction = "sunset", POSIXct.out = TRUE)$day,length(lon),length(lat)))
    
    list.ras <- list(x = lon, y = lat, z = sr.grid*24*60)
    ex <- raster::extent(list.ras)
    sr.ras <- raster::brick(list.ras$z, xmn=ex[1], xmx=ex[2], ymn=ex[3], ymx=ex[4], transpose=T, crs)
    sr.ras <- raster::flip(sr.ras, direction = 'y')
    
    list.ras <- list(x = lon, y = lat, z = ss.grid*24*60)
    ss.ras <- raster::brick(list.ras$z, xmn=ex[1], xmx=ex[2], ymn=ex[3], ymx=ex[4], transpose=T, crs)
    ss.ras <- raster::flip(ss.ras, direction = 'y')

    # need to be able to cut SRSS times from tag that aren't within limits of the grid
    min.sr <- sapply(1:365, function(i) cellStats(sr.ras[[i]],stat='min',na.rm=T))
    max.sr <- sapply(1:365, function(i) cellStats(sr.ras[[i]],stat='max',na.rm=T))
    min.ss <- sapply(1:365, function(i) cellStats(ss.ras[[i]],stat='min',na.rm=T))
    max.ss <- sapply(1:365, function(i) cellStats(ss.ras[[i]],stat='max',na.rm=T))
    
    
    # make some calculations on the tag data: yday, dtime, etc
    light <- light[,c('Day','Time','Type')]
    light$dtime <- dmy_hms(paste(light$Day, light$Time, sep = ' '))
    light$yday <- yday(light$dtime)
    light$daymins <- minute(light$dtime) + (hour(light$dtime) * 60)
    light <- light[which(light$Type != ''),]
    lightDates <- as.Date(format(light$dtime, '%Y-%m-%d'))
    
    for(t in 2:(length(dateVec)) - 1){
      # data for this time step T
      light.t <- light[which(lightDates %in% dateVec[t]),]
      
      if(length(light.t[,1]) == 0){
        
      } else{
        if(length(light.t[,1]) == 1 & any(light.t$Type == 'Dawn')){
          # if we just have a dawn measurement
          didx <- light.t$yday[1]
          sr <- light.t$daymins[which(light.t$Type == 'Dawn')]
          
          if(sr < min.sr[didx] | sr > max.sr[didx]){
            sr <- NA
          }
          
          light[which(lightDates %in% dateVec[t] & light$Type == 'Dawn'), 6] <- sr
          
        } else if(length(light.t[,1]) == 1 & any(light.t$Type == 'Dusk')){
          # if we just have a dusk measurement
          didx <- light.t$yday[1]
          ss <- light.t$daymins[which(light.t$Type == 'Dusk')]
          
          if(ss < min.ss[didx] | ss > max.ss[didx]){
            ss <- NA
          }
          
          light[which(lightDates %in% dateVec[t] & light$Type == 'Dusk'), 6] <- ss
          
        } else{
          # if we have both dawn and dusk measurements
          didx <- light.t$yday[1]
          sr <- light.t$daymins[which(light.t$Type == 'Dawn')]
          ss <- light.t$daymins[which(light.t$Type == 'Dusk')]
          
          if(length(sr) > 1){
            # we want the first SR time if there are multiple
            sr <- sr[1]
          }
          
          if(length(ss) > 1){
            # we want the last SS time if there are multiple
            ss <- ss[length(ss)]
          }
          
          # filter based on possible grid values
          if(sr < min.sr[didx] | sr > max.sr[didx]){
            sr <- NA
          }
          
          if(ss < min.ss[didx] | ss > max.ss[didx]){
            ss <- NA
          }
          
          light[which(lightDates %in% dateVec[t] & light$Type == 'Dawn'), 6] <- sr
          light[which(lightDates %in% dateVec[t] & light$Type == 'Dusk'), 6] <- ss
          
          # now for likelihood
          # get the SD for this day, T
          srf <- raster::focal(sr.ras[[didx]], w = matrix(1, nrow = 3, ncol = 3), fun = function(x) sd(x, na.rm = T))
          # the SR likelihood
          srlik <- liksrss(sr, srss = sr.ras[[didx]], srsd = srf)
          
          # and sunset
          ssf <- raster::focal(ss.ras[[didx]], w = matrix(1, nrow = 3, ncol = 3), fun = function(x) sd(x, na.rm = T))
          sslik <- liksrss(ss, srss = ss.ras[[didx]], srsd = ssf)
          
          r <- srlik * sslik
          max.lat <- xyFromCell(r, which.max(r))[2]
          cds <- rbind(c(min(lon), max.lat), c(max(lon), max.lat))#, c(40,5), c(15,-45), c(-10,-25))
          lines <- SpatialLines(list(Lines(list(Line(cds)), "1")))
          r[] <- c(unlist(extract(r, lines)))
          
          L.locs[[t]] <- r
          
        }
        
      }
      

    } # end for loop
    
  } else{
    # if light is null, we create an empty light likelihood array?
  }

  
  ### NOW WE HAVE CLEANED LIGHT DATA. NEED TO ADD LIKELIHOOD CALC
  
} # end function





  #=============
  
  # after all that we add any known locations to overwrite any light data for a given T...
  
  if(!is.null(locs)){
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
  }
  
  for(t in 2:(length(dateVec)) - 1){
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
        
      } else if(!is.null(light) & dateVec[t] %in% lightDates){
        
        didx <- light$yday[t]
        light.t <- light[which(lightDates %in% dateVec[t]),]
        # get the SD for this day, T
        f1 <- raster::focal(sr.ras[[didx]], w = matrix(1, nrow = 3, ncol = 3), fun = function(x) sd(x, na.rm = T))
        # what's the SR?
        sr <- light.t$daymins[which(light.t$Type == 'Dawn')]
        # the SR likelihood
        srlik <- liksrss(sr, srss = sr.ras[[didx]], srsd = f1)
        
        # and sunset
        f2 <- raster::focal(ss.ras[[didx]], w = matrix(1, nrow = 3, ncol = 3), fun = function(x) sd(x, na.rm = T))
        sslik <- liksrss(ss, srss = ss.ras[[didx]], srsd = f2)
        
        idx <- which(dateVec == as.Date(time))
        L.locs[,,t] <- srlik * sslik
        
      } else{
        
        stop('Error: Error ellipse unspecified.')
        
      }
      
    } else{
      
      # no data so we skip this day
      
    }
    
  }
  
  
  # add tag/pop locations as known
  ilo <- which.min(abs(locs.grid$lon[1,] - iniloc$lon[1]))
  ila <- which.min(abs(locs.grid$lat[,1] - iniloc$lat[1]))
  L.locs[ilo, ila, 1] <- 1   # Initial location is known
  
  elo <- which.min(abs(locs.grid$lon[1,] - iniloc$lon[2]))
  ela <- which.min(abs(locs.grid$lat[,1] - iniloc$lat[2]))
  L.locs[elo, ela, length(dateVec)] <- 1  # End location is known
  
  # this performs some transformations to the likelihood array to convert to useable raster
  crs <- "+proj=longlat +datum=WGS84 +ellps=WGS84"
  list.locs <- list(x = locs.grid$lon[1,], y = locs.grid$lat[,1], z = L.locs)
  ex <- raster::extent(list.locs)
  L.locs <- raster::brick(list.locs$z, xmn = ex[1], xmx = ex[2], ymn = ex[3], ymx = ex[4], transpose = T, crs)
  L.locs <- raster::flip(L.locs, direction = 'y')
  
  print(paste('Light calculations took ', Sys.time() - start.t, '.'))
  
  return(list(L.locs = L.locs, gpsIdx = which(dateVec %in% gpsDates)))
  
}
