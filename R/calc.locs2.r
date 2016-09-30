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
          
          # now for likelihood
          # get the SD for this day, T
          srf <- raster::focal(sr.ras[[didx]], w = matrix(1, nrow = 3, ncol = 3), fun = function(x) sd(x, na.rm = T))
          # the SR likelihood
          srlik <- liksrss(sr, srss = sr.ras[[didx]], srsd = srf)
          
          L.locs[[t]] <- srlik
          
        } else if(length(light.t[,1]) == 1 & any(light.t$Type == 'Dusk')){
          # if we just have a dusk measurement
          didx <- light.t$yday[1]
          ss <- light.t$daymins[which(light.t$Type == 'Dusk')]
          
          if(ss < min.ss[didx] | ss > max.ss[didx]){
            ss <- NA
          }
          
          light[which(lightDates %in% dateVec[t] & light$Type == 'Dusk'), 6] <- ss
          
          # and sunset
          ssf <- raster::focal(ss.ras[[didx]], w = matrix(1, nrow = 3, ncol = 3), fun = function(x) sd(x, na.rm = T))
          sslik <- liksrss(ss, srss = ss.ras[[didx]], srsd = ssf)
          
          L.locs[[t]] <- sslik
          
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

} # end function
