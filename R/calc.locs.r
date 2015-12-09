calc.locs <- function(locs, iniloc, g, raster = TRUE, dateVec, errEll = F){
  ## Calculate the "data" likelihood, i.e. the likelihood field for each
  ## location observation

  #' @param: locs is -Locations file output from DAP for WC tags and contains
  #'        GPS, Argos, and GPE locations as applicable
  #' @param: iniloc is 2 x 5 dataframe containing day, month, year,
  #'          lat, lon for both tag and pop locations
  #' @param: g is output from setup.grid and indicates extent and resolution
  #'        of grid used to calculate likelihoods
  #' @param errEll is logical indicating whether error ellipses should be 
  #'        generated for light-based likelihoods as given from output
  #'        of WC-GPE. False if only longitude should be used.
  #' @return: L is array of lon x lat likelihood surfaces (matrices)
  #'          for each time point (3rd dimension)
    
  T <- length(locs$Longitude)
  row <- dim(g$lon)[1]
  col <- dim(g$lon)[2]
  
  L.locs <- array(0, dim = c(col, row, length(dateVec)))
  
  # Initial location is known
  ilo <- which.min(abs(g$lon[1,]-iniloc$lon[1]))
  ila <- which.min(abs(g$lat[,1]-iniloc$lat[1]))
  L.locs[ilo, ila, 1] <- 1
  
  # Calculate data likelihood
   
  locDates <- as.Date(locs$Date, format = findDateFormat(locs$Date))
  
  for(t in 1:T){
    if(locs$Type[t] == 'GPS'){
      # if GPS exists then other forms of data for that time point are obsolete
      glo <- which.min(abs(g$lon[1,]-locs$Longitude[t]))
      gla <- which.min(abs(g$lat[,1]-locs$Latitude[t]))
      L.locs[glo, gla, which(dateVec == locDates[t])] <- 1
      
    } else if(locs$Type[t] == 'Argos'){
      # if Argos exists, GPE positions are obsolete
      alo <- which.min(abs(g$lon[1,]-locs$Longitude[t]))
      ala <- which.min(abs(g$lat[,1]-locs$Latitude[t]))
      L.locs[alo, ala, which(dateVec == locDates[t])] <- 1
      
    } else if(locs$Type[t] == 'GPE'){
     if(errEll == FALSE){
       # create longitude likelihood based on GPE data
       # for now, latitude is ignored
       # SD for light-based longitude from Musyl et al. (2001)
       slon.sd <- 35/111 # Converting from kms to degrees
       
       L.light <- dnorm(t(g$lon), locs$Longitude[t], slon.sd)
       
       L.locs[,,which(dateVec == locDates[t])] <- L.light
       
     } else if(errEll == TRUE){
       stop('Error: Error ellipse functionality is not yet available.')
       #slon.sd <- locs$Error.Semi.minor.axis[t] / 1000 / 111 #semi minor axis
       #L.light.lon <- dnorm(t(g$lon), locs$Longitude[t], slon.sd) # Longitude data
       #slat.sd <- locs$Error.Semi.major.axis[t] / 1000 / 111 #semi major axis
       #L.light.lat <- dnorm(t(g$lat), locs$Latitude[t], slat.sd)
       
       #L <- raster(L.light.lat * L.light.lon, xmn=min(lon), 
      #             xmx=max(lon),ymn=min(lat),ymx=max(lat))
       # offset
       #Lext <- extent(L)
       #Ls <- shift(L, y = -1 * (locs$Offset[t] / 1000 / 111))
       #Lsx <- extend(Ls, Lext)
       #Lsx <- crop(Lsx, Lext)
       ###Lsx <- raster(matrix(extract(Ls, Lext), nrow=nrow(L), ncol=ncol(L)), xmn=min(lon), 
       ###              xmx=max(lon),ymn=min(lat),ymx=max(lat))
       
       # rotate?
       #Lsx <- rotateProj(Lsx, locs$Error.Ellipse.orientation[t])
       
       ###xres <- lon[2] - lon[1]; yres <- lat[2] - lat[1]
       ###locs$Offset[t] / 1000 / 111 # offset in degrees
       
       #L.locs[,,which(dateVec == locDates[t])] <- L.light.lat * L.light.lon
       
     }
     
    } else{}

      }
  
  # End location is known
  elo <- which.min(abs(g$lon[1,]-iniloc$lon[2]))
  ela <- which.min(abs(g$lat[,1]-iniloc$lat[2]))
  L.locs[elo, ela, length(dateVec)] <- 1
  
  crs <- "+proj=longlat +datum=WGS84 +ellps=WGS84"
  list.locs <- list(x = g$lon[1,], y = g$lat[,1], z = L.locs)
  ex <- extent(list.locs)
  
  if(raster == 'brick'){
    L.locs <- brick(list.locs$z, xmn=ex[1], xmx=ex[2], ymn=ex[3], ymx=ex[4], transpose=T, crs)
    L.locs <- flip(L.locs, direction = 'y')
    s <- L.locs
    
  } else if(raster == 'stack'){
    for(ii in 1:length(dateVec)){
      r <- raster(t(L.locs[,,ii]), xmn=ex[1], xmx=ex[2], ymn=ex[3], ymx=ex[4], crs)
      r <- flip(r, direction = 'y')
      if(ii == 1){
        s <- stack(r, quick = T)
      } else{
        s <- stack(s, r, quick = T)
      }
    }
  } else if(raster == 'array'){s <- L.locs}
  
  
  print(class(s))
  return(s)
  
}