calc.locs <- function(locs, iniloc, g, raster = TRUE, dateVec){
  ## Calculate the "data" likelihood, i.e. the likelihood field for each
  ## location observation

  #' @param: locs is -Locations file output from DAP for WC tags and contains
  #'        GPS, Argos, and GPE locations as applicable
  #' @param: iniloc is 2 x 5 dataframe containing day, month, year,
  #'          lat, lon for both tag and pop locations
  #' @param: g is output from setup.grid and indicates extent and resolution
  #'        of grid used to calculate likelihoods
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
  # SD for light-based longitude from Musyl et al. (2001)
  sl.sd <- 35/111 # Converting from kms to degrees
  
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
     # create longitude likelihood based on GPE data
     # for now, latitude is ignored
      L.light <- dnorm(t(g$lon), locs$Longitude[t], sl.sd) # Longitude data
      L.locs[,,which(dateVec == locDates[t])] <- L.light #(L.light / max(L.light, na.rm = T)) - .05
      
    } else{}

      }
  
  # End location is known
  elo <- which.min(abs(g$lon[1,]-iniloc$lon[2]))
  ela <- which.min(abs(g$lat[,1]-iniloc$lat[2]))
  L.locs[elo, ela, length(dateVec)] <- 1
  
  if(raster){
    crs <- "+proj=longlat +datum=WGS84 +ellps=WGS84"
    list.locs <- list(x = g$lon[1,], y = g$lat[,1], z = L.locs)
    ex <- extent(list.locs)
    L.locs <- brick(list.locs$z, xmn=ex[1], xmx=ex[2], ymn=ex[3], ymx=ex[4], transpose=T, crs)
    L.locs <- flip(L.locs, direction = 'y')
  }
  
  print(class(L.locs))
  return(L.locs)
  
}