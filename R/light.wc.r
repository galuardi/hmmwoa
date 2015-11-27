light.wc <- function(light,iniloc,g){
  ## Calculate the "data" likelihood, i.e. the likelihood field for each observation
  
  T <- length(light$Longitude)
  row <- dim(g$lon)[1]
  col <- dim(g$lon)[2]
  
  L <- array(0,dim=c(col, row, T + 2))
  
  # Initial location is known
  ilo <- which.min(abs(g$lon[1,]-iniloc$lon[1]))
  ila <- which.min(abs(g$lat[,1]-iniloc$lat[1]))
  L[ilo, ila, 1] <- 1
  
  # Calculate data likelihood
  # SD for light-based longitude from Musyl et al. (2001)
  sl.sd <- 35/111 # Converting from kms to degrees
  
  for(t in 1:T){
    #time <- date2time(as.POSIXct(light$Date[t], format = findDateFormat(light$Date)))
    L[,, (t + 1)] <- dnorm(t(g$lon), light$Longitude[t], sl.sd) # Longitude data
    #L[,,t] <- Lsst*Llon
  }
  
  # End location is known
  elo <- which.min(abs(g$lon[1,]-iniloc$lon[2]))
  ela <- which.min(abs(g$lat[,1]-iniloc$lat[2]))
  L[elo, ela, T + 2] <- 1
  
  L
  
}