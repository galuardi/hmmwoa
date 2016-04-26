setup.grid <- function(locations, res){
  ## Setup the discrete spatial grid for the HMM
  #' @param locations is dataframe of -Locations from WC psat tag
  #' @param res is character indicating resolution of grid. 'hycom' is
  #'        .08 to match hycom reanalysis res. 'quarter' and 'one' are
  #'        .25 and 1 deg, respectively.
  
  T <- length(locations$Longitude)
  
  # Find longitude extents
  il <- floor(min(locations$Longitude))
  al <- ceiling(max(locations$Longitude))
  lx <- 0.1 * (al - il)
  lonl <- il - lx
  lonu <- al + lx
  
  if(res == 'hycom'){
    lo.out <- 25/2 * (lonu - lonl)
  } else if(res == 'quarter'){
    lo.out <- 4 * (lonu - lonl)
  } else if(res == 'one'){
    lo.out <- 1 * (lonu - lonl)
  }
  
  # Find latitude extents
  ila <- floor(min(locations$Latitude))
  ala <- ceiling(max(locations$Latitude))
  ly <- 0.1 * (ala - ila)
  latl <- ila - ly
  latu <- ala + ly
  
  if(res == 'hycom'){
    la.out <- 25/2 * (latu - latl) # 0.08 deg
  } else if(res == 'quarter'){
    la.out <- 4 * (latu - latl) # .25 deg
  } else if(res == 'one'){
    la.out <- 1 * (latu - latl) # 1 deg
  }
  
  #  latvec <- seq(0,90)
  #  lats <- rep(0,T)
  #  for(t in 1:T){
  #    time <- date2time(lsst$date[t])
  #    #time <- as.numeric(strftime(lsst$date[t],format='%j'))
  #    ssts <- sstdb(time,lsst$lon[t],latvec)
  #    lats[t] <- latvec[sum(lsst$sst[t]<ssts)]
  #  }
  #  lx <- 0.1*(max(lats)-min(lats))
  #  latl <- min(lats) - lx
  #  latu <- max(lats) + lx
  
  # Create grid
  lo <- seq(lonl, lonu, length.out = lo.out)
  la <- seq(latl, latu, length.out = la.out)
  g <- meshgrid(lo, la)
  dlo <- lo[2] - lo[1]
  dla <- la[2] - la[1]
  
  list(lon = g$X, lat = g$Y, dlo = dlo, dla = dla)
}
