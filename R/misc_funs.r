setup.grid <- function(locations){
  ## Setup the discrete spatial grid for the HMM
#' @param locations is dataframe of -Locations from WC psat tag
  
  T <- length(locations$Longitude)
  
  # Find longitude extents
  il <- floor(min(locations$Longitude))
  al <- ceiling(max(locations$Longitude))
  lx <- 0.1 * (al - il)
  lonl <- il - lx
  lonu <- al + lx
  lo.out <- 25/2 * (lonu - lonl)
  
  # Find latitude extents
  ila <- floor(min(locations$Latitude))
  ala <- ceiling(max(locations$Latitude))
  ly <- 0.1 * (ala - ila)
  latl <- ila - ly
  latu <- ala + ly
  la.out <- (25/2) * (latu - latl)
  
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

meshgrid <- function(x, y){
  # funtion from Pedersen et al 2011
  Y <- repmat(as.matrix(y), 1, length(x))
  X <- repmat(t(as.matrix(x)), length(y), 1)
  list(X = X, Y = Y)
}

repmat <- function(X, m, n){
  # function from Pedersen et al 2011
  mx = dim(X)[1]
  nx = dim(X)[2]
  matrix(t(matrix(X, mx, nx * n)), mx * m, nx * n, byrow = T)
}