setup.grid.raster <- function(grid.raster){
  
  ex <- raster::extent(grid.raster)
  
  # Find longitude extents
  il <- floor(ex[1])
  al <- ceiling(ex[2])
  lx <- 0.1 * (al - il)
  lonl <- il - lx
  lonu <- al + lx
  
  # Find latitude extents
  ila <- floor(ex[3])
  ala <- ceiling(ex[4])
  ly <- 0.1 * (ala - ila)
  latl <- ila - ly
  latu <- ala + ly
  
  # Create grid
  lo <- raster::xFromCol(grid.raster)
  la <- raster::yFromRow(grid.raster)
  g <- meshgrid(lo, la)
  dlo <- raster::xres(grid.raster)
  dla <- raster::yres(grid.raster)
  
  list(lon = g$X, lat = g$Y, dlo = dlo, dla = dla)
  
}