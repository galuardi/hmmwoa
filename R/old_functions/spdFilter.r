
#bathy <- raster('~/Documents/WHOI/RData/Bathymetry/ETOPO1_Ice_g_gmt4.grd',
#            crs="+proj=longlat +ellps=WGS84 +datum=WGS84")



spotFilter <- function(spot, dts, speed=7.2, bathy){
  # QC SPOT tracks
  
  # spot is data.frame of spot locations
  # dts is POSIXct object containing date/times for each position in spot
  # speed is maximum speed allowed, used to filter erroneous positions
  
  require(trip); require(raster)
  
  # filter argos location classes

  
  # remove duplicate date-times
  spot$dts <- dts
  spot <- spot[which(!duplicated(spot$dts)),] #modify this to take better loc class from duplicate times
  
  # convert track to SpatialPointsDataFrame
  #coordinates(spot) <- ~Longitude + Latitude
  coordinates(spot) <- ~lon + lat
  proj4string(spot)=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84")
  
  # filter based on bathymetry - use raster::extract to sample
  #   bathy at each point and filter >=0
  ext <- raster::extract(bathy, spot) < 0
  spot <- spot[ext,]
  
  # then to trip object
  tr <- trip(spot, c('dts','id'))
  
  # filter at maximum speed in km/hr
  # e.g. 2 m/s = 7.2 km/hr
  sf <- speedfilter(tr, max.speed = speed)
  
  # and subset track based on logical output of speed filter
  spot <- data.frame(tr[sf,])#[,c(1:5)]

}
