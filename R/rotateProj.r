rotateProj <- function(spobj, angle) {
  
  ## Function written to rotate a projected spatial dataset
  ## Accessed here: https://rpubs.com/geospacedman/rotatespatial
  ## on 7 Dec 2015
  
  #' @param spobj is a spatial object
  #' @param angle is rotation angle clockwise from 0-360
  #' @return crs used to re-project (spTransform) your spatial
  #'        object resulting in rotation
  
  # get bounding box as spatial points object
  boxpts = SpatialPoints(t(bbox(spobj)), proj4string = CRS(proj4string(spobj)))
  
  # convert to lat-long
  boxLL = bbox(spTransform(boxpts, CRS("+init=epsg:4326")))
  
  # find the centre
  llc = apply(boxLL, 1, mean)
  
  # construct the proj4 string
  prj = paste0("+proj=omerc +lat_0=", llc[2], " +lonc=", llc[1], " +alpha=", 
               angle, " +gamma=0.0 +k=1.000000 +x_0=0.000 +y_0=0.000 +ellps=WGS84 +units=m ")
  
  # return as a CRS:
  CRS(prj)
}