#' Extract temperatures from World Ocean Atlas
#' 
#' \code{extract.woa} extracts the desired temperature data from a global 
#' dataset derived from monthly gridded climatology data contained in the 2013 
#' World Ocean Atlas
#' 
#' @param nc.dir directory to load the global nc file from; specify the complete
#'   path to the nc file unless it is in your current working directory
#' @param bbox bounding box of form list(long min, long max, lat min, lat max)
#' @param resolution indicates whether oceanographic data is gridded at 
#'   'quarter' or 'one' degree resolution
#'   
#' @return a list containing: DAT is an array of temperature data with 
#'   dimensions (long, lat, depth, time) depth contains 57 standard depth levels
#'   by default and levels are defined in variable 'depth' contained here. time 
#'   is monthly and spans the entire year. LON/LAT are vectors of lon/lat bounds
#' @export
#' 
#' @examples
#' none

extract.woa <- function(nc.dir, bbox=NULL, resolution){
  
  # load global nc
  nc = ncdf::open.ncdf(nc.dir)
  
  # retrieve var bounds from global nc
  lon = ncdf::get.var.ncdf(nc, 'Longitude')
  lat = ncdf::get.var.ncdf(nc, 'Latitude')
  depth = ncdf::get.var.ncdf(nc, 'Depth')
  
  if(is.null(bbox) & resolution == 'quarter'){
    bbox <- list(lonmin = -180, lonmax = 179.75, latmin = -90, latmax = 89.75)
  } else if(is.null(bbox) & resolution == 'one'){
    bbox <- list(lonmin = -179.5, lonmax = 179.5, latmin = -89.5, latmax = 89.5)
  }
  
  # set bounds for extracting data
  xmin = which.min((bbox[[1]] - lon) ^ 2)
  xmax = which.min((bbox[[2]] - lon) ^ 2)
  ymin = which.min((bbox[[3]] - lat) ^ 2) 
  ymax = which.min((bbox[[4]] - lat) ^ 2)

  if(resolution == 'quarter'){
    xlen = floor(4 * (bbox[[2]] - bbox[[1]])) # for quarter degree
    ylen = floor(4 * (bbox[[4]] - bbox[[3]])) 
  } else if(resolution == 'one'){
    xlen = bbox[[2]] - bbox[[1]] # for one degree
    ylen = bbox[[4]] - bbox[[3]]
  } else{
    stop('Resolution of input oceanographic data not defined.')
  }  
  
  dat = ncdf::get.var.ncdf(nc, start = c(xmin, ymin, 1, 1), count = c(xlen + 1, ylen + 1, 57, 12))
  
  returnWOA = list(dat = dat, lon = lon[xmin:xmax], lat = lat[ymin:ymax], depth = depth)
  
}
