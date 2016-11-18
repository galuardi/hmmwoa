#' Download HYCOM data
#' 
#' \code{get.hycom} downloads HYbrid Coordinate Ocean Model (HYCOM) data for 
#' given temporal and spatial constraints of your data.
#' 
#' The method may return before the download is completed. It will continue to 
#' display progress  until the download completes. Change the default 
#' \code{download.file} if the download is failing on your platform.
#' 
#' @param limits A list of length 4; minlon, maxlon, minlat, maxlat. Longitude values are -180,180
#' @param time A vector of length 2 with the minimum and maximum times in form 
#'   \code{as.Date(date)}.
#' @param vars A list of variables to download. This should only contain 
#'   'surf_el', 'salinity', 'water_u', water_v','water_temp' but is not checked 
#'   for errors. We currently only use 'water_temp'.
#' @param include_latlon Should the array of latitude and longitude values be 
#'   included?
#' @param filename An optional filename. If provided, then the data is 
#'   downloaded to that file. Otherwise the data is not downloaded and the url 
#'   is returned.
#' @param download.file Logical. Should use the default \code{download.file} to 
#'   query the server and download or use the optional \code{curl}. Some users 
#'   may need to use \code{curl} in order to get this to work.
#' @param dir is local directory where ncdf files should be downloaded to. 
#'   default is current working directory. if enter a directory that doesn't 
#'   exist, it will be created.
#' @param type indicates type of HYCOM product to download. 'r' is reanalysis 
#'   and 'a' is analysis. see \url{https://hycom.org/dataserver/} for details.
#' @return The url used to extract the requested data from the NetCDF subset 
#'   service.
#' @examples 
#' lon <- c(-90, -60)
#' lat <- c(0, 30)
#' time <- as.Date('2013-03-01')
#' get.hycom(lon, lat, time, type='a', filename = '', vars = 'water_temp')
#' # only returns url because filename is unspecified
#' \dontrun{ 
#' get.hycom(lon, lat, time, type='a', filename = 'my_data.nc', vars = 'water_temp')
#' nc <- open.nc('my_data.nc')
#' hycom <- var.get.nc(nc, 'water_temp')
#' image.plot(hycom[,,1])
#' }
#'   
#' @author   Function originally written for R by Ben Jones (WHOI) and modified
#'   by Camrin Braun and Ben Galuardi.
#' @references \url{https://hycom.org/}
#'   


get.hycom <- function(limits, time, vars=c('water_temp'), include_latlon=TRUE,
                      filename='', type = 'r', download.file=TRUE, dir = getwd()) {
  
  
  
  dir.create(file.path(dir), recursive = TRUE, showWarnings = FALSE)
  setwd(dir)
  
  ## Set the base URL based on the start date. If the ending date exceeds the
  ## period for this experiment, then print a warning and truncate the output
  ## early.
  
  if(type == 'a'){
    expts = data.frame(
      start=c(as.Date('1992-10-02'), as.Date('1995-08-01'),
              as.Date('2012-01-01'), as.Date('2013-08-21'),
              as.Date('2014-04-05'), as.Date('2016-04-18')),
      end=c(as.Date('1995-07-31'), as.Date('2011-12-31'),
            as.Date('2013-08-20'), as.Date('2014-04-04'),
            as.Date('2016-04-17'), Sys.Date() + 1),
      url=c('http://ncss.hycom.org/thredds/ncss/GLBu0.08/expt_19.0?',
            'http://ncss.hycom.org/thredds/ncss/GLBu0.08/expt_19.1?',
            'http://ncss.hycom.org/thredds/ncss/GLBu0.08/expt_90.9?',
            'http://ncss.hycom.org/thredds/ncss/GLBu0.08/expt_91.0?',
            'http://ncss.hycom.org/thredds/ncss/GLBu0.08/expt_91.1?',
            'http://ncss.hycom.org/thredds/ncss/GLBu0.08/expt_91.2?'))
  } else if(type == 'r'){
    expts = data.frame(
      start=c(as.Date('1992-10-02'),
              as.Date('1995-08-01')),
      end=c(as.Date('1995-07-31'),
            as.Date('2012-12-31')),
      url=c('http://ncss.hycom.org/thredds/ncss/GLBu0.08/expt_19.0',
            'http://ncss.hycom.org/thredds/ncss/GLBu0.08/expt_19.1'))
  }
  
  if(time[1] < expts$start[1])
    stop('Data begins at %s and is not available at %s.',
         strftime(expts$start[1], '%d %b %Y'),
         strftime(time[1], '%d %b %Y'))
  if(time[1] > expts$end[nrow(expts)])
    stop('Data ends at %s and is not available at %s.',
         strftime(expts$end[nrow(expts)], '%d %b %Y'),
         strftime(time[1], '%d %b %Y'))
  for(i in seq(nrow(expts))) {
    if((time[1] >= expts$start[i]) & (time[1] <= expts$end[i]))
      url = expts$url[i]
  }
  
  if(type == 'r'){
    url = sprintf('%s/%s', url, as.numeric(format(time, '%Y')))
    url = sprintf('%s?%s', url, '')
  }
  
  ## Add the variables.
  for(var in vars)
    url = sprintf('%svar=%s&', url, var)
  ## Add the spatial domain.
  url = sprintf('%snorth=%f&west=%f&east=%f&south=%f&horizStride=1&',
                url, limits[[4]], limits[[1]], limits[[2]], limits[[3]])
  # north, west, east, south
  
  ## Add the time domain.
  if(length(time) == 2){
    url = sprintf('%stime_start=%s%%3A00%%3A00Z&time_end=%s%%3A00%%3A00Z&timeStride=1&',
                  url, strftime(time[1], '%Y-%m-%dT00', start_time),
                  strftime(time[2], '%Y-%m-%dT00', end_time))
  } else if(length(time) == 1){
    url = sprintf('%stime_start=%s%%3A00%%3A00Z&time_end=%s%%3A00%%3A00Z&timeStride=1&',
                  url, strftime(time[1], '%Y-%m-%dT00'),
                  strftime(time[1], '%Y-%m-%dT00'))
  }
  
  ## Add the lat-lon points if requested.
  if(include_latlon)
    url = sprintf('%saddLatLon=true&', url)
  ## Finish the URL.
  url = sprintf('%sdisableProjSubset=on&vertCoord=&accept=netcdf', url)
  ## Download the data if a filename was provided.
  if(filename != ''){
    if(download.file == TRUE){
      download.file(url, filename, method = 'curl')
    } else if(download.file == FALSE){
      system(sprintf('curl -o "%s" "%s"', filename, url))
    }
  }
  return(url)
}

