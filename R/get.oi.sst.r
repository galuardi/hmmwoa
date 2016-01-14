
get.oi.sst <- function(lon, lat, time, filename='', download.file=TRUE, dir = getwd()) {

  #' Downloads data from the HYCOM + NCODA Global 1/12 Degree Analysis.
  #'
  #' The method may return before the download is completed. It will continue
  #' to display progress  until the download completes. 
  #' 
  #' Ideally download.file (default method) would be used instead of curl (optional), but this does not 
  #' seem to work on some platforms.
  #'
  #' @param lon An vector of length 2 with the minimum and maximum longitude. 
  #' @param lat An vector of length 2 with the minimum and maximum latitude.
  #' @param time An vector of length 2 with the minimum and maximum times.
  #' @param vars A list of variables to download. This should only contain
  #' 'sst' 'anom' 'err' but is not checked for errors
  #' @param include_latlon Should the array of latitude and longitude values be
  #' included?
  #' @param filename An optional filename. If provided, then the data is
  #' downloaded to that file. Otherwise the data is not downloaded.
  #' @param download.file Should use the default download.file function to query 
  #' the server and download or use the optional curl function. Some users may
  #' need to use curl in order to get this to work.
  #' @param: dir is directory where nc files should be downloaded to. default is
  #' current working directory. if enter a directory that doesn't exist, it will
  #' be created.
  #' @param type indicates type of HYCOM product to download. 'r' is reanalysis
  #'        and 'a' is analysis. see https://hycom.org/dataserver for details.
  #' @return The url used to extract the requested data from the NetCDF subset
  #' service.
    
  ## Function originally written for R by Ben Jones (WHOI) and modified by Camrin
  ## Braun and Ben Galuardi.
  
  require(ncdf)
  
  dir.create(file.path(dir), recursive = TRUE)
  setwd(dir)
  
  ## Set the base URL based on the start date. If the ending date exceeds the
  ## period for this experiment, then print a warning and truncate the output
  ## early.
  
  
    expts = data.frame(
      start=c(as.Date('1981-09-01')),
      end=c(Sys.Date() + 1),
      url=c('http://coastwatch.pfeg.noaa.gov/erddap/griddap/jplL4AvhrrOIv1fv2.nc?analysed_sst'))
  
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
  
  ## Add the time domain.
    if(length(time) == 2){
      url = sprintf('%s[(%s:00:00Z):1:(%s:00:00Z)]',
                    url, strftime(time[1], '%Y-%m-%dT00'),
                    strftime(time[2], '%Y-%m-%dT00'))
    } else if(length(time)==1){
      url = sprintf('%s[(%s:00:00Z):1:(%s:00:00Z)]',
                    url, strftime(time[1], '%Y-%m-%dT00'),
                    strftime(time[1], '%Y-%m-%dT00'))
    }
    
  ## Add the spatial domain.
  url = sprintf('%s[(%s):1:(%s)][(%s):1:(%s)]',
                url, lat[1], lat[2], lon[1], lon[2])

  ## Download the data if a filename was provided.
  if(filename != ''){
    if(download.file == TRUE){
      download.file(url, filename, method = 'auto')
    } else if(download.file == FALSE){
      system(sprintf('curl -o "%s" "%s"', filename, url))
    }
  }
  return(url)
}


