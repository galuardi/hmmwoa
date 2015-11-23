
get.hycom = function(lon, lat, time, vars=c('temperature'), include_latlon=TRUE,
                          filename='',download.file=TRUE, dir = getwd()) {
  #' Downloads data from the HYCOM + NCODA Global 1/12 Degree Analysis.
  #'
  #' The method may return before the download is completed. It will continue
  #' to display a progress bar until the download completes. 
  #' 
  #' Ideally download.file (default method) would be used instead of curl (optional), but this does not 
  #' seem to work on some platforms.
  #'
  #' @param lon An vector of length 2 with the minimum and maximum longitude. 
  #' @param lat An vector of length 2 with the minimum and maximum latitude.
  #' @param time An vector of length 2 with the minimum and maximum times.
  #' @param vars A list of variables to download. This should only contain
  #' 'emp', 'mld', 'mlp', qtot', 'ssh', 'surface_salinity_trend',
  #' 'surface_temperature_trend', 'salinity', 'temperature', 'u', and 'v', but is
  #' not checked for errors.
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
  #' @return The url used to extract the requested data from the NetCDF subset
  #' service.
    
  ## Function originally written for R by Ben Jones (WHOI) and modified by Camrin
  ## Braun and Ben Galuardi.
  require(ncdf)
  
  dir.create(file.path(dir),recursive=TRUE)
  setwd(dir)
  
  ## Set the base URL based on the start date. If the ending date exceeds the
  ## period for this experiment, then print a warning and truncate the output
  ## early.
  expts = data.frame(
    start=c(as.Date('2008-09-19'), as.Date('2009-05-07'),
            as.Date('2011-01-03'), as.Date('2013-08-21'),
            as.Date('2014-04-05')),
    end=c(as.Date('2009-05-06'), as.Date('2011-01-02'),
          as.Date('2013-08-20'), as.Date('2014-04-04'),
          Sys.Date() + 1),
    url=c('http://ncss.hycom.org/thredds/ncss/GLBa0.08/expt_90.6?',
          'http://ncss.hycom.org/thredds/ncss/GLBa0.08/expt_90.8?',
          'http://ncss.hycom.org/thredds/ncss/GLBa0.08/expt_90.9?',
          'http://ncss.hycom.org/thredds/ncss/GLBa0.08/expt_91.0?',
          'http://ncss.hycom.org/thredds/ncss/GLBa0.08/expt_91.1?'))
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
  ## Add the variables.
  for(var in vars)
    url = sprintf('%svar=%s&', url, var)
  ## Add the spatial domain.
  url = sprintf('%snorth=%f&west=%f&east=%f&south=%f&horizStride=1&',
                url, lat[2], lon[1], lon[2], lat[1])
  ## Add the time domain.
  if(length(time)==2){
    url = sprintf('%stime_start=%s%%3A00%%3A00Z&time_end=%s%%3A00%%3A00Z&timeStride=1&',
                  url, strftime(time[1], '%Y-%m-%dT00', start_time),
                  strftime(time[2], '%Y-%m-%dT00', end_time))
  } else if(length(time)==1){
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
    if(download.file==TRUE){
      download.file(url,filename)
    } else if(download.file==FALSE){
      system(sprintf('curl -o "%s" "%s"', filename, url))
    }
  }
  return(url)
}

