#' Download and Read Oceanographic Data
#' 
#' \code{get.env} accesses oceanographic data like sea surface temperature from a remote server and downloads the temporal and spatial extent of interest for further use
#' 
#' @param uniqueDates is a POSIXct vector of desired dates
#' @param type is a character string indicating whether you're after sea surface temperature 'sst' or hybrid coordinate ocean model 'ohc' data
#' @param spatLim is a list of spatial limits as \code{list(xmin, xmax, ymin, ymax)}
#' @param save.dir is the directory to save the downloaded data to
#'   
#' @return nothing, just downloads the data to your local machine
#' @export
#' @examples
#' splim <- list(xmin=-60, xmax=-40, ymin=10, ymax=30)
#' get.env('2015-01-01', type = 'sst', splim, getwd())

get.env <- function(uniqueDates, ptt, type = NA, spatLim, save.dir = getwd()){
  
  if(is.na(type)){
    
    stop('Type of environmental data desired not specified.')
    
  } else if(type == 'sst'){
    
    for(i in 1:length(uniqueDates)){
      time <- as.Date(uniqueDates[i])
      repeat{
        get.oi.sst(spatLim, time, filename = paste(ptt, '_', time, '.nc', sep = ''), download.file = TRUE, dir = save.dir) # filenames based on dates from above
        tryCatch({
          err <- try(RNetCDF::open.nc(paste(save.dir,'/', ptt, '_', time, '.nc', sep = '')), silent = T)
        }, error=function(e){print(paste('ERROR: Download of data at ', time, ' failed. Trying call to server again.', sep = ''))})
        if(class(err) != 'try-error') break
      }
    }

  } else if(type == 'ohc'){
    
    for(i in 1:length(uniqueDates)){
      time <- as.Date(uniqueDates[i])
      repeat{
        get.hycom(spatLim, time, type = 'a', filename = paste(ptt, '_', time, '.nc', sep = ''),
                  download.file = TRUE, dir = save.dir, vars = 'water_temp') 
        tryCatch({
          err <- try(RNetCDF::open.nc(paste(save.dir,'/', ptt, '_', time, '.nc', sep = '')), silent = T)
        }, error=function(e){print(paste('ERROR: Download of data at ', time, ' failed. Trying call to server again.', sep = ''))})
        if(class(err) != 'try-error') break
      }
    }
    
  }
  
  
}