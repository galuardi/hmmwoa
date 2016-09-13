#' Read and format tag data
#' 
#' \code{get.env} 
#' 
#' @param 
#'   
#' @return a list containing: 
#' 
#' @export
#' 
#' @examples
#' none

get.env <- function(uniqueDates, type = NA, spatLim, save.dir = getwd()){
  
  if(is.na(type)){
    
    stop('Type of environmental data desired not specified.')
    
  } else if(type == 'sst'){
    
    for(i in 1:length(uniqueDates)){
      time <- as.Date(uniqueDates[i])
      repeat{
        get.oi.sst(spatLim, time, filename = paste(ptt, '_', time, '.nc', sep = ''), download.file = TRUE, dir = save.dir) # filenames based on dates from above
        tryCatch({
          err <- try(ncdf::open.ncdf(paste(sst.dir, ptt, '_', time, '.nc', sep = '')), silent = T)
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
          err <- try(ncdf::open.ncdf(paste(ohc.dir, ptt, '_', time, '.nc', sep = '')), silent = T)
        }, error=function(e){print(paste('ERROR: Download of data at ', time, ' failed. Trying call to server again.', sep = ''))})
        if(class(err) != 'try-error') break
      }
    }
    
  }
  
  
}