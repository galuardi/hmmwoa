#' Determine date format of vector
#' 
#' \code{findDateFormat} determines the date format of a given vector of dates
#' 
#' @param dateVec a vector of dates as in those from -Histos or  -PDTs from WC tags
#' @return dateformat is character string used as input to strptime(format = dateformat)
#' @export 
#'
#' @examples
#' none
#' 
findDateFormat <- function(dateVec){

  dateformat = '%Y-%m-%d %H:%M:%S'
  ddates = as.POSIXct(strptime(as.character(dateVec), format = dateformat)) #reads dates as dates

  if (is.na(ddates[1])){
    dateformat = '%H:%M:%S %d-%b-%Y'
    ddates = as.POSIXct(strptime(as.character(dateVec), format = dateformat)) #reads dates as dates
  
    if (is.na(ddates[1])){
      dateformat = '%m-%d-%y %H:%M'
      ddates = as.POSIXct(strptime(as.character(dateVec), format = dateformat)) #reads dates as dates
      
      if (is.na(ddates[1])){
        dateformat = '%m/%d/%y %H:%M'
        ddates = as.POSIXct(strptime(as.character(dateVec), format = dateformat)) #reads dates as dates
        
        if (is.na(ddates[1])){
          dateformat = '%m/%d/%Y'
          ddates = as.POSIXct(strptime(as.character(dateVec), format = dateformat)) #reads dates as dates
          
          if (is.na(ddates[1])){
            dateformat = '%H:%M:%S %d-%b-%Y'
            ddates = as.POSIXct(strptime(as.character(dateVec), format = dateformat)) #reads dates as dates
            
            if (is.na(ddates[1])){
              dateformat = '%H:%M:%OS %d-%b-%Y'
              ddates = as.POSIXct(strptime(as.character(dateVec), format = dateformat)) #reads dates as dates
              
              if(is.na(ddates[1])){
                stop('No correct date format was found.')
              }
              
              } else {}
            } else {}
          } else {}       
      } else {}
    } else {}
  } else {}
  
  dateformat #return dateformat variable

}
