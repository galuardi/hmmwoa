#' Read and format tag data
#' 
#' \code{read.wc} reads and formats tag data output from Wildlife Computers Data Portal
#' 
#' @param 
#'   
#' @return a list containing: 
#' 
#' @export
#' 
#' @examples
#' none


read.wc <- function(ptt, wd = getwd(), tag, pop, gpeNo=NULL, type = 'sst'){
  
  if(type == 'pdt'){
    # READ IN PDT DATA FROM WC FILES
    data <- read.table(paste(wd, ptt,'-PDTs.csv', sep=''), sep=',',header=T,blank.lines.skip=F, skip = 0)
    data <- extract.pdt(data)
    dts <- as.POSIXct(data$Date, format = findDateFormat(data$Date))
    d1 <- as.POSIXct('1900-01-02') - as.POSIXct('1900-01-01')
    didx <- dts >= (tag + d1) & dts <= (pop - d1)
    data <- data[didx,]
    udates <- unique(as.Date(data$Date))
    
  } else if(type == 'sst'){
    # READ IN TAG SST FROM WC FILES
    data <- read.table(paste(wd, ptt, '-SST.csv', sep=''), sep=',',header=T, blank.lines.skip=F)
    dts <- as.POSIXct(data$Date, format = findDateFormat(data$Date))
    d1 <- as.POSIXct('1900-01-02') - as.POSIXct('1900-01-01')
    didx <- dts >= (tag + d1) & dts <= (pop - d1)
    data <- data[didx,]
    if (length(data[,1]) <= 1){
      stop('Something wrong with reading and formatting of tags SST data. Check date format.')
    }
    dts <- as.POSIXct(data$Date, format = findDateFormat(data$Date))
    udates <- unique(as.Date(dts))
    
  } else if(type == 'light'){
    # READ IN LIGHT DATA FROM WC FILES
    data <- read.table(paste(wd, ptt, '-LightLoc.csv', sep=''), sep=',',header=T, blank.lines.skip=F,skip=2)
    data <- data[which(!is.na(data[,1])),]
    dts <- as.POSIXct(data$Day, format = '%d-%b-%Y', tz = 'UTC')
    d1 <- as.POSIXct('1900-01-02') - as.POSIXct('1900-01-01')
    didx <- (dts > (tag + d1)) & (dts < (pop - d1))
    data <- data[didx,]
    udates <- unique(as.Date(dts))
    
  }
  
  return(list(data = data, udates = udates))
           
}

