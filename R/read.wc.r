#' Read and format tag data
#' 
#' \code{read.wc} reads and formats tag data output from Wildlife Computers Data Portal
#' 
#' @param ptt is individual ID number
#' @param wd is the directory where your data lives
#' @param tag is POSIXct object of the tagging date
#' @param pop is POSIXct object of the pop-up date
#' @param type is character indicating which type of data to read. Choices are 'sst', 'pdt', 'light' corresponding to those data files output from WC Data Portal
#'   
#' @return a list containing: the data read as a data.frame and a date vector of unique dates in that data
#' 
#' @examples
#' # READ IN TAG DATA
#'ptt <- '141256'
#'
#' # TAG/POPUP DATES AND LOCATIONS (dd, mm, YYYY, lat, lon)
#' iniloc <- data.frame(matrix(c(13, 10, 2015, 41.575, -69.423, 
#'                              24, 2, 2016, 26.6798, -69.0147), nrow = 2, ncol = 5, byrow = T))
#' colnames(iniloc) = list('day','month','year','lat','lon')
#' tag <- as.POSIXct(paste(iniloc[1,1], '/', iniloc[1,2], '/', iniloc[1,3], sep=''), format = '%d/%m/%Y')
#' pop <- as.POSIXct(paste(iniloc[2,1], '/', iniloc[2,2], '/', iniloc[2,3], sep=''), format = '%d/%m/%Y')
#' myDir <- '~/Documents/WHOI/RCode/hmmwoa/inst/extdata/' # WHERE YOUR DATA LIVES, THIS IS THE EXAMPLE DATA
#' tag.sst <- read.wc(ptt, wd = myDir, type = 'sst', tag=tag, pop=pop)
#' str(tag.sst)



read.wc <- function(ptt, wd = getwd(), tag, pop, type = 'sst'){
  
  if(substr(wd, nchar(wd), nchar(wd)) == '/'){
    
  } else{
    wd <- paste(wd, '/', sep='')
  }

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
    dts <- as.POSIXct(data$Day, format = findDateFormat(data$Day), tz = 'UTC')
    d1 <- as.POSIXct('1900-01-02') - as.POSIXct('1900-01-01')
    didx <- (dts > (tag + d1)) & (dts < (pop - d1))
    data <- data[didx,]
    udates <- unique(as.Date(dts))
    
  }
  
  return(list(data = data, udates = udates))
           
}

