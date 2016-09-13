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


read.wc <- function(ptt, iniloc){
  
  # READ IN PDT DATA FROM WC FILES
  pdt <- read.table(paste(ptt,'-PDTs.csv', sep=''), sep=',',header=T,blank.lines.skip=F, skip = 0)
  pdt <- extract.pdt(pdt)
  tag <- as.POSIXct(paste(iniloc[1,1], '/', iniloc[1,2], '/', iniloc[1,3], sep=''), format = '%d/%m/%Y')
  pop <- as.POSIXct(paste(iniloc[2,1], '/', iniloc[2,2], '/', iniloc[2,3], sep=''), format = '%d/%m/%Y')
  dts <- as.POSIXct(pdt$Date, format = findDateFormat(pdt$Date))
  d1 <- as.POSIXct('1900-01-02') - as.POSIXct('1900-01-01')
  didx <- dts >= (tag + d1) & dts <= (pop - d1)
  pdt <- pdt[didx,]
  pdt.udates <- unique(as.Date(pdt$Date))
  
  #==============  
  
  # VECTOR OF DATES FROM DATA. THIS IS USED IN MANY FUNCTIONS 
  dateVec <- as.Date(seq(tag, pop, by = 'day'))

  #==============    
  
  # READ IN LIGHT DATA FROM WC FILES
  locs <- read.table(paste(ptt, '-Locations-GPE2.csv', sep=''), sep=',', header = T, blank.lines.skip = F)
  dts <- format(as.POSIXct(locs$Date, format = findDateFormat(locs$Date)), '%Y-%m-%d')
  didx <- dts > (tag + d1) & dts < (pop - d1)
  locs <- locs[didx,]
  
  #==============    
  
  # READ IN TAG SST FROM WC FILES
  tag.sst <- read.table(paste(ptt, '-SST.csv', sep=''), sep=',',header=T, blank.lines.skip=F)
  dts <- as.POSIXct(tag.sst$Date, format = findDateFormat(tag.sst$Date))
  didx <- dts >= (tag + d1) & dts <= (pop - d1)
  tag.sst <- tag.sst[didx,]
  if (length(tag.sst[,1]) <= 1){
    stop('Something wrong with reading and formatting of tags SST data. Check date format.')
  }
  dts <- as.POSIXct(tag.sst$Date, format = findDateFormat(tag.sst$Date))
  sst.udates <- unique(as.Date(dts))
  
  # need to return pdt, locs, dateVec, sst.udates, pdt.udates
  return(list(dateVec = dateVec, pdt = pdt, pdt.udates = pdt.udates,
              locs = locs, tag.sst = tag.sst, sst.udates = sst.udates))
  
}

