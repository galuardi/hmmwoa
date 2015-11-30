extract.pdt = function(data){
  
  # extract PDT data from WC tag output
  # need to convert to long format for use in R
  
  #' @param: data is data frame read from .csv output of Wildlife Computers
  #' DAP processor. File ends in "-PDTs.csv"
  #' 
  #' @return: pdt is formatted data frame of pdt data
  #' 
  #' Written by Camrin Braun (WHOI) for R on 12 Nov 2015
  #' cbraun@whoi.edu

  # format accordingly
  
  #data = data[,c(1:which(colnames(data)=='MaxTemp8'))]
  
  if(any('X.Ox1' %in% names(data))){
    dropidx <- c(grep('Ox', names(data)), grep('Disc', names(data)))
    data <- data[,-dropidx]
  }
  
  vars = names(data[,c(which(names(data)=='Depth1'):length(names(data)))])
  pdt <- reshape(data, ids = data$Date, direction='long',
                 varying = vars, times = vars, sep='', timevar = 'BinNum')
  keepNames = c('Ptt','Date','NumBins','BinNum','Depth','MinTemp','MaxTemp')
  pdt <- pdt[,c(keepNames)]
  row.names(pdt) <- NULL
  
  # need to figure out date conversion then sort
  pdt$Date <- as.POSIXct(pdt$Date, format=findDateFormat(pdt$Date))
  pdt <- pdt[order(pdt$Date,pdt$Depth),]               
  pdt <- pdt[which(!is.na(pdt$Depth)),]
  # write out / return
  return(pdt)
  
}
