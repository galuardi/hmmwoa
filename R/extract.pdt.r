extract.pdt = function(pdt){
  
  # extract PDT data from WC tag output
  # need to convert to long format for use in R
  
  #' @param: 
  #' 

  # format accordingly
  
  data = data[,c(1:which(colnames(data)=='MaxTemp8'))]
  vars = names(data[,c(which(names(data)=='Depth1'):length(names(data)))])
  pdt <- reshape(data, ids = data$Date, direction='long',
                 varying = vars, times = vars, sep='', timevar = 'BinNum')
  keepNames = c('Ptt','Date','NumBins','BinNum','Depth','MinTemp','MaxTemp')
  pdt <- pdt[,c(keepNames)]
  row.names(pdt) <- NULL
  
  # need to figure out date conversion then sort
  pdt$Date <- as.Date(pdt$Date)
  pdt <- pdt[order(pdt$Date,pdt$Depth),]               
  
  
  
  # write out / return
    
}
