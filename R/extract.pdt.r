extract.pdt = function(data){
  
  # extract PDT data from WC tag output
  # need to convert to long format for use in R
  
  #' @param: data is data frame read from .csv output of Wildlife Computers
  #' DAP processor. File ends in "-PDTs.csv"
  #' 
  #' @return: pdt is formatted data frame of pdt data

  # eliminate any oxygen data
  if(any(grep('X.Ox', colnames(data)))){
    dropidx <- c(grep('Ox', names(data)), grep('Disc', names(data)))
    data <- data[,-dropidx]
  }
  
  # convert to long format
  vars = names(data[,c(which(names(data) == 'Depth1'):length(names(data)))])
  pdt <- reshape(data, ids = data$Date, direction = 'long',
                 varying = vars, times = vars, sep='', timevar = 'BinNum')
  keepNames = c('Ptt', 'Date', 'NumBins', 'BinNum', 'Depth', 'MinTemp', 'MaxTemp')
  pdt <- pdt[,c(keepNames)]
  row.names(pdt) <- NULL
  
  # date conversion then sort
  pdt$Date <- as.POSIXct(pdt$Date, format = findDateFormat(pdt$Date))
  pdt <- pdt[order(pdt$Date, pdt$Depth),]               
  #pdt <- pdt[which(!is.na(pdt$Depth)),]
  pdt <- pdt[!is.na(pdt$Depth),]
  
  udates <- unique(format(pdt$Date, '%Y-%m-%d'))
  dates.tr <- format(pdt$Date, '%Y-%m-%d')
  dateidx <- match(dates.tr, udates)
  
  ## THIS SECTION NEEDS WORK. CURRENTLY USING MULTIPLE POINTS PER DAY
  ## TO CONSTRUCT A SINGLE PDT PROFILE. MOSTLY WORKS WELL BUT IN SOME
  ## CASES THE ANIMAL USES VERY DISTINCT WATER MASSES CAUSING THIS
  ## AVERAGING TECHNIQUE TO DO WEIRD THINGS AND CONSTRUCT UNREALISTIC
  ## PROFILES
  
  for(i in 1:max(dateidx)){
    pdt.i <- pdt[which(dateidx == i),]

    if(length(unique(pdt.i[,5])) <= 2){ 
    
      pdt.t <- pdt.i
      
      } else{
        if(length(which(pdt.i$BinNum == 1)) > 1){
          pdt.i <- pdt.i[order(pdt.i$Depth),]
          z <- unique(pdt.i$Depth)#; z <- sort(z)
          z[z < 0] = 0; z <- unique(z)
          minT <- approx(pdt.i$Depth, pdt.i$MinTemp, xout = z)$y
          maxT <- approx(pdt.i$Depth, pdt.i$MaxTemp, xout = z)$y
          pdt.t <- pdt.i[1:length(z),]
          pdt.t[,c(5:7)] <- cbind(z,minT,maxT)
          pdt.t[,4] <- seq(1, length(z), by = 1)
          pdt.t[,3] <- length(z)
          pdt.t[,2] <- paste(format(pdt.t$Date, '%Y-%m-%d'),' 00:00:00', sep = '')
        } else{
          pdt.t <- pdt.i
        }
      
      }
    
    if(i == 1){
      pdtNew <- pdt.t
    } else{
      pdtNew <- rbind(pdtNew, pdt.t)
    }
    
  }
  # write out / return
  return(pdtNew)
}
