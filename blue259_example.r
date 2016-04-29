# RUN BLUE 259 VIA HMMWOA
library(hmmwoa)

# SETWD
setwd('~/Documents/WHOI/Data/Blues/2015/141259/') 

#----------------------------------------------------------------------------------#
# ADD MAP DATA
library(rworldmap)
data("countriesLow")

#----------------------------------------------------------------------------------#
# READ IN TAG DATA
ptt <- 141259

# TAGGING LOCATION
iniloc <- data.frame(matrix(c(13, 10, 2015, 41.3, -69.27, 
                              10, 4, 2016, 40.251, -36.061), nrow = 2, ncol = 5, byrow = T))
colnames(iniloc) = list('day','month','year','lat','lon')

# READ IN PDT DATA FROM WC FILES
pdt <- read.table(paste(ptt,'-PDTs.csv', sep=''), sep=',',header=T,blank.lines.skip=F, skip = 0)
pdt <- extract.pdt(pdt)
tag <- as.POSIXct(paste(iniloc[1,1], '/', iniloc[1,2], '/', iniloc[1,3], sep=''), format = '%d/%m/%Y')
pop <- as.POSIXct(paste(iniloc[2,1], '/', iniloc[2,2], '/', iniloc[2,3], sep=''), format = '%d/%m/%Y')
dts <- as.POSIXct(pdt$Date, format = findDateFormat(pdt$Date))
d1 <- as.POSIXct('1900-01-02') - as.POSIXct('1900-01-01')
didx <- dts >= (tag + d1) & dts <= (pop - d1)
pdt <- pdt[didx,]

# VECTOR OF DATES FROM DATA. THIS IS USED IN MANY FUNCTIONS 
udates <- unique(as.Date(pdt$Date))
dateVec <- as.Date(seq(tag, pop, by = 'day'))

#----------------------------------------------------------------------------------#
# LIGHT LIKELIHOOD
# Light-based Longitude Likelihood
#----------------------------------------------------------------------------------#
# READ IN LIGHT DATA FROM WC FILES
locs <- read.table(paste(ptt, '-Locations-GPE2.csv', sep=''), sep=',', header = T, blank.lines.skip = F)
dts <- format(as.POSIXct(locs$Date, format = findDateFormat(locs$Date)), '%Y-%m-%d')
didx <- dts > (tag + d1) & dts < (pop - d1)
locs <- locs[didx,]

# SPATIAL LIMITS
#sp.lim <- list(lonmin = -82, lonmax = -25, latmin = 15, latmax = 50)

if (exists('sp.lim')){
  locs.grid <- setup.locs.grid(sp.lim)
} else{
  locs.grid <- setup.locs.grid(locs)
  sp.lim <- list(lonmin = min(locs.grid$lon[1,]), lonmax = max(locs.grid$lon[1,]),
                 latmin = min(locs.grid$lat[,1]), latmax = max(locs.grid$lat[,1]))
}

## SKIP LIGHT FOR NOW, CURRENT WC OUTPUT DOESN'T INCLUDE THIS ANYMORE
# GET THE LIKELIHOOD ELLIPSES
t <- Sys.time()
L.locs <- calc.locs(locs, iniloc, locs.grid, dateVec = dateVec, errEll = T)
Sys.time() - t # around 20 seconds with user-defined limits, up to 5 mins or so with locs limits

# something here for quick example plot to check L.locs?

#----------------------------------------------------------------------------------#
# SST LIKELIHOOD
#----------------------------------------------------------------------------------#

# READ IN TAG SST FROM WC FILES
tag.sst <- read.table(paste(ptt, '-SST.csv', sep=''), sep=',',header=T, blank.lines.skip=F)
dts <- as.POSIXct(tag.sst$Date, format = findDateFormat(tag.sst$Date))
didx <- dts >= (tag + d1) & dts <= (pop - d1)
tag.sst <- tag.sst[didx,]
if (length(tag.sst[,1]) <= 1){
  stop('Something wrong with reading and formatting of tags SST data. Check date format.')
}

# IF USING SST
{

dts <- as.POSIXct(tag.sst$Date, format = findDateFormat(tag.sst$Date))
udates <- unique(as.Date(dts))

sst.dir <- paste('~/Documents/WHOI/RData/SST/OI/', ptt, '/',sep = '')

for(i in 1:length(udates)){
  time <- as.Date(udates[i])
  repeat{
    get.oi.sst(sp.lim, time, filename = paste(ptt, '_', time, '.nc', sep = ''), download.file = TRUE, dir = sst.dir) # filenames based on dates from above
    tryCatch({
      err <- try(ncdf::open.ncdf(paste(sst.dir, ptt, '_', time, '.nc', sep = '')), silent = T)
    }, error=function(e){print(paste('ERROR: Download of data at ', time, ' failed. Trying call to server again.', sep = ''))})
    if(class(err) != 'try-error') break
  }
}

t <- Sys.time()
L.sst <- calc.sst(tag.sst, sst.dir = sst.dir, dateVec = dateVec)
Sys.time() - t
# focal for SD calc takes about .5 sec each t step, the dnorm is <.3 sec
# total 2 min for blue shk 259

}

#----------------------------------------------------------------------------------#
# OHC / HYCOM LIKELIHOOD(S)
#----------------------------------------------------------------------------------#

# IF USING OHC HYCOM
{

udates <- unique(as.Date(pdt$Date))
ohc.dir <- paste('~/Documents/WHOI/RData/HYCOM/', ptt, '/',sep = '')

for(i in 1:length(udates)){
  time <- as.Date(udates[i])
  repeat{
    get.hycom(sp.lim, time, type='a', filename = paste(ptt, '_', time, '.nc', sep = ''),
              download.file = TRUE, dir = ohc.dir, vars = 'water_temp') 
    tryCatch({
      err <- try(ncdf::open.ncdf(paste(ohc.dir,ptt,'_',time,'.nc',sep='')),silent=T)
    }, error=function(e){print(paste('ERROR: Download of data at ',time,' failed. Trying call to server again.',sep=''))})
    if(class(err) != 'try-error') break
  }
}

# calc.ohc
t <- Sys.time()
L.ohc <- calc.ohc(pdt, ohc.dir = ohc.dir, dateVec = dateVec, isotherm = '')
Sys.time() - t
# focal takes <8 secs and dnorm 2-7 secs for each t step (day)

}


