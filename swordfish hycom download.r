# swordfish env downloads
woa = ncdf::get.var.ncdf(nc, start = c(xmin, ymin, 1, 1), count = c(xlen + 1, ylen + 1, 57, 12))

dts <- c('25-Nov-2010','20-Jul-2011','15-Sep-2011',
         '13-Apr-2012','27-Sep-2013','26-Mar-2014')
dts <- as.Date(dts, format = '%d-%b-%Y')
udates <- c(seq(dts[1], dts[2], by = 'day'),
            seq(dts[3], dts[4], by = 'day'),
            seq(dts[5], dts[6], by = 'day'))
locs <- read.table('~/Documents/WHOI/RData/Swords/sword bounds.csv',sep=',',header=T)
locs[14,5] = locs[14,5] * -1
plot(locs$TagLong, locs$TagLat, pch=16, col='green', ylim = c(0,60), xlim = c(-85,-20))
points(locs$PopLong, locs$PopLat, pch=16, col='red')
plot(countriesLow, add=T)

sp.lim <- list(lonmin = -85, lonmax = -20, latmin = 0, latmax = 60)

ohc.dir <- paste('~/Documents/WHOI/RData/HYCOM/Swords/',sep = '')

for(i in 1:length(udates)){
  time <- as.Date(udates[i])
  repeat{
    get.hycom(sp.lim, time, type='a', filename = paste('Swords_', time, '.nc', sep = ''),
              download.file = TRUE, dir = ohc.dir, vars = 'water_temp') 
    tryCatch({
      err <- try(ncdf::open.ncdf(paste(ohc.dir,ptt,'_',time,'.nc',sep='')),silent=T)
    }, error=function(e){print(paste('ERROR: Download of data at ',time,' failed. Trying call to server again.',sep=''))})
    if(class(err) != 'try-error') break
  }
}



