# run it

# load and format tag data
# read in data
setwd('~/Documents/WHOI/RData/WhiteSharks/2013/121325/')
ptt = 121325
data = read.table('121325-PDTs.csv',sep=',',header=T,blank.lines.skip=F)

pdt = extract.pdt(data)

# download daily hycom products
#time <- c(as.Date(min(pdt$Date)),as.Date(max(pdt$Date)))
lon = c(-90, -40)
lat = c(10, 50)

ohc.dir <- '~/Documents/WHOI/RData/HYCOM/Lydia/'

udates <- unique(as.Date(pdt$Date))

for(i in 129:length(udates)){
  time <- as.Date(udates[i])
  repeat{
    get.hycom(lon,lat,time,filename=paste(ptt,'_',time,'.nc',sep=''),download.file=TRUE,dir=ohc.dir) # filenames based on dates from above
    #err <- try(open.ncdf(paste(ohc.dir,ptt,'_',time,'.nc',sep='')),silent=T)
    tryCatch({
      err <- try(open.ncdf(paste(ohc.dir,ptt,'_',time,'.nc',sep='')),silent=T)
    }, error=function(e){print(paste('ERROR: Download of data at ',time,' failed. Trying call to server again.',sep=''))})
    if(class(err) != 'try-error') break
  }
}


# calc.ohc
likelihood = calc.ohc(pdt,ohc.dir=ohc.dir)

tryCatch({
  print(i)
  if (i==7) stop("Urgh, the iphone is in the blender !")
}, error=function(e){})

tryCatch({
  err <- try(open.ncdf(paste(ohc.dir,ptt,'_',time,'.nc',sep='')),silent=T)
}, error=function(e){})

