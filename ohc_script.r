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
# how to paramaterize sdx? seems like it should scale relative to
# magnitude of ohc?
# maybe best way is, after testing double tag data, to use known locations
# and ohc at those positions compared to that calculated from the tags
# pdt data. that should be a great way to get sd estimates for ohc based on real data
likelihood = calc.ohc(pdt,ohc.dir=ohc.dir,ptt=121325,sdx=.5,plot=T)

pdf('lyd likelihood.pdf')
image.plot(lon,lat,likelihood[,,15])
dev.off()

########

install.packages('devtools')
library(devtools)
install_github('camrinbraun/hmmwoa')

