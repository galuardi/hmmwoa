# run it

# load and format tag data
# read in data
setwd('~/Documents/WHOI/RData/WhiteSharks/2013/121325/')
data = read.table('121325-PDTs.csv',sep=',',header=T,blank.lines.skip=F)

pdt = extract.pdt(data)

# download daily hycom products
#time <- c(as.Date(min(pdt$Date)),as.Date(max(pdt$Date)))
lon <- c(-80,-30)
lat <- c(15,40)

ohc.dir <- '~/Documents/WHOI/RData/HYCOM/Lydia/'

for(i in 1:length(pdt$Date)){
  time <- as.Date(pdt$Date[i])
  get.hycom(lon,lat,time,filename=paste(ptt,'_',time,'.nc',sep=''),download.file=TRUE,dir=ohc.dir) # filenames based on dates from above

}


# calc.ohc
likelihood = calc.ohc(pdt,ohc.dir=ohc.dir)


