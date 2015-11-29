# run it

# load and format tag data
# read in data
setwd('~/Documents/WHOI/RData/Swords/2013/106795/')
ptt = 106795
data = read.table('106795-PDTs.csv',sep=',',header=T,blank.lines.skip=F, skip = 2)

pdt = extract.pdt(data)

# download daily hycom products
#time <- c(as.Date(min(pdt$Date)),as.Date(max(pdt$Date)))
lon = c(-70, -15)
lat = c(20, 60)

ohc.dir <- paste('~/Documents/WHOI/RData/HYCOM/', ptt, '/',sep = '')

udates <- unique(as.Date(pdt$Date))

for(i in 1:5){#length(udates)){
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
likelihood <- calc.ohc(pdt, ohc.dir = ohc.dir, ptt = 106795, sdx = 23.14)

plot.ohc(lik = likelihood, filename='lyd lik.pdf', write.dir = getwd())


# getting near end of this working but need to confirm sdx value. maybe
# its actually sd of hycom ohc for that day?



pdf('lyd likelihood.pdf')
image.plot(lon,lat,likelihood[,,15])
dev.off()

########

install.packages('devtools')
library(devtools)
install_github('camrinbraun/hmmwoa')

spot <- read.table('lydia_spot.csv',sep=',',header=T)


pdf('try dat.pdf')
image.plot(lon,lat,lik)
if(length(which(sdays==time))>0){
  points(spot[which(sdays==time),c(8,7)])
}
dev.off()

pdf('try spots.pdf')
for(i in 1:length(udates)){
  time <- as.Date(udates[i])
  image.plot(lon,lat,likelihood[,,i],zlim=c(0,1))
  if(length(which(sdays==time))>0){
    points(spot[which(sdays==time),c(8,7)],col='white')
  }
  title(paste(time))
}
dev.off()


which.min(abs(lon-(-60)))
which.min(abs(lat-(20.5)))

res <- matrix(ncol=2,nrow=length(udates))

for(i in 1:length(udates)){
  # time and pdt
  time <- as.Date(udates[i])
  pdt.i <- pdt[which(as.Date(pdt$Date)==time),]
  isotherm = min(pdt.i$MidTemp,na.rm=T)
  
  # get day's spot locations and find mean position
  spot.i <- spot[which(sdays==time),c(8,7)]
  mlon <- mean(spot.i[,1])
  mlat <- mean(spot.i[,2])
  
  # load day's hycom
  nc = open.ncdf(paste(ohc.dir,ptt,'_',time,'.nc',sep=''))
  dat = get.var.ncdf(nc, 'temperature')
  dat[dat<isotherm] = NA
  dat = dat - isotherm
  ohc = cp*rho*apply(dat, 1:2, sum, na.rm=T)/10000 
  
  # sample hycom at that position
  idx <- c(which.min(abs(lon-(mlon))), which.min(abs(lat-(mlat))))
  hyc.i <- ohc[idx[1],idx[2]]
  
  # calc ohc for tag and store
  tag = approx(pdt.i$Depth,pdt.i$MidTemp,xout=depth)
  tag = tag$y - isotherm
  tag.ohc = cp*rho*sum(tag,na.rm=T)/10000
  
  # compare hycom to that day's tag-based ohc
  res[i,] <- c(hyc.i,tag.ohc)
  print(i)
}

subres <- res[which(res[,1]!='NA' & res[,1]>0),]

pdf('try.pdf',width=11,height=8)
par(mfrow=c(1,2))
image.plot(lon,lat,lik)
points(spot[which(sdays==time),c(8,7)],col='white')
plot(dat.i[1:19],depth[1:19],ylim=c(1000,0),type='l',xlab='temp',ylab='depth',lwd=2,xlim=c(8,31))
lines(pdt.i$MidTemp,pdt.i$Depth,col='blue',lwd=2)
dev.off()
#layout(matrix(c(1,2,1,3), 2, 2, byrow = TRUE))



#############################


if(plot==T){
  par(mfrow=c(1,2))
  image.plot(lon,lat,lik)
  if(length(which(sdays==time))>0){
    points(spot[which(sdays==time),c(8,7)],col='white')
  }
  title(paste(time))
  plot(dat.i[1:19],depth[1:19],ylim=c(1000,0),type='l',xlab='temp',ylab='depth',lwd=2,xlim=c(8,31))
  lines(pdt.i$MidTemp,pdt.i$Depth,col='blue',lwd=2)
  
}

if(plot == T)
  
  
  library(raster)
f <- dat
r <- raster(t(f[,,1]),xmn=-90,xmx=-40,ymx=55,ymn=10)
r <- flip(r,direction='y')
plot(r)
plot(countriesLow, add = TRUE)

## check we get sensible orientation and scaling and offset and
## interpretation etc. etc.
#library(rworldmap)
#data(countriesLow)

projection(r) <- "+proj=longlat +datum=WGS84 +ellps=WGS84"
newproj <- CRS("+init=epsg:4705")
pr1 <- projectRaster(r, crs=newproj)
plot(pr1)
plot(countriesLow, add = TRUE)

r <- raster(nrow=18, ncol=36)
m <- matrix(1:ncell(r), nrow=18)
r[] <- as.vector(t(m))
extent(r) <- extent(0, 360, -90, 90)
rr <- rotate(r)


getOHCsd <- function(pdt, isotherm){
  pdt$MidTemp <- (pdt$MaxTemp + pdt$MinTemp) / 2
  mtidx = pdt$MinTemp <= isotherm
  sum(pdt$MidTemp[mtidx]-pdt$MinTemp[mtidx], na.rm=T)/nrow(pdt)  
  # average of summed differences where minTemp <=isotherm. 
  # The idea is that the midtemps in the index would not be as representative at those depths
}


