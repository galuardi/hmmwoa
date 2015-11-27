# calculate light-based likelihood
setwd('~/Documents/WHOI/RData/Swords/2013/106795/')
ptt <- 106795
iniloc <- data.frame(matrix(c(27, 9, 2013, 46.47683333, -45.5640, 
                   2, 11, 2013, 30.92645, -39.6919), nrow = 2, ncol = 5, byrow = T))

pdt <- read.table('106795-PDTs.csv',sep=',',header=T,blank.lines.skip=F, skip = 2)
pdt <- pdt[,c(grep('X.Ox', colnames(pdt)) * -1)]
pdt <- pdt[,c(grep('Disc', colnames(pdt)) * -1)]
data  <- pdt

pdt <- extract.pdt(pdt)
pdt <- pdt[!is.na(pdt$Depth),]
tag <- as.POSIXct(paste(iniloc[1,1], '/', iniloc[1,2], '/', iniloc[1,3], sep=''), format = '%d/%m/%Y')
pop <- as.POSIXct(paste(iniloc[2,1], '/', iniloc[2,2], '/', iniloc[2,3], sep=''), format = '%d/%m/%Y')
dts <- as.POSIXct(pdt$Date, format = findDateFormat(pdt$Date))
didx <- dts >= tag & dts <= pop
pdt <- pdt[didx,]

lon = c(-70, -15)
lat = c(20, 60)

ohc.dir <- paste('~/Documents/WHOI/RData/HYCOM/', ptt, '/',sep = '')

udates <- unique(as.Date(pdt$Date))

for(i in 1:length(udates)){
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

L.ohc <- calc.ohc(pdt, ohc.dir = ohc.dir, ptt = 106795, sdx = 10)

plot.ohc(lik = L.ohc, ohc.dir = ohcdir, pdt = pdt.data, filename = '106795_lik.pdf', write.dir = getwd())


##
# PDT Match
##

# set limits of interest
limits = c(lon,lat) # (min long, max long, min lat, max lat)

nc.dir = '/Users/Cam/Documents/WHOI/RData/pdtMatch/WOA_25deg/global/'

return.woa = extract.woa(nc.dir, limits, resolution = 'quarter')
dat = return.woa$dat; lon = return.woa$lon; lat = return.woa$lat; depth = return.woa$depth

# eliminate Pacific from matching
dat = removePacific(dat, lat, lon)

# perform matching
L.pdt = calc.pdt(pdt, dat, lat, lon)

plot.woa(L.pdt, return.woa, '106795_woa.pdf', pdt = pdt, write.dir = getwd())


##
# Light-based Longitude Likelihood
##

locs <- read.table('106795-Locations.csv', sep=',', header = T, blank.lines.skip = F)
#light <- light[light$Type == 'GPE',]

ngrid <- c(limits[2] - limits[1], limits[4] - limits[3])
g <- setup.grid(locs)
lon <- g$lon[1,]
lat <- g$lat[,1]

# tag and pop locations
#iniloc = data.frame(matrix(c(-41.38018333,30.39173333,-52.24500,36.688),
#                           nrow=2,ncol=2,byrow=T))

colnames(iniloc) = list('day','month','year','lat','lon')

L.locs <- lik.locs(locs, iniloc, g)


##
# plot loess
##

#data <- read.table(file='106795-PDTs.csv',sep=',',header=T,blank.lines.skip=F, skip = 2)

depths = as.vector(as.matrix(cbind(data[,c(seq(15,ncol(data)-2,by=3))])))
mintemps = as.vector(as.matrix(cbind(data[,c(seq(16,ncol(data)-1,by=3))])))
maxtemps = as.vector(as.matrix(cbind(data[,c(seq(17,ncol(data),by=3))])))
midtemps = (maxtemps + mintemps) / 2

ddates = as.POSIXct(strptime(as.character(data$Date),format = findDateFormat(data$Date))) #reads dates as dates
year = as.numeric(format(ddates, '%Y')) #extracts year
tagyear <- year[1]
DOY = round(julian(ddates,origin=as.Date(paste(year[1],'-01-01',sep=''))),digits=0) #calculate DOY
jday = as.numeric(DOY+((year[1]-tagyear)*365)) #finish DOY calculation

#creates matrix of days
jd_for_interp <- as.vector(as.matrix(cbind(rep(jday, (ncol(data) - 14 / 3)))))

# clean data using NA values in depths
ii=which(!is.na(depths));
depths=depths[ii];midtemps=midtemps[ii];jd_for_interp=jd_for_interp[ii];

# sets depth constraint
z = 1:max(depths, na.rm = T)

# run the LOESS interp
results=grid2dloess(data=midtemps,xgrid=jd_for_interp,ygrid=depths,
                    span_x=5,span_y=150,xgrid_est=jday,ygrid_est=z)

# then plot
plot.days=seq(0,length(jday)-1,by=1) #x axis setup
plot.loess.pdt(plot.days,z,results$sm_data,year,filename='106795_pdt.pdf',DOY.first=min(jday),zissou=FALSE,tempRange=c(4,26))


#######
## END
#######


setwd('~/Documents/WHOI/RData/Swords/2011/110496/')


lsst = read.table('110496-Locations.csv',sep=',',header=T)

source('/Users/Cam/Documents/WHOI/RCode/sphmm/get.sst.r')
sst.file <- "/Users/Cam/Documents/WHOI/RData/sst.wkmean.1990-present.nc"
land.file <- '/Users/Cam/Documents/WHOI/RData/lsmask.nc'
strtDate <- as.Date(lsst$date[1])
endDate <- as.Date(lsst$date[length(lsst[,1])])
sstResults <- get.sst(date1=strtDate,date2=endDate,lon=c(270,310),lat=c(45,10),
                      sst.file=sst.file,land.file=land.file,do.plot=F)
#date1=strtDate;date2=endDate;lon=c(260,350);lat=c(60,-40);
#sst.file=sst.file;land.file=land.file

# read boundary polygon shapefile
sstBound <- readOGR('/Users/Cam/Documents/WHOI/RData/Maps/LydiaBound.shp','LydiaBound')
shp <- fortify(sstBound)

xpos <- shp$long; ypos <- shp$lat
tpos<-c('2013-03-04','2013-03-15')
sstx <- xtractogon(xpos,ypos,tpos,'erdVH2chlamday')
#sstx <- xtractogon(mbnms$Longitude,mbnms$Latitude,as.Date(tpos),'erdAAssta8day')


#for(t in 1:dim(sstResults$sst)[3]){
#  sstResults$sst[,,t] = t(sstResults$sst[,,t])
#}

# remove Pacific from sstResults
# currently sstResults formatted strangely as seems lats are stored as backwards,
# such that southern data points are at lon=65north but axis is then just flipped

# get idx of relevant sst data
sst.idx = rep(0,length(lsst$date))
for(t in 1:length(lsst$date)){
  sst.idx[t] = which.min(abs(as.Date(lsst$date[t],format='%Y-%m-%d %H:%M:%S') - sstResults$dates))
}



run.sphmm <- function(lsst,iniloc,ngrid,par0=c(8.908,10.27,1.152,0.0472,0.707,0.866),do.fit=FALSE,do.plot=TRUE,show.movie=FALSE){
}

#dt <- date2time(lsst$date[2:length(lsst$date)])-date2time(lsst$date[1:(length(lsst$date)-1)])
dt <- 1

## Setup discrete spatial grid w/desired resolution
#ngrid <- c(50,49)
ngrid <- dim(sstResults$sst)[1:2]
# ngrid format is number of cells in c(x,y)

g <- setup.grid(lsst,ngrid)
row <- dim(g$lon)[1]
col <- dim(g$lon)[2]

## Compute likelihood for observations
L <- data.lik(lsst,iniloc,g)


for(i in 1:dim(L)[1]){
  image.plot(g$lon,g$lat,L[i,,])
  invisible(readline(prompt="Press [enter] to plot next time t."))
}



## Calculate the "data" likelihood, i.e. the likelihood surface for each observation time T

data.lik <- function(lsst,iniloc,g){
  ## Calculate the "data" likelihood, i.e. the likelihood field for each observation
  
  T <- length(lsst$lon)
  row <- dim(g$lon)[1]
  col <- dim(g$lon)[2]
  
  L <- array(0,dim=c(col,row,T))
  
  # Initial location is known
  ilo <- which.min(abs(g$lon[1,]-iniloc$lon[1]))
  ila <- which.min(abs(g$lat[,1]-iniloc$lat[1]))
  L[1,ila,ilo] <- 1
  
  # Calculate data likelihood
  # SD for light-based longitude from Musyl et al. (2001)
  sl.sd <- 35/111 # Converting from kms to degrees
  # SD for SST from Pedersen et al. (2011) Oikos
  sst.sd <- 0.71
  
  for(t in 2:(T-1)){
    time <- date2time(lsst$date[t])
    # find time point in sstResults that correspond to time t
    sst <- t(sstResults$sst[,,sst.idx[t]])
    #sst <- apply(sst,2,rev)
    #image.plot(g$lat[,1],g$lon[1,],sst)
    #sstg <- sstdb(time,g$lon,g$lat)
    Lsst <- dnorm(sst,lsst$sst[t],sst.sd)  # SST data
    Llon <- dnorm(g$lon,lsst$lon[t],sl.sd) # Longitude data
    L[t,,] <- Lsst*Llon
  }
  
  # End location is known
  elo <- which.min(abs(g$lon[1,]-iniloc$lon[2]))
  ela <- which.min(abs(g$lat[,1]-iniloc$lat[2]))
  L[T,ela,elo] <- 1
  
  L
}


## Number of time steps
T <- length(lsst$lon)

## Fixed parameter values
par0=c(8.908,10.27,1.152,0.0472,0.707,0.866)
D1 <- par0[1:2]
D2 <- par0[3:4]
p <- par0[5:6]

if(do.fit){
  guess <- c(log(10),log(10),log(0.5),log(0.5),log(0.95/0.05),log(0.95/0.05))
  fit <- nlm(neg.log.lik.fun,guess,g,L,dt)
  D1 <- exp(fit$estimate[1:2])
  D2 <- exp(fit$estimate[3:4])
  p <- 1/(1+exp(-fit$estimate[5:6]))
}

## Setup transition matrices
G1 <- make.kern(D1,g)
K1 <- uniformization(G1,dt)
G2 <- make.kern(D2,g)
K2 <- uniformization(G2,dt)
P <- matrix(c(p[1],1-p[1],1-p[2],p[2]),2,2,byrow=TRUE)

## Run smoother and filter
f <- hmm.filter(g,L,K1,K2,P)
s <- hmm.smoother(f,K1,K2,P)
sphmm <- calc.track(s,g)
sphmm$date <- lsst$date
sphmm$p.resid <- apply(s,c(1,2),sum)[2,]


#plot.results <- function(save.plot=TRUE,show.movie=FALSE){
#require(raster) # if want to use extent()
## Show movement as animation
#if(show.movie) show.movie(s)
## Plot behaviour
pp <- apply(s,c(1,2),sum)[2,]
ind <- seq(1,length(tr$b),length.out=length(pp))
graphics.off();
if(save.plot) pdf('../spathmm/sphmm.pdf',width=7,height=8)
par(mfrow=c(2,1))
tr$date<-as.POSIXct(strptime(tr$date,"%Y-%m-%d %H:%M:%S",tz="GMT"))
plot(I(b-1)~date,tr,col='grey',type='h',ylab='Probability of resident',main='Estimated behaviour',xlab='Date')
lines(pp~tr$date[ind],lwd=2,col='red')
xl <- c(min(tr$lon)-5,max(tr$lon)+5)
yl <- c(min(tr$lat)-15,max(tr$lat)+5)
plot(tr$lon,tr$lat,type='n',main='Estimated movements',ylab='Latitude',xlab='Longitude',xlim=xl,ylim=yl)
rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], col = "steelblue1")
lines(tr$lon,tr$lat,col='white')
lines(sphmm$meanlon,sphmm$meanlat,col='black')
points(tr$lon[1],tr$lat[1],bg='green',pch=21)
TT <- length(tr$lon)
points(tr$lon[TT],tr$lat[TT],bg='red',pch=21)
if(save.plot) dev.off()

sst <- numeric(T)
## Simple diagnostics plot ###
## Resample SST
for(t in 1:T){
  time <- date2time(lsst$date[t])
  sst[t] <- sstdb(time,sphmm$meanlon[t],sphmm$meanlat[t])
}

if(save.plot) pdf('../plot/sphmmDiagn.pdf',width=7,height=7)
if(!save.plot) dev.new()
par(mfrow=c(2,2))
ssterr <- sst-lsst$sst;
sdsst <- sqrt(var(ssterr))
ssterr <- ssterr/sdsst
lonerr <- sphmm$meanlon-lsst$lon;
sdlon <- sqrt(var(lonerr))
lonerr <- lonerr/sdlon
plot(tr$date[ind],ssterr,xlab='Date',ylab='SST residual',main='Residuals through time',ylim=c(-3,3),pch=".",cex=3)
abline(h=c(-2,0,2),col='grey',lwd=2,lty=2)
qqnorm(ssterr,main='QQ-plot of residuals',pch=".",cex=2)
abline(a=0,b=1,col='grey',lwd=2)
plot(tr$date[ind],lonerr,xlab='Date',ylab='Longitude residual',ylim=c(-3,3),pch=".",cex=3)
abline(h=c(-2,0,2),col='grey',lwd=2,lty=2)
qqnorm(lonerr,main='',pch=".",cex=2)
abline(a=0,b=1,col='grey',lwd=2)
#  if(save.plot) dev.off()
#}


# ==========================
# END

# some old gpe formatting
# read in GPE data with all PSAT combine script
# keep tag/pop dates from this load for below
# read in all individual csv files and bind
setwd('/Users/Cam/Documents/WHOI/RData/animove/writedata/')
indivList = list.files(getwd())
for (ii in 1:length(indivList)){
  indiv = read.table(file=indivList[ii],sep='',header=T)
  if(ii==1){ combine = indiv } else{
    combine = rbind(combine, indiv)
  }
}
combine <- combine[which(combine$species=='BET'),]
by(combine[,1], combine$ptt, length)
gpe <- combine[which(combine$ptt==133028),]


#gpe <- read.table(file='~/Documents/WHOI/RData/animove/writedata/')
str(gpe)
gpe$dates = paste(as.numeric(gpe$year),'-',as.numeric(gpe$month),'-',
                  as.numeric(gpe$day),' ',as.numeric(gpe$hour),':',
                  as.numeric(gpe$minute),':',00,sep='')
gpe$dates = as.POSIXct(gpe$dates,format='%Y-%m-%d %H:%M:%S',tz='GMT')

gpe = data.frame(date=as.character(gpe$dates),sst=as.numeric(gpe$sst),lon=as.numeric(gpe$lon))

# trim NAs, usually from tag/pop which will be fixed in sphmm anyway
gpe <- gpe[which(!is.na(gpe$sst)),]

# now do sphmm

# tag/pop set as iniloc
iniloc = data.frame(matrix(c(-28.77405,38.62618,-53.61800,43.23000),
                           nrow=2,ncol=2,byrow=T))
colnames(iniloc) = list('lon','lat')
lsst <- gpe

source('/Users/Cam/Documents/WHOI/RCode/sphmm/get.sst.r')
sst.file <- "/Users/Cam/Documents/WHOI/RData/sst.wkmean.1990-present.nc"
land.file <- '/Users/Cam/Documents/WHOI/RData/lsmask.nc'
strtDate <- as.Date(lsst$date[1])
endDate <- as.Date(lsst$date[length(lsst[,1])])
sstResults <- get.sst(date1=strtDate,date2=endDate,lon=c(270,360),lat=c(60,10),
                      sst.file=sst.file,land.file=land.file,do.plot=F)
date1=strtDate;date2=endDate;lon=c(260,350);lat=c(60,-40);
sst.file=sst.file;land.file=land.file

# EVENTUALLY WANT THIS TO WORK SO WE CAN USE HIGHER RES DATA AND IN REAL-TIME
# read boundary polygon shapefile
#sstBound <- readOGR('/Users/Cam/Documents/WHOI/RData/Maps/LydiaBound.shp','LydiaBound')
#shp <- fortify(sstBound)

#xpos <- shp$long; ypos <- shp$lat
#tpos<-c('2013-03-04','2013-03-15')
#sstx <- xtractogon(xpos,ypos,tpos,'erdVH2chlamday')
#sstx <- xtractogon(mbnms$Longitude,mbnms$Latitude,as.Date(tpos),'erdAAssta8day')


#for(t in 1:dim(sstResults$sst)[3]){
#  sstResults$sst[,,t] = t(sstResults$sst[,,t])
#}

# remove Pacific from sstResults

# get idx of relevant sst data
sst.idx = rep(0,length(lsst$date))
for(t in 1:length(lsst$date)){
  sst.idx[t] = which.min(abs(as.Date(lsst$date[t],format='%Y-%m-%d %H:%M:%S') - sstResults$dates))
}



#run.sphmm <- function(lsst,iniloc,ngrid,par0=c(8.908,10.27,1.152,0.0472,0.707,0.866),do.fit=FALSE,do.plot=TRUE,show.movie=FALSE){
#}

#dt <- date2time(lsst$date[2:length(lsst$date)])-date2time(lsst$date[1:(length(lsst$date)-1)])
dt <- .5

## Setup discrete spatial grid w/desired resolution
#ngrid <- c(50,49)
ngrid <- dim(sstResults$sst)[1:2]
# ngrid format is number of cells in c(x,y)

g <- setup.grid(lsst,ngrid)
row <- dim(g$lon)[1]
col <- dim(g$lon)[2]

## Compute likelihood for observations
L <- data.lik(lsst,iniloc,g)


for(i in 1:dim(L)[1]){
  #image.plot(g$lon[1,],g$lat[,1],L[i,,])
  image(L[i,,])
  invisible(readline(prompt="Press [enter] to plot next time t."))
}

# remove NA values from L



## Calculate the "data" likelihood, i.e. the likelihood surface for each observation time T

data.lik <- function(lsst,iniloc,g){
  ## Calculate the "data" likelihood, i.e. the likelihood field for each observation
  
  T <- length(lsst$lon)
  row <- dim(g$lon)[1]
  col <- dim(g$lon)[2]
  
  L <- array(0,dim=c(T,row,col))
  
  # Initial location is known
  ilo <- which.min(abs(g$lon[1,]-iniloc$lon[1]))
  ila <- which.min(abs(rev(g$lat[,1])-iniloc$lat[1]))
  L[1,ila,ilo] <- 1
  
  # Calculate data likelihood
  # SD for light-based longitude from Musyl et al. (2001)
  sl.sd <- 35/111 # Converting from kms to degrees
  # SD for SST from Pedersen et al. (2011) Oikos
  sst.sd <- 0.71
  
  for(t in 2:(T-1)){
    time <- date2time(lsst$date[t])
    # find time point in sstResults that correspond to time t
    sst <- t(sstResults$sst[,,sst.idx[t]])
    #sst <- apply(sst,2,rev)
    #image.plot(g$lat[,1],g$lon[1,],sst)
    #sstg <- sstdb(time,g$lon,g$lat)
    Lsst <- dnorm(sst,lsst$sst[t],sst.sd)  # SST data
    Llon <- dnorm(g$lon,lsst$lon[t],sl.sd) # Longitude data
    L[t,,] <- Lsst*Llon
  }
  
  # End location is known
  elo <- which.min(abs(g$lon[1,]-iniloc$lon[2]))
  ela <- which.min(abs(rev(g$lat[,1])-iniloc$lat[2]))
  L[T,ela,elo] <- 1
  
  L
}


## Number of time steps
T <- length(lsst$lon)

## Fixed parameter values
par0=c(8.908,10.27,1.152,0.0472,0.707,0.866)
D1 <- par0[1:2]
D2 <- par0[3:4]
p <- par0[5:6]

if(do.fit){
  guess <- c(log(10),log(10),log(0.5),log(0.5),log(0.95/0.05),log(0.95/0.05))
  fit <- nlm(neg.log.lik.fun,guess,g,L,dt)
  D1 <- exp(fit$estimate[1:2])
  D2 <- exp(fit$estimate[3:4])
  p <- 1/(1+exp(-fit$estimate[5:6]))
}

## Setup transition matrices
G1 <- make.kern(D1,g)
K1 <- uniformization(G1,dt)
G2 <- make.kern(D2,g)
K2 <- uniformization(G2,dt)
P <- matrix(c(p[1],1-p[1],1-p[2],p[2]),2,2,byrow=TRUE)

## Run smoother and filter
f <- hmm.filter(g,L,K1,K2,P)
s <- hmm.smoother(f,K1,K2,P)
sphmm <- calc.track(s,g)
sphmm$date <- lsst$date
sphmm$p.resid <- apply(s,c(1,2),sum)[2,]

sphmm$pxdate <- as.POSIXct(sphmm$date,format=findDateFormat(sphmm$date))


xlims = c(min(sphmm$meanlon)-5,max(sphmm$meanlon)+5)
ylims = c(min(sphmm$meanlat)-5,max(sphmm$meanlat)+5)
countries <- readOGR('~/Documents/WHOI/RData/Countries/countries.shp','countries.shp')
plot('worldHires',xlim=xlims,ylim=ylims)
