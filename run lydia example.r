# RUN LYDIA EXAMPLE

# calculate light-based likelihood
setwd('~/Documents/WHOI/RData/WhiteSharks/2013/121325/')
ptt <- 121325

iniloc <- data.frame(matrix(c(3, 3, 2013, 30.3917, -81.3802, 
                              31, 8, 2013, 30.668, -79.972), nrow = 2, ncol = 5, byrow = T))
colnames(iniloc) = list('day','month','year','lat','lon')

pdt <- read.table(paste(ptt,'-PDTs.csv', sep=''), sep=',',header=T,blank.lines.skip=F, skip = 0)

pdt <- extract.pdt(pdt)
tag <- as.POSIXct(paste(iniloc[1,1], '/', iniloc[1,2], '/', iniloc[1,3], sep=''), format = '%d/%m/%Y')
pop <- as.POSIXct(paste(iniloc[2,1], '/', iniloc[2,2], '/', iniloc[2,3], sep=''), format = '%d/%m/%Y')
dts <- as.POSIXct(pdt$Date, format = findDateFormat(pdt$Date))
didx <- dts >= tag & dts <= pop
pdt <- pdt[didx,]

lon = c(-90, -40)
lat = c(10, 55)

udates <- unique(as.Date(pdt$Date))
dateVec <- as.Date(seq(tag, pop, by = 'day'))

##
# OHC / HYCOM
##

## LET'S IGNORE HYCOM / OHC FOR NOW. CURRENTLY LYDIA'S TIMESPAN ISN'T
## AVAILABLE IN THE UNIFORM SPATIAL PROJECTION. THIS ISN'T A HUGE
## ISSUE AS I'VE MANAGED TO DEAL WITH THAT BUT I'D RATHER GET
## THIS TO YOU NOW AND JUST USE WOA. ONCE THE REST OF THE ROUTINE
## IS WORKING, OHC IS PIECE OF CAKE TO DROP IN.

ohc = FALSE
if (ohc){
  ohc.dir <- paste('~/Documents/WHOI/RData/HYCOM/', ptt, '/',sep = '')
  
  for(i in 1:length(udates)){
    time <- as.Date(udates[i])
    repeat{
      get.hycom(lon,lat,time,filename=paste(ptt,'_-',time,'.nc',sep=''),download.file=TRUE,dir=ohc.dir, vars = 'water_temp') # filenames based on dates from above
      #err <- try(open.ncdf(paste(ohc.dir,ptt,'_',time,'.nc',sep='')),silent=T)
      tryCatch({
        err <- try(open.ncdf(paste(ohc.dir,ptt,'_',time,'.nc',sep='')),silent=T)
      }, error=function(e){print(paste('ERROR: Download of data at ',time,' failed. Trying call to server again.',sep=''))})
      if(class(err) != 'try-error') break
    }
  }
  
  # calc.ohc
  L.ohc <- calc.ohc(pdt, ohc.dir = ohc.dir)
  
  plot.ohc(lik = L.ohc, ohc.dir = ohcdir, pdt = pdt.data, 
           filename = paste(ptt,'_ohclik.pdf', sep = ''), write.dir = getwd())
}

##
# PDT / WOA
##

# set limits of interest
limits = c(lon, lat) # (min lon, max lon, min lat, max lat)

woa.dir = '/Users/Cam/Documents/WHOI/RData/pdtMatch/WOA_25deg/global/'

return.woa = extract.woa(woa.dir, limits, resolution = 'quarter')
dat = return.woa$dat; lon = return.woa$lon; lat = return.woa$lat; depth = return.woa$depth

# eliminate Pacific from woa data
dat = removePacific(dat, lat, lon)

# check woa data
graphics.off()
image.plot(lon,lat,dat[,,1,1])

# perform matching
# 'stack' makes the end of this routine much slower than 'brick' or 'array'
# but is only 10 extra seconds or so
L.pdt <- calc.pdt(pdt, dat, lat, lon, raster = 'stack', dateVec = dateVec)

# try quick plot to check, if raster = 'stack' or 'brick' above
data(countriesLow)
plot(L.pdt[[100]])
plot(countriesLow, add = T)

plot = FALSE
if(plot){
  plot.woa(as.array(L.pdt), return.woa, paste(ptt, '_woalik.pdf', sep=''), pdt = pdt, write.dir = getwd())
}

##
# Light-based Longitude Likelihood (ellipse error is a work in progress)
##

locs <- read.table(paste(ptt, '-Locations.csv', sep=''), sep=',', header = T, blank.lines.skip = F)
dts <- format(as.POSIXct(locs$Date, format = findDateFormat(locs$Date)), '%Y-%m-%d')
didx <- dts > tag & dts < pop
locs <- locs[didx,]

g <- setup.grid(locs, res = 'quarter') # make sure loading function from misc_funs.r
ngrid <- rev(dim(g$lon))
lon <- g$lon[1,]
lat <- g$lat[,1]

L.locs <- calc.locs(locs, iniloc, g, raster = 'stack', dateVec = dateVec)
# try quick plot to check, if raster = 'stack' or 'brick' above
plot(L.locs[[2]])
plot(countriesLow, add = T)

# sync resolutions of pdt to locs to match grid, g
L.pdt <- spatial_sync_raster(L.pdt, L.locs)

plot(L.pdt[[4]])
plot(countriesLow, add = T)

# multiply daily likelihood matrices
T <- dim(L.pdt)[3]
idx.pdt <- vector('logical', length = T)
idx.locs <- vector('logical', length = T)

for(i in 2:(T-1)){
  idx.pdt[i] <- any(as.matrix(L.pdt[[i]]) != 0)
  idx.locs[i] <- any(!is.na(as.matrix(L.locs[[i]])))
}

dateIdx <- 0
for(i in 2:(T-1)){
  idx <- which(c(idx.locs[i], idx.pdt[i]))
  if(sum(idx) == 3){
    r <- L.pdt[[i]] * L.locs[[i]]
  } else if(sum(idx) == 2){
    r <- L.pdt[[i]]
  } else if(sum(idx) == 1){
    r <- L.locs[[i]]
  } else if(sum(idx) == 0){
    dateIdx = c(dateIdx, i)
  }
 # r <- L.pdt[[i]] * L.locs[[i]]
  if(i == 2){
    s <- stack(r)
  } else{
    s <- stack(s, r)
  }
}

# add known tag/pop locations
tagL <- spatial_sync_raster(L.locs[[1]], s)
popL <- spatial_sync_raster(L.locs[[T]], s)
s <- stack(tagL, s, popL)

# cut out days for which no pdt/loc data exists
#dateIdx <- sort(unique(c(which(dateVec %in% as.Date(pdt$Date)), c(1,which(dateVec %in% as.Date(locs$Date)),T))))
s.sub <- subset(s, which(!(seq(1:T) %in% dateIdx)))

T <- dim(s.sub)[3]

# trim out matrices from time array that don't have any likelihood 
# data right now. this will be improved later.
#L.locs <- L.locs.a[,,c(1,which(dateVec %in% as.Date(locs$Date)),T)]

# now need to re-format the array to match dims in sphmm (time, lat, lon)
L <- aperm(as.array(flip(s.sub, direction = 'y')), c(3,2,1))
# check that it worked ok
image.plot(lon, lat, L[3,,])

# try sphmm
## Number of time steps
T <- dim(L)[1]

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
dt <- 1
G1 <- make.kern(D1,g)
K1 <- uniformization(G1,dt)
##
# [1] Error in diag(A) : no method for coercing this S4 class to a vector

G2 <- make.kern(D2,g)
K2 <- uniformization(G2,dt)
P <- matrix(c(p[1],1-p[1],1-p[2],p[2]),2,2,byrow=TRUE)

## Run smoother and filter
f <- hmm.filter(g,L,K1,K2,P)
s <- hmm.smoother(f,K1,K2,P)
sphmm <- calc.track(s,g)
sphmm$date <- lsst$date
sphmm$p.resid <- apply(s,c(1,2),sum)[2,]

# can see known positions from her SPOT tag
tr <- read.table(paste(ptt,'-SPOT.csv', sep=''), sep=',',header=T,blank.lines.skip=F, skip = 0)
tr <- tr[,c(4, 7, 8)]
tr$b <- 1 # set an arbitrary state at this point
colnames(tr) <- list('date','lat','lon','b')
dts <- as.POSIXct(tr$date, format = findDateFormat(tr$date))
didx <- dts >= tag & dts <= pop
tr <- tr[didx,]

plot.results(save.plot = F)


