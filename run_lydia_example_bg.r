
# RUN LYDIA EXAMPLE

library(fields)
library(raster)
library(imager)
library(ncdf)
library(plyr)
#library(abind) # don't need this with modified funciton
library(reshape2)
library(rworldmap)
library(spatial.tools)

# calculate light-based likelihood
setwd('C:/Users/ben/Google Drive/Camrin-WOA/hmmwoa_files/')
#setwd('~/Documents/WHOI/RData/WhiteSharks/2013/121325/')

source('C:\\Users\\ben\\Documents\\GitHub\\hmmwoa\\R\\findDateFormat.r')
source('C:\\Users\\ben\\Documents\\GitHub\\hmmwoa\\R\\extract.pdt.r')
source('C:\\Users\\ben\\Documents\\GitHub\\hmmwoa\\R\\extract.woa.r')
source('C:\\Users\\ben\\Documents\\GitHub\\hmmwoa\\R\\removePacific.r')
source('C:\\Users\\ben\\Documents\\GitHub\\hmmwoa\\R\\calc.pdt.r')
source('C:\\Users\\ben\\Documents\\GitHub\\hmmwoa\\R\\plot.woa.r')
source('C:\\Users\\ben\\Documents\\GitHub\\hmmwoa\\R\\calc.locs.r')
source('C:\\Users\\ben\\Documents\\GitHub\\hmmwoa\\R\\misc_funs.r')


ptt <- 121325

iniloc <- data.frame(matrix(c(3, 3, 2013, 30.3917, -81.3802, 
                              31, 8, 2013, 30.668, -79.972), nrow = 2, ncol = 5, byrow = T))
colnames(iniloc) = list('day','month','year','lat','lon')

pdt <- read.table(paste(ptt,'-PDTs.csv', sep=''), sep=',',header=T,blank.lines.skip=F, skip = 0)

pdt <- extract.pdt(pdt)
tag <- as.POSIXct(paste(iniloc[1,1], '/', iniloc[1,2], '/', iniloc[1,3], sep=''), format = '%d/%m/%Y')
pop <- as.POSIXct(paste(iniloc[2,1], '/', iniloc[2,2], '/', iniloc[2,3], sep=''), format = '%d/%m/%Y')
dts <- as.POSIXct(pdt$Date, format = findDateFormat(pdt$Date))
didx <- dts >= tag+1 & dts <= pop-1
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

# woa.dir = '/Users/Cam/Documents/WHOI/RData/pdtMatch/WOA_25deg/global/'
woa.dir = "C:/Users/ben/Google Drive/Camrin-WOA/hmmwoa_files/"

return.woa = extract.woa(woa.dir, limits, resolution = 'quarter')
dat = return.woa$dat; 
lon = as.numeric(return.woa$lon); 
lat = as.numeric(return.woa$lat); 
depth = as.numeric(return.woa$depth)

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
plot(L.pdt[[2]])
plot(countriesLow, add = T)

# plot = FALSE
# if(plot){
  # plot.woa(as.array(L.pdt), return.woa, paste(ptt, '_woalik.pdf', sep=''), pdt = pdt, write.dir = getwd())
# }

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

L.locs <- calc.locs(locs, iniloc, g, raster = T, dateVec = dateVec)
# try quick plot to check, if raster = 'stack' or 'brick' above
plot(L.locs[[181]])
plot(countriesLow, add = T)

#plot WOA and light likelihood together
plot(L.locs[[4]]*L.pdt[[4]])
plot(countriesLow, add = T)

# sync resolutions of pdt to locs to match grid, g
L.pdt <- spatial_sync_raster(L.pdt, L.locs)

plot(L.pdt[[4]])
plot(countriesLow, add = T)

# multiply daily likelihood matrices

# check sums of L components

# need an index of where likelihoods are zeros.. for each L component

L.locs = as.array(L.locs)
L.pdt = as.array(L.pdt)
L.locs[is.na(L.locs)] = 0 # turn NA to 0
L.pdt[is.na(L.pdt)] = 0

# you're the king of apply(). it's so handy!
nalocidx = apply(L.locs,3, sum, na.rm=T)==0 # does sum of likelihood surface
# at each time point == 0?
napdtidx = apply(L.pdt,3, sum, na.rm=T)==0

naLidx = nalocidx+napdtidx # where both are zeros. These will be interpolted in the filter
dateIdx = naLidx==0 # may not need this but here for now..

Lmat = L.pdt*0
# where naLidx==0, both likelihoods are zero
#       naLidx==1, one has data
#       naLidx==2, both have data
idx1 = naLidx==1
idx2 = naLidx==2

Lmat[,,idx1] = L.pdt[,,idx1]+L.locs[,,idx1] # when only 1 has data
Lmat[,,idx2] = L.pdt[,,idx2]*L.locs[,,idx2] # when both have data

crs <- "+proj=longlat +datum=WGS84 +ellps=WGS84"
list.pdt <- list(x = lon, y = lat, z = L.pdt)
ex <- extent(list.pdt)

for(i in 1:T){
  L.i <- raster(Lmat[,,i], xmn=ex[1], xmx=ex[2], ymn=ex[3], ymx=ex[4], crs)
  if(i==1) L <- L.i else L <- stack(L, L.i)
}
#Lmatr <- raster(Lmat[,,1], xmn=ex[1],xmx=ex[2],ymn=ex[3],ymx=ex[4])
# now need to re-format the array to match dims in sphmm (time, lat, lon)
L.arr <- aperm(as.array(flip(L, direction = 'y')), c(3,2,1)) # this is transposed and is still time, lon, lat
#L <- aperm(Lmat, c(3,2,1))  # using arrays..
# L <- aperm(as.array(flip(s.sub, direction = 'y')), c(3,1,2))  # this looks better.. or not!
# check that it worked ok
ex <- extent(L)
lon <- seq(ex[1], ex[2], length=dim(Lras)[2])
lat <- seq(ex[3], ex[4], length=dim(Lras)[1])
image.plot(lon, lat, L.arr[1,,])

L <- L.arr

## ******
## now insert portions of imager script
library(imager)

par0=c(8.908,10.27,3,1,0.707,0.866) # what units are these?
D1 <- par0[1:2] # parameters for kernel 1. this is behavior mode transit
D2 <- par0[3:4] # parameters for kernel 2. resident behavior mode
p <- par0[5:6] # not sure what these parameters are.. look like the diagonal of a 2x2 transition matrix.  


# Probably need to express kernel movement in terms of pixels per time step. The sparse matrix work likely renders this unnecessary, but going back to gausskern, it is. For example, if we have .25 degree and daily time step, what would the speed of the fish be when moving fast? 4 pixels/day?

K1 = as.cimg(gausskern(D1[1], D1[2], muadv = 0))
K2 = as.cimg(gausskern(D2[1], D2[2], muadv = 0))
P <- matrix(c(p[1],1-p[1],1-p[2],p[2]),2,2,byrow=TRUE)

# make all NA's very tiny for the convolution
# the previous steps may have taken care of this...
# L[L==0] = 1e-15
# L[is.na(L)] = 1e-15

# add a 'skip' index for missing days in the L.. 

#### Function based on Pedersen 2011, originally using sparse matrix multiplication.. editing for gaussian kernel convolution
# watch out for matrix dimensionality. even though dimensions lined up, data did not.. 
hmm.filter2 <- function(g,L,K1,K2,P){
  require(imager) # convolution function
  require(magic) # has a rotate function.. and isn't matlab
  ## Filter data to estimate locations and behaviour
  
  T <- dim(L)[1] # dimension of time 
  row <- dim(g$lon)[1] # nrows
  col <- dim(g$lon)[2] # ncols
  m <- 2 # Number of behavioural states
  
  pred <- array(0,dim=c(m,T,col,row)) # empty array for prediction step. ordering in col before row emulates lon before lat
  phi  <- array(0,dim=c(m,T,col,row)) # posterior (final) step array
  # Start in resident state at the known initial location
  phi[2,1,,]  <- L[1,,] # first position is known
  pred[2,1,,] <- L[1,,] # first position is known
  psi <- rep(0,T-1) # sum of the probability of both states at each step
  # Start filter iterations
  for(t in 2:T){
    # replace this part with older workflow using a gaussian kernel.. 
    # p1 <- as.vector(phi[1,t-1,,])
    # p2 <- as.vector(phi[2,t-1,,])
    # q1 <- as.vector(p1%*%K1)
    # q2 <- as.vector(p2%*%K2)
    
    p1 = as.cimg(t(phi[1,t-1,,]))
    p2 = as.cimg(t(phi[2,t-1,,]))
    q1 = convolve(p1, K1)
    q2 = convolve(p2, K2)
    
    # q1 = arot(t(as.matrix(q1)),3)
    # q2 = arot(t(as.matrix(q2)),3)
    q1 = t(as.matrix(q1))
    q2 = t(as.matrix(q2))
    
    # 	par(mfrow=c(1,2))
    # 	image(q1)
    # 	image(q2)
    
    # pred[1,t,,] <- matrix(P[1,1]*q1+P[2,1]*q2,row,col)
    # pred[2,t,,] <- matrix(P[1,2]*q1+P[2,2]*q2,row,col)
    
    # multiply by transition probability 
    pred[1,t,,] <- P[1,1]*q1+P[2,1]*q2
    pred[2,t,,] <- P[1,2]*q1+P[2,2]*q2
    
    sumL = sum(L[t,,])  
    if(sumL > 0){
      post1 <- pred[1,t,,]*L[t,,]
      post2 <- pred[2,t,,]*L[t,,]
    }else{
      post1 <- pred[1,t,,]
      post2 <- pred[2,t,,]
    }
    
    psi[t-1] <- sum(as.vector(post1), na.rm=T) + sum(as.vector(post2), na.rm=T)
    
    # remove NaNs... 
    # normalise (divide here by sum, not max)
    # 	post1 <- normalise(post1)
    # 	post2 <- normalise(post2)
    # 	post1[is.nan(post1)] = 0
    # 	post2[is.nan(post2)] = 0
    
    phi[1,t,,] <- post1/(psi[t-1]+1e-15)
    phi[2,t,,] <- post2/(psi[t-1]+1e-15)
  }
  list(phi=phi,pred=pred,psi=psi)
}

f = hmm.filter2(g,L,K1,K2,P)
res = apply(f$phi[1,,,],2:3,sum, na.rm=T)
image.plot(lon, lat, res/max(res), zlim = c(.05,1))

# Next up.... 
hmm.smoother <- function(f,K1,K2,P,plot=TRUE){
  ## Smoothing the filtered estimates
  ## The equations for smoothing are presented in Pedersen et al. 2011, Oikos, Appendix
  T <- dim(f$phi)[2]
  row <- dim(f$phi)[3]
  col <- dim(f$phi)[4]
  
  smooth <- array(0,dim=dim(f$phi))
  smooth[,T,,] <- f$phi[,T,,]
  for(t in T:2){
    RAT <- smooth[,t,,]/(f$pred[,t,,]+1e-15)
    #     Rp1 <- as.vector(K1 %*% as.vector(RAT[1,,]))
    #     Rp2 <- as.vector(K2 %*% as.vector(RAT[2,,]))
    
    
    p1 = as.cimg(t(RAT[1,,]))
    Rp1 <- convolve(p1, K1)
    p2 = as.cimg(t(RAT[2,,]))
    Rp2 <- convolve(p2, K2)
    
    Rp1 = t(as.matrix(Rp1))
    Rp2 = t(as.matrix(Rp2))
    
    if(plot){
      par(mfrow=c(1,2))
      image.plot(Rp1)
      image.plot(Rp2)
    }
    
    post1 <- matrix(P[1,1]*Rp1 + P[1,2]*Rp2,row,col)
    post2 <- matrix(P[2,1]*Rp1 + P[2,2]*Rp2,row,col)
    post1 <- post1 * f$phi[1,t-1,,]
    post2 <- post2 * f$phi[2,t-1,,]
    fac <- sum(as.vector(post1)) + sum(as.vector(post2))
    smooth[1,t-1,,] <- post1/fac
    smooth[2,t-1,,] <- post2/fac
  }
  smooth
}

s = hmm.smoother(f, K1, K2, P, plot = T)

sres = apply(s[1,,,],2:3,sum, na.rm=T)
image.plot(lon, lat, sres/max(sres), zlim = c(.05,1))

sres = apply(s[2,,,],2:3,sum, na.rm=T)
image.plot(lon, lat, sres/max(sres), zlim = c(.05,1))

# calculate track
# don't use this
calc.track(s, g)  # dimensions flipped...

# switch the dimensions in the calc.track.r... gives a weird output.. ON FIN LAND!
# this is either 1) right and we have to deal with the L and K elements or 2) the dimensions need adjusting..
distr = s
meanlat <- apply(apply(distr,c(2,4),sum)*as.vector(repmat(t(as.matrix(g$lat[,1])),T,1)),1,sum)
meanlon <- apply(apply(distr,c(2,3),sum)*as.vector(repmat(t(as.matrix(g$lon[1,])),T,1)),1,sum)


plot(meanlon, meanlat)
plot(countriesLow, add = T)

sres = apply(s,c(3,4), sum, na.rm=T)
image.plot(lon, lat, sres/max(sres), zlim = c(.01,1))
lines(meanlon, meanlat, pch=19, col=2)
plot(countriesLow, add = T)




##########
## END
##########


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
#CDB: [1] Error in diag(A) : no method for coercing this S4 class to a vector

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


