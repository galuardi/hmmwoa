library(analyzepsat)
library(lubridate)
library(raster)

load('C:/Users//benjamin.galuardi/Google Drive//SCDNR/allfits-alltracks-allpar-psat.Rdata')

# work with tag # 37078

df = cbind(allfits[[2]]$date, allfits[[2]]$nominal.track, allfits[[2]]$SST[,1])
names(df) = c('year','month','day','lon','lat','SST')

dates = as.POSIXct(strptime(paste(df[,1], df[,2], df[,3], sep='-'), '%Y-%m-%d', tz='UTC'))

lsst = data.frame(dates, lon = df$lon, sst = df$SST)

iniloc = df[1, 4:5]

par0 = c(8.908,10.27,1.152,0.0472,0.707,0.866)  # from Martins sphmm.demo function
# should be kernel 1 mean, kernel 1 sd, kernel 2 mean, kernel 2 sd, ?, ? Bmode??

sstf = open.ncdf('C:/Users/benjamin.galuardi/Google Drive/SCDNR/ANALYSIS/SSTNEW (1)/37078-SST/oisst.nc')

sst = get.var.ncdf(sstf, 'sst')

sstr = stack('C:/Users/benjamin.galuardi/Google Drive/SCDNR/ANALYSIS/SSTNEW (1)/37078-SST/oisst.nc')

#-----------------------------------------------------#
# slight rewrite of Cam's (rewrite) function


lik.locs <- function(locs,iniloc,g){
  ## Calculate the "data" likelihood, i.e. the likelihood field for each observation
  
  T <- length(locs$lon)
  row <- dim(g$lon)[1]
  col <- dim(g$lon)[2]
  
  L <- array(0,dim=c(col, row, T + 2))
  
  # Initial location is known
  ilo <- which.min(abs(g$lon[1,]-iniloc$lon[1]))
  ila <- which.min(abs(g$lat[,1]-iniloc$lat[1]))
  L[ilo, ila, 1] <- 1
  
  # Calculate data likelihood
  # SD for light-based longitude from Musyl et al. (2001)
  sl.sd <- 35/111 # Converting from kms to degrees
  
  for(i in 2:(T-1)){
#     if(locs$Type[t] == 'GPS'){
#       
#       glo <- which.min(abs(g$lon[1,]-locs$Longitude[t]))
#       gla <- which.min(abs(g$lat[,1]-locs$Latitude[t]))
#       L[glo, gla, (t+1)] <- 1
#       
#     } else if(locs$Type[t] == 'Argos'){
#       
#       alo <- which.min(abs(g$lon[1,]-locs$Longitude[t]))
#       ala <- which.min(abs(g$lat[,1]-locs$Latitude[t]))
#       L[alo, ala, (t+1)] <- 1
#       
#     } else if(locs$Type[t] == 'GPE'){
      L.light <- dnorm(t(g$lon), locs$lon[i], sl.sd) # Longitude data
      L[,, i] <- L.light
    }

  # End location is known
  elo <- which.min(abs(g$lon[1,]-iniloc$lon[2]))
  ela <- which.min(abs(g$lat[,1]-iniloc$lat[2]))
  L[elo, ela, T] <- 1
  
  L
  
}


# use my truncated normal likelihood for one layer.. 
# and rasters..
sst1 = sstr[[1]]
sst1@extent[1] = sst1@extent[1]-360
sst1@extent[2] = sst1@extent[2]-360
# tmat = (matrix(dtruncnorm(as.matrix(sst1), a=10, b=40, mean = 30.07, sd = 2.5), 20, 21))
tmat = (matrix(dnorm(as.matrix(sst1),mean = 30.07, sd = 1), 20, 21))
rmat = sst1
values(rmat) = tmat
plot(rmat, zlim = c(.05,1))
contour(sst1, add=T)
world(add=T)

#-----------------------------------------------------#
# begin SPHMM stuff with modifications
#-----------------------------------------------------#

# source('../sim/sstdb.R')  ## Need an SST matrix... 
require(Matrix) # For expm, matrix exponential and sparsity
# 
# lsst <- read.table('../sim/lsstdata.csv',header=TRUE,sep=',')
# tr <- read.table('../sim/truedata.csv',header=TRUE,sep=',')
# iniloc <- read.table('../sim/iniloc.csv',header=TRUE,sep=',')

## Calculate time step
dt <- date2time(lsst$dates[2])-date2time(lsst$dates[1])

## Setup discrete spatial grid
ngrid <- c(50,49)

# need to rewrite setup.grid! SPHMM uses the sstdb.r for Latitude limits.. 
g <- setup.grid(lsst,ngrid)
row <- dim(g$lon)[1]
col <- dim(g$lon)[2]

gr = raster(nrows = row, ncols = col, xmn = min(g$lon), xmx = max(g$lon), ymn = min(g$lat), ymx = max(g$lat))

# make liklihood with Cam/Ben modified function
myL = lik.locs(df, iniloc, g)

## Compute likelihood for observations
L <- data.lik(lsst,iniloc,g)

# plot the longitude likelihood
par(mfrow=c(1,2))
image(dates, g$lat[,1], apply(L, 1:2, sum, na.rm=T), col = terrain.colors(100))  # latitude/time
image(dates, g$lon[1,], apply(L, c(1,3), sum, na.rm=T), col = terrain.colors(100))  # longitude/time

# plot from modified function
image(g$lon[1,]-360, g$lat[,1], apply(lik.locs(df, iniloc, g), 1:2, sum))
world(add=T)

## Number of time steps
T <- length(lsst$lon)

## Fixed parameter values
D1 <- par0[1:2]
D2 <- par0[3:4]
p <- par0[5:6]


## Setup transition matrices
G1 <- make.kern(D1,g,dt)
K1 <- uniformization(G1,dt)
G2 <- make.kern(D2,g,dt)
K2 <- uniformization(G2,dt)
P <- matrix(c(p[1],1-p[1],1-p[2],p[2]),2,2,byrow=TRUE)

