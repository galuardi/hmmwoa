#================================================================#
# HMM BFT script
#================================================================#
basedir=('C:/Users/ben/Documents/HMM-DGLM')
setwd(basedir)

load('C:/Users/ben/Documents/HMM-DGLM/DGLM_and_CONV.Rdata')
source('C:/Ben/SOFTWARE/RDRR/make.sstmat.r')
blsst = make.sstmat('C:/Users/ben/Documents/HMM-DGLM/SST/BLSST_08132007_09102008.nc', type="Blended")
blsst$DATA[is.na(blsst$DATA)] = NaN
load('C:/Users/ben/Documents/HMM-DGLM/SST/MWIR_0708_3day.RData')
source('C:/Users/ben/Documents/HMM-DGLM/HMMfunctions.r')
source('C:/Users/ben/Documents/HMM-DGLM/hmmfilter.r')
source('C:/Users/ben/Documents/HMM-DGLM/Viterbi.r')
source('C:/Users/ben/Documents/HMM-DGLM/resampleBathy.r')
# source('C:/Ben/toolbox/ls.objects.r')
source('C:/Users/ben/Documents/HMM-DGLM/make.btrack.alldates.r')

library(analyzepsat)
library(mvtnorm)
library(fields)
library(EBImage)
library(raster)

xtracks.hmm = NULL
for(i in 1:21){
 track = make.btrack.alldates(allfmat[[i]], psat[[i]])
 track$SST = apply(psat[[i]]$T,1,max,na.rm=T)
# track$SSTsd = apply(psat[[i]]$T,1,sd,na.rm=T)
 zidx = is.na(track$max_depth)
 track$max_depth[zidx] = 0
 ntrack = .fill.vals(track$SST)
 track[,11] = ntrack
 track$sstsd = getSSTsd(psat[[i]], surf = -11)
 names(track)[1:3] = c('Day','Month','Year')
 xtracks.hmm[[i]] = track
 rm(track, ntrack)
}


k=4

track = xtracks.hmm[[k]]
# track[,1:3] =rev(track[,1:3]) 
#track[290:294,11]=23.0489 

lims= c((range(track[,8],finite=T)),range(track[,9],finite=T))+c(-3,3,-3,3)

 L = make.like(blsst, bathy, track, lims, uselog = F)
 # LL = make.like(blsst, bathy, track, lims, uselog = T)
 
allpost = hmmfilter(L, track, ks = 29, sigma = 3)
allpost2 = hmmsmooth(allpost, track, 29, 3)

allpost2[allpost2<0] = 0

pidx = 107

par(mfrow=c(1,3))
image.plot(L$lon, L$lat,L$L[,,pidx])
# image.plot(LL$L[,,10])
image.plot(L$lon, L$lat,allpost[,,pidx])
image.plot(L$lon, L$lat,allpost2[,,pidx])

# most probable track
mpt = mpt.viterbi(allpost2,L, track, blsst, D = c(29,29), D2s = .09)

mpt2 = mpt.viterbi(allpost2,L, track, blsst, D = c(39,39), D2s = .09)


#=============================================================#
# Run the Viterbi Code and then do the following::::
#=============================================================#

vmpt = allpost2[,,1]
for(i in 2:(dim(allpost2)[3])){
 vmpt = allpost2[,,i]+vmpt
}

par(mfrow=c(1,2))
ftrack = xtracks.hmm[[k]]
plot.btrack(allnewbtracks[[k]], add=F, ci = F)
title('ukfsst/bathymetry')
ftrack[,8:9] = mpt
ftrack[1,8:9] = xtracks.hmm[[k]][1,8:9]
len = nrow(xtracks.hmm[[k]])
ftrack[len,8:9] = xtracks.hmm[[k]][len,8:9]

plot.btrack(ftrack, add=F)
title('HMM w/sst/bathymetry')

vmpt[log(vmpt)<as.numeric(quantile(log(vmpt),.5))] = 0
image(L$lon, L$lat, log(vmpt), col=tim.colors(100), add=T)
plot.btrack(ftrack, add=T)
contour(L$lon, L$lat, log(vmpt),  add=T)

par(mfrow=c(3,3))

for(i in 1:9){ 
image.plot(lon, lat, allpost2[,,i])
if (i==1|i==9)plot(map2$SP, add=T, col = 'khaki')
}

for(i in 20:28){ 
image.plot(L$lon, L$lat, allpost2[,,i])
if (i==20|i==28)plot(map2$SP, add=T, col = 'khaki')
}

save.image('Viterbi_tag19.Rdata')

# Contours for CI/day
makeCIhmm <- function(allpost2, L){
 out = list()
 for(i in 1:dim(allpost2)[3]){
	temp = allpost[,,i]
    temp[log(temp)<as.numeric(quantile(log(temp),.5,na.rm=T))] = 0
   out[[i]] = contourLines(L$lon,L$lat, temp)
 }
 out
}
 image.plot(L$lon-360,L$lat,allpost[,,2])
 lapply(makeCIhmm(allpost2,L)[[2]],lines)


