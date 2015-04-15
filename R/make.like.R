make.like <-
function(sst, bathy, fmat, lims = c(-85,-65,30,45), uselog = F){
require(raster)
#=======================================================#
# Make the Likelihood array for SST
#=======================================================#

lon = sst$lon
lat = sst$lat

xmin = which.min(((lon)-lims[1])^2)
xmax = which.min(((lon)-lims[2])^2)
ymin = which.min(((lat)-lims[3])^2)
ymax = which.min(((lat)-lims[4])^2)

tempL = numeric(length = c(dim(sst$DATA)[1]*dim(sst$DATA)[2]*nrow(fmat)))
dim(tempL) = c(dim(sst$DATA)[1], dim(sst$DATA)[2], nrow(fmat))
tempL = tempL[xmin:xmax, ymin:ymax,]
land = sst$DATA[xmin:xmax, ymin:ymax, 1]*0+2
land[is.nan(land)] = 0
land[land==2] = 1

time1 = Sys.time()


for(i in 2:(nrow(fmat)-1)){   
	sigma  = fmat$sstsd[i]
   if(class(sst$sstdates)[1]=="POSIXct"){
	tdate = ISOdate(fmat$Year[i],fmat$Month[i], fmat$Day[i])
   	zidx = which.min((as.numeric(tdate) - as.numeric(sst$sstdates))^2)
     }else{
    	tdate = mdy.date(fmat$Month[i], fmat$Day[i], fmat$Year[i])
    	zidx = which.min((tdate - sst$sstdates)^2)
	}
   if(uselog==T){
		tempL[,,i] = normalise(log.like.sst(sst$DATA[xmin:xmax, ymin:ymax,zidx], as.numeric(fmat[i,'SST']), sigma = sigma, log = F))
   }else{
		# tempL[,,i] = like.sst(as.numeric(fmat[i,'SST']),sst$DATA[xmin:xmax, ymin:ymax,zidx], sigma = sigma)		
		tempL[,,i] = like.sst.tnorm(ingrid = sst$DATA[xmin:xmax, ymin:ymax,zidx], datax = as.numeric(fmat[i,'SST']), sigma = sigma)			
		}
  
}

sstL = list(lon = sst$lon[xmin:xmax], lat = sst$lat[ymin:ymax], sstL = tempL)
print(paste(Sys.time()-time1, 'seconds to make sst likelihood'))
rm(tempL)

#=======================================================#
# Make the Likelihood array for Depth
#=======================================================#
bath2 = resamp.bath(sstL, bathy)
tempB = sstL$sstL*NaN

for(i in 2:(nrow(fmat)-1)){
   # print(i)
   zidx = bath2$data<=as.numeric(fmat[i,10])   
   bb = tempB[,,i]
   bb[zidx==T] = 1
   bb[zidx==F] = 0
   tempB[,,i] = bb#  *sstL$sstL[,,i]  # do this one at a time.. something wrong with 3Dmatrix multiplication
   }
#bathL = list(lon = bath2$lon, lat = bath2$lat, bathL = tempB)

# tempB = tempB*land


print(paste(Sys.time()-time1, 'seconds to make depth likelihood'))
#rm(tempB)
list(lon = sstL$lon, lat = sstL$lat, L = sstL$sstL*tempB)
# bathL
}
