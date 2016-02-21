# create netcdf of all months combined for global dataset
# need monthly mean/sd of woa temperature climatology
require(raster); require(ncdf)

# set dir containing monthly climatology files
woa.dir = '/Users/Cam/Documents/WHOI/RData/pdtMatch/WOA_25deg/'
setwd(woa.dir)
ncfiles = dir(woa.dir, pattern = '.nc')

# read one
nc = open.ncdf(paste(woa.dir, ncfiles[2], sep = ''))

# retrieve var bounds from global nc
lon = get.var.ncdf(nc, 'lon')
lat = get.var.ncdf(nc, 'lat')
depth = get.var.ncdf(nc, 'depth')

# define new vars for ncdf
crs <- "+proj=longlat +datum=WGS84 +ellps=WGS84"
dimX <- dim.def.ncdf('Longitude', 'degrees E', lon) # x dimension
dimY <- dim.def.ncdf('Latitude', 'degrees N', lat) # y dimension
dimZ <- dim.def.ncdf('Depth', 'meters', depth) # z dimension
dimT <- dim.def.ncdf('Time', 'months', c(1:12), unlim=FALSE) # time dimension

mv <- -9999
# define variable of interest. here it's temp and sd of temp in 4d
temp <- var.def.ncdf('temp', 'degrees C', list(dimX,dimY,dimZ,dimT), mv, longname = 'Mean Temperature') 
sd <- var.def.ncdf('std dev', 'degrees C', list(dimX,dimY,dimZ,dimT), mv, longname = 'Standard Deviation of Temperature',prec='double') 

# create new ncdf file
ncnew <- create.ncdf('../WOA_25deg/global/woa13_25deg_global_new.nc', sd, verbose=T) # create new ncdf on disk
#ncnew <- create.ncdf('../WOA_25deg/global/woa13_25deg_global_new.nc', list(temp,sd), verbose=T) # create new ncdf on disk

# loop through monthly climatology
for(i in 1:12){
  nc = open.ncdf(paste(woa.dir, ncfiles[i], sep = ''))
  
  # extract mean temps and write to temp var
  mu.temp = get.var.ncdf(nc, 't_an', start = c(1, 1, 1, 1), count = c(length(lon), length(lat), length(depth), 1))
  # put the extracted temp data into newly created 4d netcdf (i = month)
  t <- Sys.time()
  put.var.ncdf(ncnew, temp, mu.temp[,,1], start = c(1,1,1,i), count = c(-1,-1,1,1), verbose=T)
  Sys.time() - t # takes about 6.5 mins on the mac
  
  # focal calc on mean temp and write to sd var
  list.r <- list(x = lon, y = lat, z = mu.temp)
  ex <- extent(list.r)
  
  # and loop through depth levels
  f.arr <- array(NA, dim=c(length(lon),length(lat),57))
  t <- Sys.time()
  for(ii in 1:57){
    # create and orient raster
    r <- raster(t(list.r$z[,,ii]), xmn=ex[1], xmx=ex[2], ymn=ex[3], ymx=ex[4], crs)
    r <- flip(r, direction = 'y')
    # focal matrix and calculation
    w = matrix(1/9, nrow = 3, ncol = 3)
    f <- focal(r, w, function(x) sd(x))
    
    # put results in an array
    f.arr[,,ii] <- t(as.matrix(flip(f,direction='y')))
    put.var.ncdf(ncnew, sd, f.arr, start = c(1,1,1,1), count = c(-1,-1,-1,1),verbose=T)
    print(paste(round(ii/57*100),'%',sep=''))
  }
  Sys.time() - t # about 16 mins on the mac w/o put.var inside
  
  
  # write the sd var
  put.var.ncdf(ncnew, sd, f.arr[,,1], start = c(1,1,1,i), count = c(-1,-1,1,1))
  
  
}

close.ncdf(ncnew)



#===========================
#===========================

# create netcdf of all months combined for global dataset
# need monthly mean/sd of woa temperature climatology
lonVector = seq (-180, 179.75, by = .25) # define x bounds vector
latVector = seq(-90, 89.75, by = .25) # define y bounds vector
dimX <- dim.def.ncdf('Longitude', 'degrees E', lonVector) # x dimension
dimY <- dim.def.ncdf('Latitude', 'degrees N', latVector) # y dimension
dimZ <- dim.def.ncdf('Depth', 'meters', depth) # z dimension
dimT <- dim.def.ncdf('Time', 'months', c(1:12), unlim=FALSE) # time dimension

# define variable of interest. here it's temp in 4d
temp <- var.def.ncdf('temp', 'degrees C', list(dimX,dimY,dimZ,dimT), 'NA', longname = 'Mean Temperature') 
sd <- var.def.ncdf('std dev', 'degrees C', list(dimX,dimY,dimZ,dimT), 'NA', longname = 'Standard Deviation of Temperature') 

ncnew <- create.ncdf('woa13_25deg_global_new.nc', temp, verbose=T) # create new ncdf on disk

# loop through the new ncdf and put in data from monthly ncdf files
for(i in 1:12){
  nc = open.ncdf(ncfiles[i]) # open monthly nc file
  # extract temp data from monthly nc according to dims x,y,z,month(i)
  dat = get.var.ncdf(nc, 't_an', start = c(xmin,ymin,1,1), count = c(xlen,ylen,57,1))
  dat[dat<(-3)] = mv
  dat[dat>45] = mv
  # put the extracted temp data into newly created 4d netcdf (i = month)
  put.var.ncdf(ncnew, temp, dat, start = c(xmin,ymin,1,i), count = c(xlen,ylen,57,1))
}

close.ncdf(ncnew) # write out new netcdf



# smooth/interpolate sd for WOA13

woa.dir = '/Users/Cam/Documents/WHOI/RData/pdtMatch/WOA_25deg/'

# load global nc
ncfiles = dir(woa.dir, pattern = '.nc')

for(i in 1:12){
  nc = open.ncdf(paste(woa.dir, ncfiles[2], sep = ''))
  
  if(i==1){
    # retrieve var bounds from global nc
    lon = get.var.ncdf(nc, 'lon')
    lat = get.var.ncdf(nc, 'lat')
    depth = get.var.ncdf(nc, 'depth')
    crs <- "+proj=longlat +datum=WGS84 +ellps=WGS84"
    
  }
  
  temp = get.var.ncdf(nc, 't_an', start = c(1, 1, 1, 1), count = c(1440, 720, 57, 1))
  
  list.r <- list(x = lon, y = lat, z = temp)
  ex <- extent(list.r)
  r <- stack(t(list.r$z), xmn=ex[1], xmx=ex[2], ymn=ex[3], ymx=ex[4], crs)
  r <- flip(r, direction = 'y')
  
  f = focal(r, w = matrix(1/9, nrow = 3, ncol = 3), fun = sd)
  
  # put it in an array?
  
}

image.plot(lon,lat,dat[,,1])


#r = flip(raster(t(dat[,,1])),2)
plot(r, col = tim.colors(100))
f1 = focal(r, w=matrix(1/9,nrow=3,ncol=3), fun=sd)
plot(f1)

fmat = list()

for(i in 1:57){
  r = flip(raster(t(dat[,,i])),2)
  plot(r, col = tim.colors(100))
  f1 = focal(r, w=matrix(1/9,nrow=3,ncol=3), fun=sd)
  plot(f1, add=T)
  fmat[[i]] = f1
}

s <- raster::as.array(r, transpose = T)



# INTERPOLATE GAPS OUT OF WOA SD USING LOCFIT

#xp <- seq(min(lon),max(lon),length.out=(length(lon)-1)*4+1)
#yp <- seq(min(lat),max(lat),length.out=(length(lat)-1)*4+1)
dat1 <- dat[,,1]
m <- matrix(t(dat1), nrow=length(lat), dimnames=list(lat,lon))
mm <- melt(m)
colnames(mm) = list('lat','lon','sd')

new <- matrix(NA, nrow=length(lat), ncol=length(lon), dimnames=list(lat,lon))
new <- melt(new)
colnames(new) = list('lat','lon','sd.m')

# do the regression
fit = locfit(sd~lon:lat, data=mm, alpha = .025)  # try different values of alpha

# predict from the model
z.pred<-predict(fit, newdata=new)

# then format and plot
z.mat = t(matrix(z.pred, nrow = length(lat), ncol = length(lon)))
image.plot(lon, lat, z.mat,zlim=c(0,max(z.mat,na.rm=T)),xlim=c(-100,-30),ylim=c(10,50))
#contour(xp-360, yp, z.mat, add=T,col='black')
world(add=T, col='grey', fill=T)

image.plot(lon, lat, dat1,xlim=c(-100,-30),ylim=c(10,50))

#======================

# create netcdf of all months combined for global dataset
# need monthly mean/sd of woa temperature climatology
lonVector = seq (-180, 179.75, by = .25) # define x bounds vector
latVector = seq(-90, 89.75, by = .25) # define y bounds vector
dimX <- dim.def.ncdf('Longitude', 'degrees E', lonVector) # x dimension
dimY <- dim.def.ncdf('Latitude', 'degrees N', latVector) # y dimension
dimZ <- dim.def.ncdf('Depth', 'meters', depth) # z dimension
dimT <- dim.def.ncdf('Time', 'months', c(1:12), unlim=FALSE) # time dimension

# define variable of interest. here it's temp in 4d
temp <- var.def.ncdf('temp', 'degrees C', list(dimX,dimY,dimZ,dimT), 'NA', longname = 'Mean Temperature') 
sd <- var.def.ncdf('std dev', 'degrees C', list(dimX,dimY,dimZ,dimT), 'NA', longname = 'Standard Deviation of Temperature') 

ncnew <- create.ncdf('woa13_25deg_global_new.nc', temp, verbose=T) # create new ncdf on disk

# loop through the new ncdf and put in data from monthly ncdf files
for(i in 1:12){
  nc = open.ncdf(ncfiles[i]) # open monthly nc file
  # extract temp data from monthly nc according to dims x,y,z,month(i)
  dat = get.var.ncdf(nc, 't_an', start = c(xmin,ymin,1,1), count = c(xlen,ylen,57,1))
  dat[dat<(-3)] = mv
  dat[dat>45] = mv
  # put the extracted temp data into newly created 4d netcdf (i = month)
  put.var.ncdf(ncnew, temp, dat, start = c(xmin,ymin,1,i), count = c(xlen,ylen,57,1))
}

close.ncdf(ncnew) # write out new netcdf

# now try it
filename <- 'woa13_25deg_global.nc'
nc.try <- open.ncdf( filename )
print(paste("File",filename,"contains",nc$ndims,"dimensions"))
dat = get.var.ncdf(nc.try, 'temp', start = c(xmin,ymin,1,2), count = c(xlen,ylen,57,1))



dat = get.var.ncdf(nc, 'Temperature', start = c(xmin,ymin,1,1), count = c(xlen,ylen,57,1))



close.ncdf(nc)

