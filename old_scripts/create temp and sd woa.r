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
sd <- var.def.ncdf('std dev', 'degrees C', list(dimX,dimY,dimZ,dimT), mv, longname = 'Standard Deviation of Temperature') 

# create new ncdf file
ncnew <- create.ncdf('../WOA_25deg/global/woa13_25deg_global_combine_again.nc', sd, verbose=T) # create new ncdf on disk
#ncnew <- create.ncdf('../WOA_25deg/global/woa13_25deg_global_new.nc', list(temp,sd), verbose=T) # create new ncdf on disk

# loop through monthly climatology
for(i in 1:12){
  nc = open.ncdf(paste(woa.dir, ncfiles[i], sep = ''))
  
  # extract mean temps and write to temp var
  sd.temp = get.var.ncdf(nc, 't_sd', start = c(1, 1, 1, 1), count = c(-1, -1, -1, 1))
  
  # put the extracted sd temp data into newly created 4d netcdf (i = month)
  t <- Sys.time()
  put.var.ncdf(ncnew, sd, sd.temp, start = c(1,1,1,i), count = c(-1,-1,-1,1), verbose=T)
  Sys.time() - t # takes about 6.5 mins on the mac
  
}
close.ncdf(ncnew)

# errors in here somewhere. looks like all but the first month are messed up

