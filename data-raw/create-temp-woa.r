#=================================
## HOW WE COMPILED THE GLOBAL MONTHLY WOA DATA
## AND ADDED THEM TO PACKAGE AS EXAMPLE DATA
#=================================

# create netcdf of all months combined for global dataset
# need monthly mean/sd of woa temperature climatology
require(raster); require(ncdf)

# set dir containing monthly climatology files
woa.dir = '/Users/Cam/Documents/WHOI/RData/pdtMatch/WOA_1deg/'
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
#sd <- var.def.ncdf('std dev', 'degrees C', list(dimX,dimY,dimZ,dimT), mv, longname = 'Standard Deviation of Temperature') 

# create new ncdf file
ncnew <- create.ncdf('../WOA_1deg/global/woa13_1deg_global_meantemp.nc', temp, verbose=T) # create new ncdf on disk
#ncnew <- create.ncdf('../WOA_25deg/global/woa13_25deg_global_new.nc', list(temp,sd), verbose=T) # create new ncdf on disk

# loop through monthly climatology
for(i in 1:12){
  nc = open.ncdf(paste(woa.dir, ncfiles[i], sep = ''))
  
  # extract mean temps and write to temp var
  dat = get.var.ncdf(nc, 't_an', start = c(1, 1, 1, 1), count = c(-1, -1, -1, 1))
  
  # put the extracted temp data into newly created 4d netcdf (i = month)
  t <- Sys.time()
  put.var.ncdf(ncnew, 'temp', dat, start = c(1,1,1,i), count = c(-1,-1,-1,1))
  Sys.time() - t # takes about 6.5 mins on the mac
  
}
close.ncdf(ncnew)

# finally, i ran the resulting ncdf through extract.woa to subset it relative to a shark's example data to distribute with the package
woa.dir <- '/Users/Cam/Documents/WHOI/RData/pdtMatch/WOA_1deg/global/woa13_1deg_global_meantemp.nc'
return.woa <- extract.woa(woa.dir, sp.lim, resolution = 'one')
woa.one <- list(watertemp = woa.one$dat, lon = woa.one$lon, lat = woa.one$lat, depth = woa.one$depth)
setwd('~/Documents/WHOI/RCode/hmmwoa/')
devtools::use_data(woa.one)


