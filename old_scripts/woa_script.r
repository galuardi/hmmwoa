
#--------------------------------------------------------------------
# format PDT for ML
#

#load pdtHeader names
#load('/Users/Cam/Documents/WHOI/RCode/pdtMatch/pdtHeader.RData')
#load('/Users/Cam/Documents/WHOI/RCode/pdtMatch/removePacific.RData')

#load convert.pdt function
source('/Users/Cam/Documents/WHOI/RCode/pdtMatch/convert.pdt.r')
#source('/Users/Cam/Documents/WHOI/RCode/pdtMatch/extract.woa.r')
source('/Users/Cam/Documents/WHOI/RCode/pdtMatch/createGrid.r')
source('/Users/Cam/Documents/WHOI/RCode/pdtMatch/matchProfiles.r')
source('/Users/Cam/Documents/WHOI/RCode/pdtMatch/findDateFormat.r')
source('/Users/Cam/Documents/WHOI/RCode/pdtMatch/removePacific.r')

#type = read.table('/Users/Cam/Documents/WHOI/RData/BaskingSharks/tagtype.csv',sep=',',header=T)

## architecture for matching basko pdt profiles
#tagyear = 2011
ptt = 121325
species = 'WhiteSharks'
#27:28, no 29, skipped 30,skipped 39
#for (i in 40:45){
#  ptt = type[i,1]; tagyear = type[i,8]
#  setwd(paste('/Users/Cam/Documents/WHOI/Data/', species, '/', tagyear, '/', ptt, '/', sep=''))
  setwd('~/Documents/WHOI/RData/WhiteSharks/2013/121325/')  
  
  # read in the WC PSAT data
  tagfile = paste(ptt,'-PDTs.csv', sep='')
  pdt.data = read.table(tagfile, sep=',', skip=0, header = T, blank.lines.skip = F)
  
  #pdt = convert.pdt(pdt.data)
  pdt <- extract.pdt(pdt.data)
  
  # format woa data
  # basically function to open connection to each month's woa data
  # and extract the temp data based on user's required lat/lon
  
  # set limits of interest
  limits = c(-70,-15,20,60) # (min long, max long, min lat, max lat)
  #limits = c(-100,-97.75,-30,-27.75) # (min long, max long, min lat, max lat)
  
  
  nc.dir = '/Users/Cam/Documents/WHOI/RData/pdtMatch/WOA_25deg/global/'
  
  return.woa = extract.woa(nc.dir, limits, resolution = 'quarter')
  dat = return.woa$dat; lon = return.woa$lon; lat = return.woa$lat; depth = return.woa$depth
  
  # eliminate Pacific from matching
  dat = removePacific(dat, lat, lon)
  
  # generates grid of lat x lon at 1 or 1/4 deg resolution
  #returnGrid = createGrid(dat, lat, lon)
  #gridPts = returnGrid$gridPts; cellidx = returnGrid$cellidx
  
  # perform matching
  lik = calc.pdt(pdt, dat, lat, lon)

  #spot <- read.table('lydia_spot.csv',sep=',',header=T)
  
  plot.woa(lik, return.woa, '106795_woa.pdf', pdt = pdt, write.dir = getwd())

  
  lon <- limits[c(1, 2)]
  lat <- limits[c(3, 4)]
  
  ohc.dir <- paste('~/Documents/WHOI/RData/HYCOM/', ptt, '/',sep = '')
  
  for(i in 1:length(udates)){
    time <- as.Date(udates[i])
    repeat{
      get.hycom(lon,lat,time,filename=paste(ptt,'_',time,'.nc',sep=''),download.file=TRUE,dir=ohc.dir) # filenames based on dates from above
      tryCatch({
        err <- try(open.ncdf(paste(ohc.dir,ptt,'_',time,'.nc',sep='')),silent=T)
      }, error=function(e){print(paste('ERROR: Download of data at ',time,' failed. Trying call to server again.',sep=''))})
      if(class(err) != 'try-error') break
    }
  }
  
  likelihood <- calc.ohc(pdt, ohc.dir = ohc.dir, ptt = 106795, sdx = 10)
  
  plot.ohc(lik = likelihood, filename='lyd lik.pdf', write.dir = getwd())
  
  