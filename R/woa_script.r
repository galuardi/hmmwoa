
#--------------------------------------------------------------------
# format PDT for ML
#

#load pdtHeader names
#load('/Users/Cam/Documents/WHOI/RCode/pdtMatch/pdtHeader.RData')
#load('/Users/Cam/Documents/WHOI/RCode/pdtMatch/removePacific.RData')

#load convert.pdt function
source('/Users/Cam/Documents/WHOI/RCode/pdtMatch/convert.pdt.r')
source('/Users/Cam/Documents/WHOI/RCode/pdtMatch/extract.woa.r')
source('/Users/Cam/Documents/WHOI/RCode/pdtMatch/createGrid.r')
source('/Users/Cam/Documents/WHOI/RCode/pdtMatch/matchProfiles.r')
source('/Users/Cam/Documents/WHOI/RCode/pdtMatch/findDateFormat.r')
source('/Users/Cam/Documents/WHOI/RCode/pdtMatch/removePacific.r')

#type = read.table('/Users/Cam/Documents/WHOI/RData/BaskingSharks/tagtype.csv',sep=',',header=T)

## architecture for matching basko pdt profiles
#tagyear = 2011
#ptt = 96002
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
  limits = c(-90,-30,10,50) # (min long, max long, min lat, max lat)
  #limits = c(-100,-97.75,-30,-27.75) # (min long, max long, min lat, max lat)
  
  
  nc.dir = '/Users/Cam/Documents/WHOI/RData/pdtMatch/WOA_25deg/global/'
  
  return.woa = extract.woa(nc.dir, limits, resolution = 'quarter', pdt)
  dat = return.woa$dat; lon = return.woa$lon; lat = return.woa$lat; depth = return.woa$depth
  
  # eliminate Pacific from matching
  dat = removePacific(dat, lat, lon)
  
  # generates grid of lat x lon at 1 or 1/4 deg resolution
  returnGrid = createGrid(dat, lat, lon)
  gridPts = returnGrid$gridPts; cellidx = returnGrid$cellidx
  
  #Rprof('test.txt')
  # perform matching
  #system.time(
  mse = matchProfiles(pdt, dat, lat, lon, gridPts,cellidx)
  #)
  
  save.image(paste('/Users/Cam/Documents/WHOI/RData/',species,'/',tagyear,'/',ptt,'/',ptt,'-match.RData',sep=''))
  #Rprof(NULL)
  #summaryRprof('test.txt')
  
  #=================
  # END
  #=================
  cat('Finished ',ptt)
  
}



library(fields)
library(maptools)
pdf(file='test.pdf',height=11,width=8)
#par(mfrow=c(3,1))
#par(mar=c(1,3,2,3.5))
#for(i in 1:3){  
i = 1
image.plot(lon, lat, mse[,,5])
world(col='black', fill=T, add=T)
title(paste('MSE ',ptt,' day ',i))
#  grid (); box(); axis(1); axis(2)
#legend('topleft', paste('mse for uid ', i,sep=''), cex=1.3, bg='white',text.font=2)
#}


dev.off()

