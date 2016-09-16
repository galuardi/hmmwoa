## SRSS likelihood example - Blue shark 256
load('srss_example.RData')

# you have all the variables as loaded and formatted below
# of interest will be 'light' and 'spot'
# for the raw light data and the 'true' spot positions










#=================
## loaded and formatted this data as below
#=================

setwd('~/Documents/WHOI/Data/Blues/2015/141256/')
light <- read.table('141256-LightLoc.csv',sep=',',header=T,blank.lines.skip=F,skip=2)

# LOAD TAG/POP 
iniloc <- data.frame(matrix(c(13, 10, 2015, 41.575, -69.423, 
                              24, 2, 2016, 26.6798, -69.0147), nrow = 2, ncol = 5, byrow = T))
colnames(iniloc) = list('day','month','year','lat','lon')
tag <- as.POSIXct(paste(iniloc[1,1], '/', iniloc[1,2], '/', iniloc[1,3], sep=''), format = '%d/%m/%Y')
pop <- as.POSIXct(paste(iniloc[2,1], '/', iniloc[2,2], '/', iniloc[2,3], sep=''), format = '%d/%m/%Y')

# VECTOR OF DATES FROM DATA. THIS WILL BE THE TIME STEPS, T, IN THE LIKELIHOODS
dateVec <- as.Date(seq(tag, pop, by = 'day')) 

# LOAD SPOT LOCATIONS FOR "TRUE" TRACK
all <- read.table('~/Documents/WHOI/RData/sharkSiteData/AllArgosData.csv', sep=',', header=T)
all <- all[which(all$ptt == 141268),]
alldts <- as.POSIXct(all$date,format=findDateFormat(all$date))
all <- all[which(alldts >= tag & alldts <= pop),]

