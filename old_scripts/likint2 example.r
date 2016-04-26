# calculate sd on a grid
load('C:/Users/benjamin.galuardi/Google Drive/Camrin-WOA/hmmwoa_files/HMM_WORK_LYDIA.Rdata')
load('~/Documents/WHOI/RData/WhiteSharks/2013/121325/lydia_hmm.RData')
## hopefully this has pdt,dat, and the other variables defined in the run script?

# code below taken from within calc.pdt.int()
udates <- unique(pdt$Date)
T <- length(udates)

pdt$MidTemp <- (pdt$MaxTemp + pdt$MinTemp) / 2

#L.pdt <- array(0, dim = c(dim(dat)[1:2], length(dateVec)))

#for(i in 1:T){
 i=10

  # define time based on tag data
  time <- udates[i]
  pdt.i <- pdt[which(pdt$Date == time),]
  
  #extracts depth from tag data for day i
  y <- pdt.i$Depth[!is.na(pdt.i$Depth)] 
  y[y<0] <- 0
  
  #extract temperature from tag data for day i
  x <- pdt.i$MidTemp[!is.na(pdt.i$Depth)]  
  
  # use the which.min
  depIdx = apply(as.data.frame(pdt.i$Depth), 1, FUN=function(x) which.min((x-depth)^2))
  woaDep <- depth[depIdx] 
  
  # make predictions based on the regression model earlier for the temperature at standard WOA depth levels for low and high temperature at that depth
  fit.low <- locfit(pdt.i$MinTemp ~ pdt.i$Depth)
  fit.high <- locfit(pdt.i$MaxTemp ~ pdt.i$Depth)
  n = length(depth[depIdx])
  
  pred.low = predict(fit.low, newdata = depth[depIdx], se = T, get.data = T)
  pred.high = predict(fit.high, newdata = depth[depIdx], se = T, get.data = T)
  
  # data frame for next step
  df = data.frame(low=pred.low$fit-pred.low$se.fit*sqrt(n)
                  , high=pred.high$fit+pred.high$se.fit*sqrt(n)
                  , depth = depth[depIdx])
  
 pdtMonth <- as.numeric(format(as.Date(pdt.i$Date), format='%m'))[1]
 dat.i = dat[,,,pdtMonth] #extract months climatology

    # calculate sd using Le Bris neighbor method and focal()
    sd.i = array(NA,dim=dim(dat.i))
    for(ii in 1:57){
      r = flip(raster(t(dat.i[,,ii])),2)
      #plot(r, col = tim.colors(100))
      f1 = focal(r, w=matrix(1/9,nrow=3,ncol=3), fun=sd)
      f1 = t(as.matrix(flip(f1,2)))
      #plot(f1, add=T)
      sd.i[,,ii] = f1
    } 
    
  # setup the likelihood array for each day. Will have length (dim[3]) = n depths
  #lik.pdt = array(NA, dim=c(dim(dat)[1], dim(dat)[2], length(depIdx)))
  
  #for (b in 1:length(depIdx)) {
   b=1 
    #calculate the likelihood for each depth level, b
    #lik.pdt[,,b] = likint(dat.i[,,b], df[b,1], df[b,2], sd.i[,,b])
   
   likint2 <- function(woa, woasd, minT, maxT){
     wlist = array(1e-6, dim=c(dim(woa)[1], dim(woa)[2], 2))
     wlist[,,1] = woa
     wlist[,,2] = woasd
     wlist[is.na(wlist)] = 1e-6
     as.matrix(aaply(wlist, 1:2, .fun = function(x) integrate(dnorm, lower = minT, upper = maxT , mean = x[1], sd = x[2])$value))
   }
   
   
   pdf('likint2 output.pdf',height=12,width=8)
   par(mfrow=c(3,1))
   image.plot(dat.i[,,b])
   title('Mean Temp')
   image.plot(sd.i[,,b])
   title('SD from Focal')
   image.plot(likint2(dat.i[,,b], sd.i[,,b], df[b,1], df[b,2]))
   title('likint2 output')
   dev.off()
   
<<<<<<< HEAD
   
   #-----------------------------------------------------------------------------#
   # Try using dplyr... which is supposedly much faster
   
   library(dplyr)
   
   likint3 <- function(woa, woasd, minT, maxT){
     # wlist = array(1e-6, dim=c(dim(woa)[1], dim(woa)[2], 2))
     wdf = data.frame(woa = as.vector(woa), sd = as.vector(woasd))
     wdf[is.na(wdf)] = 1e-6
#      wdf = add_rownames(wdf)
#      wrow = group_by(wdf, rowname)
     #      wlist[,,1] = woa
     #      wlist[,,2] = woasd
     #      wlist[is.na(wlist)] = 1e-6
     # matrix(do(wdf, .fun = function(x) integrate(dnorm, lower = minT, upper = maxT , mean = x[1], sd = x[2])$value), dim(woa)[1], dim(woa)[2])
     # res = wdf %>% rowwise() %>% do(integrate(dnorm, lower = minT, upper = maxT, mean = woa, sd = sd)$value)
     res = wdf %>% rowwise() %>% do(function(x) dnorm(10, mean = woa, sd = sd))
   }
   
   t1 = Sys.time()
   image.plot(likint3(dat.i[,,b], sd.i[,,b], df[b,1], df[b,2]))
   t2 = Sys.time() 
   t2-t1
   
   
=======
#--------------------------------------------------------#
# version where we only integrate within the range   
#   
# fuckin A.. this is correct after all. It's the likelihood of the cell given the temperature observed, not the other way round.. limits of integration are set according to the possible values of the cell, i.e., the obsevred value and calculated sd. The dnorm meanand sd are set from the tag measured mean (midTemp) and sd (adhoc sd)    
#--------------------------------------------------------#   
likint3 <- function(w, wsd, minT, maxT){
  midT = (maxT+minT)/2
  Tsd = (maxT-minT)/4
  widx = w>=minT&w<=maxT&!is.na(w)
  widxv = as.vector(widx)
  wdf = data.frame(w = as.vector(w[widx]), wsd = as.vector(wsd[widx]))
  wdf$wsd[is.na(wdf$wsd)] = 0
  # wint = apply(wdf, 1, function(x) pracma::integral(dnorm, minT, maxT, mean = x[1], sd = x[2]))
  wint = apply(wdf, 1, function(x)integrate(dnorm, x[1]-x[2], x[1]+x[2], mean = midT, sd = Tsd*2)$value) 
  w = w*0  
  w[widx] = wint
  w
} 
    
w = dat.i[,,b]
wsd = sd.i[,,b]


t1 = Sys.time()
lik2 = likint2(dat.i[,,b], sd.i[,,b], df[b,1], df[b,2])
t2 = Sys.time()
l2time = t2-t1

t1 = Sys.time()
lik3 = likint3(dat.i[,,b], sd.i[,,b], df[b,1], df[b,2])
t2 = Sys.time()
l3time = t2-t1

pdf('likint update.pdf',height=12,width=12)
par(mfrow=c(2,2))
image.plot(dat.i[,,b])
title('Mean Temp')
image.plot(sd.i[,,b])
title('SD from Focal')
image.plot(lik2)
title(paste0('likint2 output ', round(l2time),3), ' seconds')
image.plot(lon,lat,lik3[[1]])
points(-80,30)
title(paste0('likint3 output ', round(l3time,3), ' seconds'))
dev.off()

pdf('compare likint.pdf',height=12,width=8)
par(mfrow=c(2,1))
for(b in 1:8){
  image.plot(lon,lat,lik2.pdt[,,b])
  plot(countriesLow,add=T)
  title(paste('lik2 - depth ',woaDep[b]))
  image.plot(lon,lat,lik3.pdt[,,b])
  plot(countriesLow,add=T)
  title('lik3')
}
lik2.pdt.prod <- apply(lik2.pdt, 1:2, prod, na.rm=F)
lik3.pdt.prod <- apply(lik3.pdt, 1:2, prod, na.rm=F)
image.plot(lon,lat,lik2.pdt.prod)
plot(countriesLow,add=T)
title('product of lik2 layers')
image.plot(lon,lat,lik3.pdt.prod)
plot(countriesLow,add=T)
title('product of lik3 layers')
dev.off()




>>>>>>> origin/filter-2.0
