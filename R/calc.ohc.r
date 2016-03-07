calc.ohc <- function(pdt, isotherm = '', ohc.dir, g, dateVec, raster = 'stack'){
  # compare tag data to ohc map and calculate likelihoods
  
  #' @param: pdt is variable containing tag-collected PDT data
  #' @param: isotherm is default '' in which isotherm is calculated
  #' on the fly based on daily tag data. Otherwise, numeric isotherm
  #' constraint can be specified (e.g. 20).
  #' @param: ohc.dir is local directory where get.hycom downloads are
  #' stored.
  #' @return: likelihood is array of likelihood surfaces representing
  #' matches between tag-based ohc and hycom ohc maps
  
  # constants for OHC calc
  cp <- 3.993 # kJ/kg*C <- heat capacity of seawater
  rho <- 1025 # kg/m3 <- assumed density of seawater
  
  # calculate midpoint of tag-based min/max temps
  pdt$MidTemp <- (pdt$MaxTemp + pdt$MinTemp) / 2
  
  # get unique time points
  udates <- unique(pdt$Date)
  T <- length(udates)
  
  if(isotherm != '') iso.def <- TRUE else iso.def <- FALSE
  
  for(i in 1:T){
    
    time <- udates[i]
    pdt.i <- pdt[which(pdt$Date == time),]
    
    # open day's hycom data
    nc <- open.ncdf(paste(ohc.dir, 'Lyd_', as.Date(time), '.nc', sep=''))
    dat <- get.var.ncdf(nc, 'water_temp')
    
    if(i == 1){
      depth <- get.var.ncdf(nc, 'depth')
      lon <- get.var.ncdf(nc, 'lon')
      lat <- get.var.ncdf(nc, 'lat')
      #f.arr <- array(NA, dim=c(length(lon),length(lat),T))
    }
    
    #extracts depth from tag data for day i
    y <- pdt.i$Depth[!is.na(pdt.i$Depth)] 
    y[y<0] <- 0
    
    #extract temperature from tag data for day i
    x <- pdt.i$MidTemp[!is.na(pdt.i$Depth)]  
    
    # use the which.min
    depIdx = unique(apply(as.data.frame(pdt.i$Depth), 1, FUN=function(x) which.min((x-depth)^2)))
    hycomDep <- depth[depIdx]
    
    # make predictions based on the regression model earlier for the temperature at standard WOA depth levels for low and high temperature at that depth
    fit.low <- locfit(pdt.i$MinTemp ~ pdt.i$Depth)
    fit.high <- locfit(pdt.i$MaxTemp ~ pdt.i$Depth)
    n = length(hycomDep)
    
    pred.low = predict(fit.low, newdata = hycomDep, se = T, get.data = T)
    pred.high = predict(fit.high, newdata = hycomDep, se = T, get.data = T)
    
    # data frame for next step
    df = data.frame(low=pred.low$fit-pred.low$se.fit*sqrt(n)
                    , high=pred.high$fit+pred.high$se.fit*sqrt(n)
                    , depth = hycomDep)
    
    # isotherm is minimum temperature recorded for that time point
    if(iso.def == FALSE) isotherm <- min(df$low, na.rm = T)
    
    # perform tag data integration at limits of model fits
    minT.ohc <- cp * rho * sum(df$low - isotherm, na.rm = T) / 10000
    maxT.ohc <- cp * rho * sum(df$high - isotherm, na.rm = T) / 10000
    midT.ohc <- cp * rho * sum(pdt.i$MidTemp - isotherm, na.rm = T) / 10000
    datProf <- dat[lonIdx,latIdx,]
    #dat.ohc <- cp * rho * sum(datProf[depIdx] - isotherm, na.rm = T) / 10000
    
    # Perform hycom integration
    dat[dat<isotherm] <- NA
    dat <- dat - isotherm
    ohc <- cp * rho * apply(dat[,,depIdx], 1:2, sum, na.rm = T) / 10000 
    ohc[ohc == 0] <- NA
    
    # calc sd of OHC
    # focal calc on mean temp and write to sd var
    r = flip(raster(t(ohc)),2)
    sdx = focal(r, w=matrix(1/9,nrow=3,ncol=3), fun=function(x) sd(x, na.rm = T))
    sdx = t(as.matrix(flip(sdx,2)))
    
    # compare hycom to that day's tag-based ohc
    lik.ohc <- likint2(ohc, sdx, minT.ohc, maxT.ohc)
    
    spot = c(26.3810005, -70.6969986)
    
    pdf('test ohc.pdf',height=12,width=8)
     par(mfrow=c(3,1))
     image.plot(lon-360,lat,ohc)
     plot(countriesLow,add=T)
     points(spot[2], spot[1],pch=16,col='white')
     title('ohc')
     image.plot(lon-360,lat,sdx)
     plot(countriesLow,add=T)
     points(spot[2], spot[1],pch=16,col='white')
     title('sdx')
     image.plot(lon-360,lat,lik.ohc)
     plot(countriesLow,add=T)
     points(spot[2], spot[1],pch=16,col='white')
     title(paste('lik.ohc - ',minT.ohc,', ',maxT.ohc))
    dev.off()

    latIdx <- which.min((spot[1]-lat)^2)
    lonIdx <- which.min((spot[2]-(lon-360))^2)
    plot(depth~c(dat[lonIdx,latIdx,]+isotherm),type='l',ylim=c(1000,0),xlim=c(0,25))
    lines(pdt.i$Depth~pdt.i$MidTemp,col='red')
    lines(df[,3]~df[,2],lty=2,col='red')
    lines(df[,3]~df[,1],lty=2,col='red')
    abline(v=isotherm)
    
    ohc[lonIdx,latIdx]
    
    ####
    dnorm(tag.ohc, mean=ohc, sd=sdx)
    
    lik.pdt[,,b] = likint2(dat.i[,,depIdx[b]], sd.i[,,depIdx[b]], df[b,1], df[b,2])
    
    wlist = array(1e-6, dim=c(dim(ohc)[1], dim(ohc)[2], 2))
    wlist[,,1] = ohc
    wlist[,,2] = sdx
    wlist[is.na(wlist)] = 1e-6
    res<-as.matrix(aaply(wlist, 1:2, .fun = function(x) dnorm(tag.ohc, mean = x[1], sd = .7)))#, lower = minT, upper = maxT , )$value))
    
    if(i == 1){
      # result will be array of likelihood surfaces
      L.ohc <- array(0, dim = c(dim(lik), length(dateVec)))
    }
    
    idx <- which(dateVec == as.Date(time))
    L.ohc[,,idx] = lik
    
    print(paste(time, ' finished.', sep=''))
    
  }
  
  
  crs <- "+proj=longlat +datum=WGS84 +ellps=WGS84"
  list.ohc <- list(x = lon-360, y = lat, z = L.ohc)
  ex <- extent(list.ohc)
  L.ohc <- brick(list.ohc$z, xmn=ex[1], xmx=ex[2], ymn=ex[3], ymx=ex[4], transpose=T, crs)
  L.ohc <- flip(L.ohc, direction = 'y')
  s <- stack(L.ohc)
  return(s)
}
  # make L.pdt match resolution/extent of g
  #row <- dim(g$lon)[1]
  #col <- dim(g$lon)[2]
  #ex <- extent(c(min(g$lon[1,]), max(g$lon[1,]), min(g$lat[,1]), max(g$lat[,1])))
  #crs <- "+proj=longlat +datum=WGS84 +ellps=WGS84"
  #rasMatch <- raster(ex, nrows=row, ncols=col, crs = crs)
  #L.ohc <- spatial_sync_raster(L.ohc, rasMatch)
  
  #if(raster == 'brick'){
   # s <- L.ohc
  #} else if(raster == 'stack'){
  #  s <- stack(L.ohc)
  #} else if(raster == 'array'){
  #  s <- raster::as.array(L.ohc, transpose = T)
  #}
  
#  print(class(L.ohc))
  # return ohc likelihood surfaces
 # return(L.ohc)
  
#}
