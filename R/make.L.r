#' Read and format tag data
#' 
#' \code{make.L} reads and formats tag data output from Wildlife Computers Data Portal
#' 
#' @param 
#'   
#' @return a list containing: 
#' 
#' @export
#' 
#' @examples
#' none


make.L <- function(L.res, idx = c(1,2,3,4), plot = TRUE){
  # where 1==pdt, 2==ohc, 3==locs, 4==sst
  
  # MAKE AN ARRAY OF ZEROS
  Lmat = L.res[[1]][[1]] * 0
  
  # RASTER TO ARRAY
  L.res <- lapply(L.res[[1]], FUN=raster::as.array)
  
  # GET RID OF NA VALUES - not working
  L.res[is.na(L.res)] = 0
  
  napdtidx = apply(L.res$L.pdt, 3, sum, na.rm=T) != 0
  
  # need to sum rows and cols at each layer as in apply(x, 3, sum) but for a raster
  pdtidx <- sum(L.pdt, na.rm=T) != 0
  
  if(length(idx) == 1){
    
    # only one input then there's nothing to do here
    L.res[[1]][idx]
    
  } else if(length(idx) == 2){
    
    # two inputs get added/multiplied together
    Lmat[,,idx1] = L.locs[,,idx1] + L.sst[,,idx1] # when only 1 has data
    Lmat[,,idx2] = L.sst[,,idx2] * L.locs[,,idx2] # when both have data
    
  } else if(length(idx) == 3){
    
    Lmat[,,idx1] = L.locs[,,idx1] + L.sst[,,idx1] + L.ohc[,,idx1] # when only 1 has data
    Lmat[,,idx3] = L.locs[,,idx3] * L.sst[,,idx3] * L.ohc[,,idx3] # when all have data
    
    # USE THE INDICES TO POPULATE L
    for(b in which(idx2)){
      if(nasstidx[b] & nalocidx[b]){
        Lmat[,,b] = L.sst[,,b] * L.locs[,,b]
      } else if(nasstidx[b] & naohcidx[b]){
        Lmat[,,b] = L.sst[,,b] * L.ohc[,,b]
      } else if(nalocidx[b] & naohcidx[b]){
        Lmat[,,b] = L.locs[,,b] * L.ohc[,,b]
      }
    }
    
  } else if(length(idx) == 4){
    
    stop('This many input indices is not yet functional. Make sure you are not trying to use both pdt and ohc inputs.')
    
  }
  
  # MAKE A LIST OF LIKELIHOOD TO SET RASTER EXTENT
  # define projection
  crs <- "+proj=longlat +datum=WGS84 +ellps=WGS84"
  list.pdt <- list(x = lon, y = lat, z = Lmat)
  ex <- raster::extent(list.pdt)
  
  # MAKE A RASTER OUT OF IT
  T <- dim(Lmat)[3]
  for(i in 1:T){
    L.i <- raster::raster(Lmat[,,i], xmn=ex[1], xmx=ex[2], ymn=ex[3], ymx=ex[4], crs)
    if(i==1) L <- L.i else L <- stack(L, L.i)
  }
  
  # CREATE A MORE COARSE RASTER FOR PARAMETER ESTIMATION LATER
  L.mle <- raster::resample(L, L.mle.res)
  
  # MAKE BOTH RASTERS (COARSE AND FINE RES L's) INTO AN ARRAY
  #L <- aperm(raster::as.array(raster::flip(L, direction = 'y')), c(3, 2, 1))
  #L.mle <- aperm(raster::as.array(raster::flip(L.mle, direction = 'y')), c(3, 2, 1))
  
  if(plot){
    # CHECK THAT IT WORKED OK
    plot(L[[1]])
    plot(countriesLow,add=T)
    
  }
  
  #----------------------------------------------------------------------------------#
  # MAKE ALL NA'S VERY TINY FOR THE CONVOLUTION
  # the previous steps may have taken care of this...
  #----------------------------------------------------------------------------------#
  L[L == 0] = 1e-15
  L[is.na(L)] = 1e-15
  L.mle[L.mle == 0] = 1e-15
  L.mle[is.na(L.mle)] = 1e-15
  
  
  
}

## END

# are all cells in a given likelihood surface == 0?
#na1idx = lapply(L.res[[1]], FUN = (apply(L.res[[1]][[1]], 3, sum, na.rm=T) != 0))# does sum of likelihood surface
#na2idx = apply(L.2, 3, sum, na.rm=T) != 0
#na3idx = apply(L.3, 3, sum, na.rm=T) != 0
#na4idx = apply(L.4, 3, sum, na.rm=T) != 0


# HERE, WE CHOOSE WHICH L's TO USE
# INDICATES WHICH L LAYERS, IF ANY, ARE ALL ZEROS FOR EACH DAY


# WHERE BOTH ARE ZEROS. THESE WILL BE INTERPOLTED IN THE FILTER
#naLidx = nalocidx + nasstidx + naohcidx

# MAKE AN ARRAY OF ZEROS
#Lmat = L.pdt * 0
# where naLidx==0, both likelihoods are zero
#       naLidx==1, one has data
#       naLidx==2, both have data
#idx1 = naLidx == 1
#idx2 = naLidx == 2
#idx3 = naLidx == 3

#Lmat[,,idx1] = L.locs[,,idx1] + L.sst[,,idx1] + L.ohc[,,idx1] # when only 1 has data
#Lmat[,,idx2] = L.sst[,,idx2] * L.locs[,,idx2] # when both have data
#Lmat[,,idx3] = L.locs[,,idx3] * L.sst[,,idx3] * L.ohc[,,idx3] # when all have data

# USE THE INDICES TO POPULATE L
#for(b in which(idx2)){
#  if(nasstidx[b] & nalocidx[b]){
#    Lmat[,,b] = L.sst[,,b] * L.locs[,,b]
#  } else if(nasstidx[b] & naohcidx[b]){
#    Lmat[,,b] = L.sst[,,b] * L.ohc[,,b]
#  } else if(nalocidx[b] & naohcidx[b]){
#    Lmat[,,b] = L.locs[,,b] * L.ohc[,,b]
#  }
#}

# MAKE A LIST OF LIKELIHOOD TO SET RASTER EXTENT
# define projection
#crs <- "+proj=longlat +datum=WGS84 +ellps=WGS84"
#list.pdt <- list(x = lon, y = lat, z = Lmat)
#ex <- raster::extent(list.pdt)

# MAKE A RASTER OUT OF IT
#T <- dim(Lmat)[3]
#for(i in 1:T){
#  L.i <- raster::raster(Lmat[,,i], xmn=ex[1], xmx=ex[2], ymn=ex[3], ymx=ex[4], crs)
#  if(i==1) L <- L.i else L <- stack(L, L.i)
#}

# CREATE A MORE COARSE RASTER FOR PARAMETER ESTIMATION LATER
#L.mle <- raster::resample(L, L.mle.res)

# MAKE BOTH RASTERS (COARSE AND FINE RES L's) INTO AN ARRAY
#L <- aperm(raster::as.array(raster::flip(L, direction = 'y')), c(3, 2, 1))
#L.mle <- aperm(raster::as.array(raster::flip(L.mle, direction = 'y')), c(3, 2, 1))

#if(plot){
#  # CHECK THAT IT WORKED OK
#  lon <- g$lon[1,]
#  lat <- g$lat[,1]
#  fields::image.plot(lon, rev(lat), L[1,,])
#  par(mfrow=c(2,1))
#  fields::image.plot(lon, lat, L2[1,,])
#  image.plot(L1[1,,])
#  plot(countriesLow,add=T)
#  
#}

#----------------------------------------------------------------------------------#
# MAKE ALL NA'S VERY TINY FOR THE CONVOLUTION
# the previous steps may have taken care of this...
#----------------------------------------------------------------------------------#
#L[L == 0] = 1e-15
#L[is.na(L)] = 1e-15
#L.mle[L.mle == 0] = 1e-15
#L.mle[is.na(L.mle)] = 1e-15


#*** end function

# INDEX WHERE LIKELIHOODS ARE ZEROS.. FOR EACH L COMPONENT
#L.1 = raster::as.array(L.res[[1]][[1]])
#L.2 = raster::as.array(L.res[[1]][2])
#L.3 = raster::as.array(L.res[[1]][3])
#L.4 = raster::as.array(L.res[[1]][4])

# turn NA to 0
#L.1[is.na(L.1)] = 0
#L.2[is.na(L.2)] = 0
#L.3[is.na(L.3)] = 0 
#L.4[is.na(L.4)] = 0

cellStats yields sum of each layer in raster stack
