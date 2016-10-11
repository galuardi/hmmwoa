#' Combine individual source likelihoods
#' 
#' \code{make.L} combines individual likelihoods from various data sources (e.g. SST, OHC) to make overall combined likelihoods for each time point
#' 
#' @param L1 a likelihood array
#' @param L2 a likelihood array
#' @param L3 a likelihood array
#' @param L.locs is array of dim(L1) that contains any known locations during the deployment that will supersede any other likelihood at that time point
#' @param L.mle.res is a coarse resolution array of dim(L1) that speeds up the parameter estimation step later on
#' @param plot is logical indicating whether you want an example plot
#'   
#' @return a list containing: L, the overall likelihood array and L.mle, a more coarse version of L used later for parameter estimation
#' 
#' @examples
#' none

make.L <- function(L1, L2 = NULL, L3 = NULL, known.locs = NULL, L.mle.res, dateVec = NULL, locs.grid = NULL, iniloc = NULL){

  if(!is.null(known.locs)){
    print('known locs are being used')
    # convert input date, lat, lon to likelihood surfaces with dim(L1)
    L.locs <- L1 * 0
    known.locs$date <- as.Date(known.locs$date)
    
    if(is.null(dateVec)){stop('Error: dateVec is null.')}
    if(is.null(locs.grid)){stop('Error: locs.grid is null.')}
    
    # need lat/lon vectors from locs.grid
    lon <- locs.grid$lon[1,]
    lat <- locs.grid$lat[,1]
    
    idx <- which(dateVec %in% known.locs$date)
    for(i in idx){
      known.locs.i <- known.locs[which(known.locs$date %in% dateVec[i]),]
      
      if(length(known.locs.i[,1]) > 1){
        # if multiple known locations are provided for a given day, only the first is used
        known.locs.i <- known.locs.i[1,]
      }
      
      x = which.min((known.locs.i$lon - lon) ^ 2)
      y = which.min((known.locs.i$lat - lat) ^ 2)

      # assign the known location for this day, i, as 1 in likelihood raster
      L.locs[[i]][cellFromXY(L.locs[[idx]], known.locs.i[,c(3,2)])] <- 1
      
    }
    
  }
  
  if(is.null(L2) & is.null(L3)){
    print('entering L1 loop')
    
    L <- L1
    
    # GET RID OF NA VALUES
    L[is.na(L)] <- 0
    
    # ALL CELLS IN A LIKELIHOOD SURFACE == 0?
    naLidx = cellStats(L, sum, na.rm=T) != 0
    
    L <- L[naLidx]
    
  } else if(!is.null(L2) & is.null(L3)){
    print('entering L2 loop')
    
    # then there were 2
    
    # MAKE AN ARRAY OF ZEROS
    L <- L1 * 0
    
    # GET RID OF NA VALUES
    L1[is.na(L1)] <- 0
    L2[is.na(L2)] <- 0
    
    # ALL CELLS IN A LIKELIHOOD SURFACE == 0?
    naL1idx = cellStats(L1, sum, na.rm=T) != 0
    naL2idx = cellStats(L2, sum, na.rm=T) != 0
    
    # WHERE BOTH ARE ZEROS. THESE WILL BE INTERPOLTED IN THE FILTER
    naLidx = naL1idx + naL2idx
    
    # where naLidx==0, both likelihoods are zero
    #       naLidx==1, one has data
    #       naLidx==2, both have data
    idx1 = which(naLidx == 1)
    idx2 = which(naLidx == 2)

    
    # INPUTS GET ADDED/MULTIP TOGETHER FOR FINAL L
    # THESE DONT WORK
    #L[[idx1]] <- L1[[idx1]] + L2[[idx1]] 
    #L[[idx1]] <- overlay(L1[[idx1]], L2[[idx1]], fun=function(a,b){a + b})
    
    # BUT THIS DOES
    for(ii in idx1){
      L[[ii]] = L1[[ii]] + L2[[ii]] # when only 1 has data
    }
    
    for(ii in idx2){
      L[[ii]] = L1[[ii]] * L2[[ii]] # when both have data
    }

    

  } else if(!is.null(L2) & !is.null(L3)){
    print('entering L2 and L3 loop')
    
    # then there were 3
    
    # MAKE AN ARRAY OF ZEROS
    L <- L1 * 0
    
    # GET RID OF NA VALUES
    L1[is.na(L1)] <- 0
    L2[is.na(L2)] <- 0
    L3[is.na(L3)] <- 0
    
    # ALL CELLS IN A LIKELIHOOD SURFACE == 0?
    naL1idx = cellStats(L1, sum, na.rm=T) != 0
    naL2idx = cellStats(L2, sum, na.rm=T) != 0
    naL3idx = cellStats(L3, sum, na.rm=T) != 0
    
    # WHERE ALL ARE ZEROS. THESE WILL BE INTERPOLTED IN THE FILTER
    naLidx = naL1idx + naL2idx + naL3idx
    
    # where naLidx==0, both likelihoods are zero
    #       naLidx==1, one has data
    #       naLidx==2, both have data
    idx1 = which(naLidx == 1)
    idx2 = which(naLidx == 2)
    idx3 = which(naLidx == 3)
    
    # COMBINING LIKELIHOODS
    for(ii in idx1){
      L[[ii]] = L1[[ii]] + L2[[ii]] + L3[[ii]]# when only 1 has data
    }
    
    for(ii in idx3){
      L[[ii]] = L1[[ii]] * L2[[ii]] * L3[[ii]] # when all have data
    }
    
    for(ii in idx2){
      if(naL1idx[ii] & naL2idx[ii]){
        L[[ii]] <- L1[[ii]] * L2[[ii]]
      } else if(naL1idx[ii] & naL3idx[ii]){
        L[[ii]] <- L1[[ii]] * L3[[ii]]
      } else if(naL2idx[ii] & naL3idx[ii]){
        L[[ii]] <- L2[[ii]] * L3[[ii]]
      }
    }
    
  }
  
  if(!is.null(known.locs)){
    print('entering known.locs loop at the end')
    print(known.locs)
    print(!is.null('known.locs'))
     for(bb in idx){
       L[[bb]] <- L.locs[[bb]]
     }
  }
  
  if(!is.null(iniloc)){
    print('entering iniloc loop at the end')
    
    if(!exists('L.locs')){
      L.locs <- L1 * 0
      idx <- c(1, length(dateVec))
    } else{
      idx <- c(1, idx, length(dateVec))
    }
    
    # need lat/lon vectors from locs.grid
    lon <- locs.grid$lon[1,]
    lat <- locs.grid$lat[,1]
    
    # tag location
    x = which.min((iniloc$lon[1] - lon) ^ 2)
    y = which.min((iniloc$lat[1] - lat) ^ 2)
    print(paste('tag',x,y))
    # assign the known location for this day, i, as 1 in likelihood raster
    L.locs[[1]][cellFromXY(L.locs[[1]], iniloc[1,c(5,4)])] <- 1
    
    # pop up location
    x = which.min((iniloc$lon[2] - lon) ^ 2)
    y = which.min((iniloc$lat[2] - lat) ^ 2)
    print(paste('pop',x,y))
    
    # assign the known location for this day, i, as 1 in likelihood raster
    L.locs[[length(dateVec)]][cellFromXY(L.locs[[length(dateVec)]], iniloc[2,c(5,4)])] <- 1
    
    # add known to L
    for(bb in idx){
      L[[bb]] <- L.locs[[bb]]
    }
  }
  
  # CREATE A MORE COARSE RASTER FOR PARAMETER ESTIMATION LATER
  L.mle <- raster::resample(L, L.mle.res)
  
  #----------------------------------------------------------------------------------#
  # MAKE ALL NA'S VERY TINY FOR THE CONVOLUTION
  # the previous steps may have taken care of this...
  #----------------------------------------------------------------------------------#
  L[L == 0] <- 1e-15
  L[is.na(L)] <- 1e-15
  L.mle[L.mle == 0] <- 1e-15
  L.mle[is.na(L.mle)] <- 1e-15
  
  # MAKE BOTH RASTERS (COARSE AND FINE RES L's) INTO AN ARRAY
  L <- aperm(raster::as.array(raster::flip(L, direction = 'y')), c(3, 2, 1))
  L.mle <- aperm(raster::as.array(raster::flip(L.mle, direction = 'y')), c(3, 2, 1))
  
  return(list(L = L, L.mle = L.mle))
}

## END
