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

make.L <- function(L1, L2 = NULL, L3 = NULL, L.mle.res, plot = TRUE){

  if(is.null(L2) & is.null(L3)){
    L <- L1
    
    # GET RID OF NA VALUES
    L[is.na(L)] <- 0
    
    # ALL CELLS IN A LIKELIHOOD SURFACE == 0?
    naLidx = cellStats(L, sum, na.rm=T) != 0
    
    L <- L[naLidx]
    
  } else if(!is.null(L2) & is.null(L3)){
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
    
    # WHERE BOTH ARE ZEROS. THESE WILL BE INTERPOLTED IN THE FILTER
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
  
  # CREATE A MORE COARSE RASTER FOR PARAMETER ESTIMATION LATER
  L.mle <- raster::resample(L, L.mle.res$L.locs)

  if(plot){
    # CHECK THAT IT WORKED OK
    require(fields)
    image.plot(L[1,,])
    plot(countriesLow,add=T)
    
  }
  
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
