# Next up.... 
#' Title
#'
#' @param f 
#' @param K1 
#' @param K2 
#' @param P 
#' @param plot 
#'
#' @return
#' @export
#'
#' @examples
hmm.smoother_test <- function(f,K1,K2,P,plot=TRUE){
  ## Smoothing the filtered estimates
  ## The equations for smoothing are presented in Pedersen et al. 2011, Oikos, Appendix
  T <- dim(f$phi)[2]
  row <- dim(f$phi)[3]
  col <- dim(f$phi)[4]
  K1 <- imager::as.cimg(K1)
  K2 <- imager::as.cimg(K2)
  
#   smooth <- array(0,dim=dim(f$phi))
#   smooth[,T,,] <- f$phi[,T,,]  # end point
  smooth <- f$phi  #default; fill in as the prediction step.
  # smooth[,T,,] <- L[T,,]/sum(L[T,,], na.rm = T) # last position is known; use the likelihood matrix, same as the filter step for starting point.. 
  # smooth[,1,,] <- f$phi[,1,,]  # starting point
  for(t in T:2){
    RAT <- smooth[,t,,] / (f$pred[,t,,] + 1e-15)
    #     Rp1 <- as.vector(K1 %*% as.vector(RAT[1,,]))
    #     Rp2 <- as.vector(K2 %*% as.vector(RAT[2,,]))
    
    
    p1 = imager::as.cimg(t(RAT[1,,]))
    Rp1 <- imager::convolve(p1, K1)
    p2 = imager::as.cimg(t(RAT[2,,]))
    Rp2 <- imager::convolve(p2, K2)
    
    Rp1 = t(as.matrix(Rp1))
    Rp2 = t(as.matrix(Rp2))
    
    if(plot){
      par(mfrow = c(1,2))
      image.plot(Rp1)
      #plot(countriesLow,add=T)
      image.plot(Rp2)
      #plot(countriesLow,add=T)
    }
    
    post1 <- matrix(P[1,1] * Rp1 + P[1,2] * Rp2, row, col)
    post2 <- matrix(P[2,1] * Rp1 + P[2,2] * Rp2, row, col)
    if(T == t){
      post1 <- f$phi[1,t,,]
      post2 <- f$phi[2,t,,]
      fac <- sum(as.vector(post1)) + sum(as.vector(post2))
      smooth[1,t,,] <- post1 / fac
      smooth[2,t,,] <- post2 / fac 
      post1 <- post1 * f$phi[1,t-1,,]
      post2 <- post2 * f$phi[2,t-1,,]
      fac <- sum(as.vector(post1)) + sum(as.vector(post2))
      smooth[1,t-1,,] <- post1 / fac
      smooth[2,t-1,,] <- post2 / fac
    }else{
      post1 <- post1 * f$phi[1,t-1,,]
      post2 <- post2 * f$phi[2,t-1,,]
      fac <- sum(as.vector(post1)) + sum(as.vector(post2))
      smooth[1,t-1,,] <- post1 / fac
      smooth[2,t-1,,] <- post2 / fac
    }
  }
  
  smooth
  
}

