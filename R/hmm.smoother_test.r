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
hmm.smoother_test <- function(f, K1, K2, P){
  ## Smoothing the filtered estimates
  ## The equations for smoothing are presented in Pedersen et al. 2011, Oikos, Appendix
  T <- dim(f$phi)[2]
  row <- dim(f$phi)[3]
  col <- dim(f$phi)[4]
  
  # convert movement kernel from matrix to cimg for convolution
  K1 <- imager::as.cimg(K1)
  K2 <- imager::as.cimg(K2)
  
  smooth <- array(0, dim = dim(f$phi))
  smooth[,T,,] <- f$phi[,T,,]
  
  #smooth <- f$phi  #default; fill in as the prediction step.

  for(t in T:2){
    RAT <- smooth[,t,,] / (f$pred[,t,,] + 1e-15)
    
    # convolve today's smoother prediction with movement kernel
    p1 = imager::as.cimg(t(RAT[1,,]))
    Rp1 <- imager::convolve(p1, K1)
    p2 = imager::as.cimg(t(RAT[2,,]))
    Rp2 <- imager::convolve(p2, K2)
    Rp1 = t(as.matrix(Rp1))
    Rp2 = t(as.matrix(Rp2))

    post1 <- matrix(P[1,1] * Rp1 + P[1,2] * Rp2, row, col)
    post2 <- matrix(P[2,1] * Rp1 + P[2,2] * Rp2, row, col)
    
    #if(T == t){
    #  post1 <- f$phi[1,t,,]
    #  post2 <- f$phi[2,t,,]
    #  fac <- sum(as.vector(post1)) + sum(as.vector(post2))
    #  smooth[1,t,,] <- post1 / fac
    #  smooth[2,t,,] <- post2 / fac 
    #  post1.try <- post1 * f$phi[1,t-1,,]
    #  post2.try <- post2 * f$phi[2,t-1,,]
    #  fac.try <- sum(as.vector(post1.try)) + sum(as.vector(post2.try))
    #  smooth[1,t-1,,] <- post1.try / fac.try
    #  smooth[2,t-1,,] <- post2.try / fac.try
    #}else{
      post1 <- post1 * f$phi[1,t-1,,]
      post2 <- post2 * f$phi[2,t-1,,]
      fac <- sum(as.vector(post1)) + sum(as.vector(post2))
      smooth[1,t-1,,] <- post1 / fac
      smooth[2,t-1,,] <- post2 / fac
    #}
  }
  
  smooth
  
}

