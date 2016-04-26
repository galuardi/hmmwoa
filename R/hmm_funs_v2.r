# New HMM and SMoother functions
# add documentation...

#### Function based on Pedersen 2011, originally using sparse matrix multiplication.. editing for gaussian kernel convolution
# watch out for matrix dimensionality. even though dimensions lined up, data did not.. 
hmm.filter2 <- function(g,L,K1,K2,P){
  
  ## Filter data to estimate locations and behaviour
  
  T <- dim(L)[1] # dimension of time 
  row <- dim(g$lon)[1] # nrows
  col <- dim(g$lon)[2] # ncols
  m <- 2 # Number of behavioural states
  
  pred <- array(0,dim=c(m,T,col,row)) # empty array for prediction step. ordering in col before row emulates lon before lat
  phi  <- array(0,dim=c(m,T,col,row)) # posterior (final) step array
  # Start in resident state at the known initial location
  phi[2,1,,]  <- L[1,,] # first position is known
  pred[2,1,,] <- L[1,,] # first position is known
  psi <- rep(0,T-1) # sum of the probability of both states at each step
  # Start filter iterations
  for(t in 2:T){
    # replace this part with older workflow using a gaussian kernel.. 
    # p1 <- as.vector(phi[1,t-1,,])
    # p2 <- as.vector(phi[2,t-1,,])
    # q1 <- as.vector(p1%*%K1)
    # q2 <- as.vector(p2%*%K2)
    
    p1 = as.cimg(t(phi[1,t-1,,]))
    p2 = as.cimg(t(phi[2,t-1,,]))
    q1 = convolve(p1, K1)
    q2 = convolve(p2, K2)
    
    # q1 = arot(t(as.matrix(q1)),3)
    # q2 = arot(t(as.matrix(q2)),3)
    q1 = t(as.matrix(q1))
    q2 = t(as.matrix(q2))
    
    # 	par(mfrow=c(1,2))
    # 	image(q1)
    # 	image(q2)
    
    # pred[1,t,,] <- matrix(P[1,1]*q1+P[2,1]*q2,row,col)
    # pred[2,t,,] <- matrix(P[1,2]*q1+P[2,2]*q2,row,col)
    
    # multiply by transition probability 
    pred[1,t,,] <- P[1,1]*q1+P[2,1]*q2
    pred[2,t,,] <- P[1,2]*q1+P[2,2]*q2
    
    sumL = sum(L[t,,])  
    if(sumL > 0){
      post1 <- pred[1,t,,]*L[t,,]
      post2 <- pred[2,t,,]*L[t,,]
    }else{
      post1 <- pred[1,t,,]
      post2 <- pred[2,t,,]
    }
    
    psi[t-1] <- sum(as.vector(post1), na.rm=T) + sum(as.vector(post2), na.rm=T)
    
    # remove NaNs... 
    # normalise (divide here by sum, not max)
    # 	post1 <- normalise(post1)
    # 	post2 <- normalise(post2)
    # 	post1[is.nan(post1)] = 0
    # 	post2[is.nan(post2)] = 0
    
    phi[1,t,,] <- post1/(psi[t-1]+1e-15)
    phi[2,t,,] <- post2/(psi[t-1]+1e-15)
  }
  list(phi=phi,pred=pred,psi=psi)
}

# Next up.... 
hmm.smoother2 <- function(f,K1,K2,P,plot=TRUE){
  ## Smoothing the filtered estimates
  ## The equations for smoothing are presented in Pedersen et al. 2011, Oikos, Appendix
  T <- dim(f$phi)[2]
  row <- dim(f$phi)[3]
  col <- dim(f$phi)[4]
  
  smooth <- array(0,dim=dim(f$phi))
  smooth[,T,,] <- f$phi[,T,,]
  for(t in T:2){
    RAT <- smooth[,t,,]/(f$pred[,t,,]+1e-15)
    #     Rp1 <- as.vector(K1 %*% as.vector(RAT[1,,]))
    #     Rp2 <- as.vector(K2 %*% as.vector(RAT[2,,]))
    
    
    p1 = as.cimg(t(RAT[1,,]))
    Rp1 <- convolve(p1, K1)
    p2 = as.cimg(t(RAT[2,,]))
    Rp2 <- convolve(p2, K2)
    
    Rp1 = t(as.matrix(Rp1))
    Rp2 = t(as.matrix(Rp2))
    
    if(plot){
      par(mfrow=c(1,2))
      image.plot(Rp1)
      #plot(countriesLow,add=T)
      image.plot(Rp2)
      #plot(countriesLow,add=T)
    }
    
    post1 <- matrix(P[1,1]*Rp1 + P[1,2]*Rp2,row,col)
    post2 <- matrix(P[2,1]*Rp1 + P[2,2]*Rp2,row,col)
    post1 <- post1 * f$phi[1,t-1,,]
    post2 <- post2 * f$phi[2,t-1,,]
    fac <- sum(as.vector(post1)) + sum(as.vector(post2))
    smooth[1,t-1,,] <- post1/fac
    smooth[2,t-1,,] <- post2/fac
  }
  smooth
}

# unmodified from sphmm
calc.track <- function(distr,g){
  ## Calculate track from probability distribution of location
  
  T <- dim(distr)[2]
  # Track calculated from mean
  meanlat <- apply(apply(distr,c(2,3),sum)*repmat(t(as.matrix(g$lat[,1])),T,1),1,sum)
  meanlon <- apply(apply(distr,c(2,4),sum)*repmat(t(as.matrix(g$lon[1,])),T,1),1,sum)
  # Track calculated from mode
  row <- dim(g$lon)[1]
  col <- dim(g$lon)[2]
  modelat <- rep(0,T)
  modelon <- rep(0,T)
  for(t in 1:T){
    asd <- apply(distr[,t,,],c(2,3),sum)
    ind <- which.max(asd)
    x <- ceiling(ind/col)
    y <- ind %% row
    modelat[t] <- g$lat[y,x]
    modelon[t] <- g$lon[y,x]
  }
  # Track calculated with Viterbi
  # --- not included in this script
  list(meanlon=meanlon,meanlat=meanlat,modelon=modelon,modelat=modelat)
}
