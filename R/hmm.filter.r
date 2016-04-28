#' HMM filter and smoother functions
#'
#' @param g grid from \code{\link{setup.grid}}
#' @param L final likelihood (2D)
#' @param K1 first movement (diffusion) kernel see \code{\link{gausskern}}
#' @param K2 second movement (diffusion) kernel see \code{\link{gausskern}}
#' @param P 2x2 probability matrix for tranisitons between states (K1 and K2)
#'
#' @return a list: list(phi = phi, pred = pred, psi = psi) where 
#' \item phi is the probability for each state at each trime step 
#' \item pred is ....
#' \item psi is.... 
#' @export
#'
#' @examples
#' P <- matrix(c(p[1], 1-p[1], 1-p[2], p[2]), 2, 2, byrow=TRUE)
hmm.filter <- function(g,L,K1,K2,P){
  
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
  list(phi = phi, pred = pred, psi = psi)
}


