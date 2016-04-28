#' Negative Log likelihood of parameters
#'
#' @param parvec vector of length 6 contianing
#' * kernel 1
#' * kernel 2
#' * diagonal of 2x2 matrix 
#' @param g grid from \code{\link{setup.grid}}
#' @param L final likelihood (2D)
#' @param dt time step
#'
#' @return parameter values
#' @export
#' @references Pedersen MW, Patterson TA, Thygesen UH, Madsen H (2011) Estimating animal behavior and residency from movement data. Oikos 120:1281-1290. doi: 10.1111/j.1600-0706.2011.19044.x
#'
#' @examples
#' # NOT RUN
#' par0 <- c(log(10), log(10), log(0.5), log(0.5), log(0.95/0.05), log(0.95/0.05))
#' fit <- nlm(get.nll.fun, par0, g, L, dt)
#' D1 <- exp(fit$estimate[1:2])
#' D2 <- exp(fit$estimate[3:4])
#' p <- 1/(1+exp(-fit$estimate[5:6])) 


get.nll.fun <- function(parvec=c(D1, D2, p), g, L){
  K1 = as.cimg(gausskern(parvec[1], parvec[2], muadv = 0))
  K2 = as.cimg(gausskern(parvec[3], parvec[4], muadv = 0))
  P <- matrix(c(parvec[5], 1-parvec[5], 1-parvec[6], parvec[6]), 2, 2, byrow = TRUE)

  # make all NA's very tiny for the convolution
  # the previous steps may have taken care of this...
#   L[L==0] = 1e-15
#   L[is.na(L)] = 1e-15
  
  # filter - moved function to sphmmfuns_hmm
  f = hmm.filter2(g,L,K1,K2,P)
  nllf <- -sum(log(f$psi))
  ##print(nllf)
  cat("\r HMM -log(L):",nllf)
  #flush.console()
  nllf
}




# 
# neg.log.lik.fun <- function(parvec,g,L,dt){
#   ## Calculating the likelihood of the parameters
#   
#   ## Transform parameters
#   D1 <- exp(parvec[1:2])
#   D2 <- exp(parvec[3:4])
#   p1 <- exp(parvec[5])/(1+exp(parvec[5])) # inv logit
#   p2 <- exp(parvec[6])/(1+exp(parvec[6])) # inv logit
#   P <- matrix(c(p1,1-p1,1-p2,p2),2,2,byrow=TRUE)
#   ## Calculate transition matrices for the two behaviours
#   G1 <- make.kern(D1,g,dt)
#   K1 <- uniformization(G1,dt)
#   G2 <- make.kern(D2,g,dt)
#   K2 <- uniformization(G2,dt)
#   ## Evaluate filter, returns likelihood value
#   f <- hmm.filter(g,L,K1,K2,P)
#   nllf <- -sum(log(f$psi))
#   ##print(nllf)
#   cat("\r HMM -log(L):",nllf); flush.console()
#   nllf
# }
# 
# guess <- c(log(10),log(10),log(0.5),log(0.5),log(0.95/0.05),log(0.95/0.05))
# fit <- nlm(neg.log.lik.fun,guess,g,L,dt)
# D1 <- exp(fit$estimate[1:2])
# D2 <- exp(fit$estimate[3:4])
# p <- 1/(1+exp(-fit$estimate[5:6]))
# 
