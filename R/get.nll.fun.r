#' Negative Log likelihood of parameters
#' 
#' @param parvec vector of length 6 contianing * kernel 1 * kernel 2 * diagonal
#'   of 2x2 matrix
#' @param g grid from \code{\link{setup.grid}}
#' @param L final likelihood (2D)
#'   
#' @return parameter values
#' @export
#' @references Pedersen MW, Patterson TA, Thygesen UH, Madsen H (2011)
#'   Estimating animal behavior and residency from movement data. Oikos
#'   120:1281-1290. doi: 10.1111/j.1600-0706.2011.19044.x
#'   
#' @examples
#' # NOT RUN
#' par0 <- c(log(10), log(10), log(0.5), log(0.5), log(0.95/0.05), log(0.95/0.05))
#' fit <- nlm(get.nll.fun, par0, g, L, dt)
#' D1 <- exp(fit$estimate[1:2])
#' D2 <- exp(fit$estimate[3:4])
#' p <- 1/(1+exp(-fit$estimate[5:6])) 


get.nll.fun <- function(parvec = c(10, 30, 5, 2, .707, .8), g, L){ #c(D1, D2, p)c(10, 30, 5, 2, .707, .8)
  K1 = imager::as.cimg(gausskern(parvec[1], parvec[2], muadv = 0))
  K2 = imager::as.cimg(gausskern(parvec[3], parvec[4], muadv = 0))
  P <- matrix(c(parvec[5], 1-parvec[5], 1-parvec[6], parvec[6]), 2, 2, byrow = TRUE)

  f = hmm.filter_test(g, L, K1, K2, P)
  nllf <- -sum(log(f$psi[f$psi>0]))
  print(paste0("\n HMM -log(L):", nllf))
  #flush.console()
  nllf
  
}

