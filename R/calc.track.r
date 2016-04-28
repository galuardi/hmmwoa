#' HMM Smoother
#' 
#' @param distr 
#' @param g 
#'
#' @return
#' @export
#'
#' @examples
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