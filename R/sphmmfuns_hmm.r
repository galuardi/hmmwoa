# Spatial HMM functions
# 18.10.2011
#  Author: Martin W. Pedersen (map@aqua.dtu.dk)
#
# Demonstration of the spatial HMM method for estimating location and behaviour
# simultaneously from PSAT data, i.e. data transmited from a satellite tag with
# two components: daylight inferred longitude coordinate, and sea-surface temperature
# logged at (in this example case) 12 hour intervals.
#
# The method is explained in Pedersen et al. (2011) Oikos, and is essentially a
# hidden Markov model with a three-dimensional state vector: two spatial coordinates
# and the behavioural coordinate. The spatial components are discretized i.e.
# partitioned into a finite number of states.
#
# The code can estimate model parameter values by setting do.fit = TRUE. For the
# default grid (50x49 cells) this procedure is relatively slow. For faster results
# the grid can be coarsened, for better results the grid can be refined.
#================================================================================

## I modified Martin's original setup.grid function (commented out here)
## and have the newer version stored in and loaded from misc_funs.r

#setup.grid <- function(lsst,ngrid){
  ## Setup the discrete spatial grid for the HMM
  
  #T <- length(lsst$lon)
  
  # Find longitude extents
 # il <- min(lsst$lon)
#  al <- max(lsst$lon)
#  lx <- 0.1*(al-il)
#  lonl <- il - lx
#  lonu <- al + lx
  
  # Find latitude extents
#  latvec <- seq(0,90)
#  lats <- rep(0,T)
#  for(t in 1:T){
#    time <- date2time(lsst$date[t])
#    #time <- as.numeric(strftime(lsst$date[t],format='%j'))
##    ssts <- sstdb(time,lsst$lon[t],latvec)
#    lats[t] <- latvec[sum(lsst$sst[t]<ssts)]
#  }
#  lx <- 0.1*(max(lats)-min(lats))
#  latl <- min(lats) - lx
#  latu <- max(lats) + lx
#  
#  # Create grid
#  lo <- seq(lonl,lonu,length.out=ngrid[1])
#  la <- seq(latl,latu,length.out=ngrid[2])
#  g <- meshgrid(lo,la)
#  dlo <- lo[2]-lo[1]
#  dla <- la[2]-la[1]
#  
#  list(lon=g$X,lat=g$Y,dlo=dlo,dla=dla)
#}
#
data.lik <- function(lsst,iniloc,g){
  ## Calculate the "data" likelihood, i.e. the likelihood field for each observation
  
  T <- length(lsst$lon)
  row <- dim(g$lon)[1]
  col <- dim(g$lon)[2]
  
  L <- array(0,dim=c(T,row,col))
  
  # Initial location is known
  ilo <- which.min(abs(g$lon[1,]-iniloc$lon[1]))
  ila <- which.min(abs(g$lat[,1]-iniloc$lat[1]))
  L[1,ila,ilo] <- 1
  
  # Calculate data likelihood
  # SD for light-based longitude from Musyl et al. (2001)
  sl.sd <- 35/111 # Converting from kms to degrees
  # SD for SST from Pedersen et al. (2011) Oikos
  sst.sd <- 0.71
  
  for(t in 2:(T-1)){
    time <- date2time(lsst$date[t])
    #time <- as.numeric(strftime(lsst$date[t],format='%j'))
    sst <- sstdb(time,g$lon,g$lat)
    Lsst <- dnorm(sst,lsst$sst[t],sst.sd)  # SST data
    Llon <- dnorm(g$lon,lsst$lon[t],sl.sd) # Longitude data
    L[t,,] <- Lsst*Llon
  }
  # End location is known
  elo <- which.min(abs(g$lon[1,]-iniloc$lon[2]))
  ela <- which.min(abs(g$lat[,1]-iniloc$lat[2]))
  L[T,ela,elo] <- 1
  
  L
}

make.kern <- function(D,g,dt){
  ## Calculate the infinitisemal generator matrices of the process
  
  row <- dim(g$lon)[1]
  col <- dim(g$lon)[2]
  
  NN <- row*col # Size of generator
  H <- matrix(0,NN,NN)
  Jx <- D[1]
  Jy <- D[2]
  
  # Jump left and right (x-direction)
  j = (row+1):NN
  print(paste('1'))
  for(i in 1:(NN-row)){
    H[i,j[i]] <- Jx # Jump right
    H[j[i],i] <- Jx # Jump left
  }
  print(paste('2'))
  # Jump up and down (y-direction)
  j = 2:NN;
  for(i in 1:(NN-1)){
    H[i,j[i]] <- Jy # Jump down
    H[j[i],i] <- Jy # Jump up
  }
  print(paste('3'))
  # Take care of top of domain
  top <- seq(row+1,NN,by=row)
  H[top,top-1] <- 0
  # Take care of bottom of domain
  bot <- seq(row,NN-row,by=row)
  H[bot,bot+1] <- 0
  # Sum rows
  H <- H + diag(-apply(H,1,sum))
  print(paste('4'))
  Matrix(H,sparse=TRUE,forceCheck=TRUE)
}

uniformization <- function(A,dt){
  ## Uniformization is an efficient way to compute the matrix exponential for a generator matrix
  ## See Grassmann 1977 - Transient solutions in Markovian queueing systems
  
  N <- dim(A)[1]
  # Find the numerical largest rate
  F <- -min(diag(A))
  # Calculate number of iterations based expression in Grassmann 1977, eq 10
  m <- ceiling(F*dt + 4*sqrt(F*dt) + 5)
  # Insert warning if m>140 ??
  #I <- spMatrix(N,N,i=1:N,j=1:N,rep(1,N))
  I <- Diagonal(N)
  #I <- diag(N)
  P <- A/F + I # Create sub-stochastic matrix (eq 8)
  
  S <- I
  pt <- I
  FPdt <- F*P*dt
  for(i in 1:m){
    S <- S %*% FPdt
    fact <- exp(lgamma(i+1))
    pt <- pt + S/fact
  }
  # Multiply by the constant
  pt <- pt * exp(-F*dt)
  pt
}

hmm.filter <- function(g,L,K1,K2,P){
  ## Filter data to estimate locations and behaviour
  
  T <- dim(L)[1]
  row <- dim(g$lon)[1]
  col <- dim(g$lon)[2]
  m <- 2 # Number of behavioural states
  
  pred <- array(0,dim=c(m,T,row,col))
  phi  <- array(0,dim=c(m,T,row,col))
  # Start in resident state at the known initial location
  phi[2,1,,]  <- L[1,,]
  pred[2,1,,] <- L[1,,]
  psi <- rep(0,T-1)
  # Start filter interations
  for(t in 2:T){
    p1 <- as.vector(phi[1,t-1,,])
    p2 <- as.vector(phi[2,t-1,,])
    q1 <- as.vector(p1%*%K1)
    q2 <- as.vector(p2%*%K2)
    
    pred[1,t,,] <- matrix(P[1,1]*q1+P[2,1]*q2,row,col)
    pred[2,t,,] <- matrix(P[1,2]*q1+P[2,2]*q2,row,col)
    
    post1 <- pred[1,t,,]*L[t,,]
    post2 <- pred[2,t,,]*L[t,,]
    psi[t-1] <- sum(as.vector(post1)) + sum(as.vector(post2))
    phi[1,t,,] <- post1/(psi[t-1]+1e-15)
    phi[2,t,,] <- post2/(psi[t-1]+1e-15)
  }
  list(phi=phi,pred=pred,psi=psi)
}

neg.log.lik.fun <- function(parvec,g,L,dt){
  ## Calculating the likelihood of the parameters
  
  ## Transform parameters
  D1 <- exp(parvec[1:2])
  D2 <- exp(parvec[3:4])
  p1 <- exp(parvec[5])/(1+exp(parvec[5])) # inv logit
  p2 <- exp(parvec[6])/(1+exp(parvec[6])) # inv logit
  P <- matrix(c(p1,1-p1,1-p2,p2),2,2,byrow=TRUE)
  ## Calculate transition matrices for the two behaviours
  G1 <- make.kern(D1,g,dt)
  K1 <- uniformization(G1,dt)
  G2 <- make.kern(D2,g,dt)
  K2 <- uniformization(G2,dt)
  ## Evaluate filter, returns likelihood value
  f <- hmm.filter(g,L,K1,K2,P)
  nllf <- -sum(log(f$psi))
  ##print(nllf)
  cat("\r HMM -log(L):",nllf); flush.console()
  nllf
}

hmm.smoother <- function(f,K1,K2,P){
  ## Smoothing the filtered estimates
  ## The equations for smoothing are presented in Pedersen et al. 2011, Oikos, Appendix
  T <- dim(f$phi)[2]
  row <- dim(f$phi)[3]
  col <- dim(f$phi)[4]
  
  smooth <- array(0,dim=dim(f$phi))
  smooth[,T,,] <- f$phi[,T,,]
  for(t in T:2){
    RAT <- smooth[,t,,]/(f$pred[,t,,]+1e-15)
    Rp1 <- as.vector(K1 %*% as.vector(RAT[1,,]))
    Rp2 <- as.vector(K2 %*% as.vector(RAT[2,,]))
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

show.movie <- function(distr){
  ## Plot the estimated probability distribution of the location as an animation
  
  graphics.off()
  nt <- dim(distr)[2]
  tr <- read.table('../sim/truedata.csv',header=TRUE,sep=',')
  dd <- apply(distr,c(2,3,4),sum)
  for(t in 1:nt){
    image(g$lon[1,],g$lat[,1],t(dd[t,,]),xlab='Longitude',ylab='Latitude')
    lines(tr$lon,tr$lat)
    Sys.sleep(0.1)
  }
}

sphmm.demo <- function(par0=c(8.908,10.27,1.152,0.0472,0.707,0.866),do.fit=FALSE,do.plot=TRUE,show.movie=FALSE){
  source('../sim/sstdb.R')
  require(Matrix) # For expm, matrix exponential and sparsity
  
  lsst <- read.table('../sim/lsstdata.csv',header=TRUE,sep=',')
  tr <- read.table('../sim/truedata.csv',header=TRUE,sep=',')
  iniloc <- read.table('../sim/iniloc.csv',header=TRUE,sep=',')
  
  ## Calculate time step
  dt <- date2time(lsst$date[2])-date2time(lsst$date[1])
  
  ## Setup discrete spatial grid
  ngrid <- c(50,49)
  g <- setup.grid(lsst,ngrid)
  row <- dim(g$lon)[1]
  col <- dim(g$lon)[2]
  
  ## Compute likelihood for observations
  L <- data.lik(lsst,iniloc,g)
  
  ## Number of time steps
  T <- length(lsst$lon)
  
  ## Fixed parameter values
  D1 <- par0[1:2]
  D2 <- par0[3:4]
  p <- par0[5:6]
  
  if(do.fit){
    guess <- c(log(10),log(10),log(0.5),log(0.5),log(0.95/0.05),log(0.95/0.05))
    fit <- nlm(neg.log.lik.fun,guess,g,L,dt)
    D1 <- exp(fit$estimate[1:2])
    D2 <- exp(fit$estimate[3:4])
    p <- 1/(1+exp(-fit$estimate[5:6]))
  }
  
  ## Setup transition matrices
  G1 <- make.kern(D1,g,dt)
  K1 <- uniformization(G1,dt)
  G2 <- make.kern(D2,g,dt)
  K2 <- uniformization(G2,dt)
  P <- matrix(c(p[1],1-p[1],1-p[2],p[2]),2,2,byrow=TRUE)
  
  ## Run smoother and filter
  f <- hmm.filter(g,L,K1,K2,P)
  s <- hmm.smoother(f,K1,K2,P)
  sphmm <- calc.track(s,g)
  sphmm$date <- lsst$date
  sphmm$p.resid <- apply(s,c(1,2),sum)[2,]
  
  if (do.plot){
    environment(plot.results) <- environment()
    plot.results(show.movie=show.movie)
  }
  
  invisible(sphmm)
}


plot.results <- function(save.plot=TRUE,show.movie=FALSE){
  ## Show movement as animation
  if(show.movie) show.movie(s)
  ## Plot behaviour
  pp <- apply(s,c(1,2),sum)[2,]
  ind <- seq(1,length(tr$b),by=72)
  graphics.off();
  if(save.plot) pdf('../spathmm/sphmm.pdf',width=7,height=8)
  par(mfrow=c(2,1))
  tr$date<-as.POSIXct(strptime(tr$date,"%Y-%m-%d %H:%M:%S",tz="GMT"))
  plot(I(b-1)~date,tr,col='grey',type='h',ylab='Probability of resident',main='Estimated behaviour',xlab='Date')
  lines(pp~tr$date[ind],lwd=2,col='red')
  xl <- c(-45,-30)
  yl <- c(11,32)
  plot(tr$lon,tr$lat,type='n',main='Estimated movements',ylab='Latitude',xlab='Longitude',xlim=xl,ylim=yl)
  rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], col = "steelblue1")
  lines(tr$lon,tr$lat,col='white')
  lines(sphmm$meanlon,sphmm$meanlat,col='black')
  points(tr$lon[1],tr$lat[1],bg='green',pch=21)
  TT <- length(tr$lon)
  points(tr$lon[TT],tr$lat[TT],bg='red',pch=21)
  if(save.plot) dev.off()
  
  sst <- numeric(T)
  ## Simple diagnostics plot ###
  ## Resample SST
  for(t in 1:T){
    time <- date2time(lsst$date[t])
    sst[t] <- sstdb(time,sphmm$meanlon[t],sphmm$meanlat[t])
  }
  
  if(save.plot) pdf('../plot/sphmmDiagn.pdf',width=7,height=7)
  if(!save.plot) dev.new()
  par(mfrow=c(2,2))
  ssterr <- sst-lsst$sst;
  sdsst <- sqrt(var(ssterr))
  ssterr <- ssterr/sdsst
  lonerr <- sphmm$meanlon-lsst$lon;
  sdlon <- sqrt(var(lonerr))
  lonerr <- lonerr/sdlon
  plot(tr$date[ind],ssterr,xlab='Date',ylab='SST residual',main='Residuals through time',ylim=c(-3,3),pch=".",cex=3)
  abline(h=c(-2,0,2),col='grey',lwd=2,lty=2)
  qqnorm(ssterr,main='QQ-plot of residuals',pch=".",cex=2)
  abline(a=0,b=1,col='grey',lwd=2)
  plot(tr$date[ind],lonerr,xlab='Date',ylab='Longitude residual',ylim=c(-3,3),pch=".",cex=3)
  abline(h=c(-2,0,2),col='grey',lwd=2,lty=2)
  qqnorm(lonerr,main='',pch=".",cex=2)
  abline(a=0,b=1,col='grey',lwd=2)
  if(save.plot) dev.off()
}
