\name{make.like}
\alias{make.like}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
make.like(sst, bathy, fmat, lims = c(-85, -65, 30, 45), uselog = F)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{sst}{
%%     ~~Describe \code{sst} here~~
}
  \item{bathy}{
%%     ~~Describe \code{bathy} here~~
}
  \item{fmat}{
%%     ~~Describe \code{fmat} here~~
}
  \item{lims}{
%%     ~~Describe \code{lims} here~~
}
  \item{uselog}{
%%     ~~Describe \code{uselog} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (sst, bathy, fmat, lims = c(-85, -65, 30, 45), uselog = F) 
{
    require(raster)
    lon = sst$lon
    lat = sst$lat
    xmin = which.min(((lon) - lims[1])^2)
    xmax = which.min(((lon) - lims[2])^2)
    ymin = which.min(((lat) - lims[3])^2)
    ymax = which.min(((lat) - lims[4])^2)
    tempL = numeric(length = c(dim(sst$DATA)[1] * dim(sst$DATA)[2] * 
        nrow(fmat)))
    dim(tempL) = c(dim(sst$DATA)[1], dim(sst$DATA)[2], nrow(fmat))
    tempL = tempL[xmin:xmax, ymin:ymax, ]
    land = sst$DATA[xmin:xmax, ymin:ymax, 1] * 0 + 2
    land[is.nan(land)] = 0
    land[land == 2] = 1
    time1 = Sys.time()
    for (i in 2:(nrow(fmat) - 1)) {
        sigma = fmat$sstsd[i]
        if (class(sst$sstdates)[1] == "POSIXct") {
            tdate = ISOdate(fmat$Year[i], fmat$Month[i], fmat$Day[i])
            zidx = which.min((as.numeric(tdate) - as.numeric(sst$sstdates))^2)
        }
        else {
            tdate = mdy.date(fmat$Month[i], fmat$Day[i], fmat$Year[i])
            zidx = which.min((tdate - sst$sstdates)^2)
        }
        if (uselog == T) {
            tempL[, , i] = normalise(log.like.sst(sst$DATA[xmin:xmax, 
                ymin:ymax, zidx], as.numeric(fmat[i, "SST"]), 
                sigma = sigma, log = F))
        }
        else {
            tempL[, , i] = like.sst.tnorm(ingrid = sst$DATA[xmin:xmax, 
                ymin:ymax, zidx], datax = as.numeric(fmat[i, 
                "SST"]), sigma = sigma)
        }
    }
    sstL = list(lon = sst$lon[xmin:xmax], lat = sst$lat[ymin:ymax], 
        sstL = tempL)
    print(paste(Sys.time() - time1, "seconds to make sst likelihood"))
    rm(tempL)
    bath2 = resamp.bath(sstL, bathy)
    tempB = sstL$sstL * NaN
    for (i in 2:(nrow(fmat) - 1)) {
        zidx = bath2$data <= as.numeric(fmat[i, 10])
        bb = tempB[, , i]
        bb[zidx == T] = 1
        bb[zidx == F] = 0
        tempB[, , i] = bb
    }
    print(paste(Sys.time() - time1, "seconds to make depth likelihood"))
    list(lon = sstL$lon, lat = sstL$lat, L = sstL$sstL * tempB)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
