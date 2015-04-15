\name{like.sst.crop}
\alias{like.sst.crop}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
like.sst.crop(zgrid = zgrid, xmin = NULL, xmax = NULL, ymin = NULL, ymax = NULL, lonvec = NULL, latvec = NULL, datax = vec, sigma = 2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{zgrid}{
%%     ~~Describe \code{zgrid} here~~
}
  \item{xmin}{
%%     ~~Describe \code{xmin} here~~
}
  \item{xmax}{
%%     ~~Describe \code{xmax} here~~
}
  \item{ymin}{
%%     ~~Describe \code{ymin} here~~
}
  \item{ymax}{
%%     ~~Describe \code{ymax} here~~
}
  \item{lonvec}{
%%     ~~Describe \code{lonvec} here~~
}
  \item{latvec}{
%%     ~~Describe \code{latvec} here~~
}
  \item{datax}{
%%     ~~Describe \code{datax} here~~
}
  \item{sigma}{
%%     ~~Describe \code{sigma} here~~
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
function (zgrid = zgrid, xmin = NULL, xmax = NULL, ymin = NULL, 
    ymax = NULL, lonvec = NULL, latvec = NULL, datax = vec, sigma = 2) 
{
    vec = as.numeric(datax)
    data.lon = as.numeric(vec[8])
    data.lat = as.numeric(vec[9])
    data.sst = as.numeric(vec[11])
    lat = latvec
    lon = lonvec
    xmin = data.lon - 3
    xmax = data.lon + 3
    ymin = data.lat - 3
    ymax = data.lat + 3
    bxlim = c(xmin, xmax)
    bylim = c(ymin, ymax)
    blon = which(lon >= bxlim[1] & lon <= bxlim[2])
    blat = which(lat >= bylim[1] & lat <= bylim[2])
    bout = zgrid[blon, blat]
    lout = zgrid * 0
    lout[blon, blat] = like.sst(bout, datax = data.sst, sigma = sigma)
    lout
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
