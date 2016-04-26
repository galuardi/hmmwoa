\name{crop.grid}
\alias{crop.grid}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
crop.grid(xmin, xmax, ymin, ymax, ingrid)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
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
  \item{ingrid}{
%%     ~~Describe \code{ingrid} here~~
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
function (xmin, xmax, ymin, ymax, ingrid) 
{
    bxlim = c(xmin, xmax)
    bylim = c(ymin, ymax)
    blon = which(ingrid$lon >= bxlim[1] & ingrid$lon <= bxlim[2])
    blat = which(ingrid$lat >= bylim[1] & ingrid$lat <= bylim[2])
    bath.out = list()
    bath.out$lon = ingrid$lon[blon]
    bath.out$lat = ingrid$lat[blat]
    bath.out$data = t(ingrid$data[blat, blon])
    bath.out
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
