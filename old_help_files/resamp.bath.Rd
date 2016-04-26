\name{resamp.bath}
\alias{resamp.bath}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
resamp.bath(sstL, bathy)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{sstL}{
%%     ~~Describe \code{sstL} here~~
}
  \item{bathy}{
%%     ~~Describe \code{bathy} here~~
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
function (sstL, bathy) 
{
    require(raster)
    xmin = min(sstL$lon)
    xmax = max(sstL$lon)
    ymin = min(sstL$lat)
    ymax = max(sstL$lat)
    b2 = crop.grid(xmin, xmax, ymin, ymax, bathy)
    s = raster(nrow = nrow(sstL$sstL), ncol = ncol(sstL$sstL))
    bath2 = raster(ncol = ncol(b2[[3]]), nrow = nrow(b2[[3]]))
    bath2 = setValues(bath2, as.vector(t(b2[[3]])))
    s <- resample(bath2, s, method = "bilinear")
    extent(s) = extent(c(min(sstL$lon) - 360, max(sstL$lon) - 
        360, min(sstL$lat), max(sstL$lat)))
    bath2 = list(lon = sstL$lon - 360, lat = sstL$lat, data = matrix(s@data@values, 
        nrow = nrow(sstL$sstL), ncol = ncol(sstL$sstL), byrow = T))
    bath2
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
