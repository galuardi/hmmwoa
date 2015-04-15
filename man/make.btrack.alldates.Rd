\name{make.btrack.alldates}
\alias{make.btrack.alldates}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
make.btrack.alldates(track, psat)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{track}{
%%     ~~Describe \code{track} here~~
}
  \item{psat}{
%%     ~~Describe \code{psat} here~~
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
function (track, psat) 
{
    require(date)
    alldates = seq(psat$day0, psat$dayT)
    alldates.dmy = date.mdy(alldates)
    alldates.dmy = cbind(alldates.dmy[[2]], alldates.dmy[[1]], 
        alldates.dmy[[3]])
    bdates = mdy.date(track[, 2], track[, 3], track[, 1])
    bidx = match(bdates, alldates)
    output = as.data.frame(matrix(NA, nrow = length(alldates), 
        ncol = ncol(track)))
    names(output) = names(track)
    output[, 1:3] = alldates.dmy
    output[bidx, 4:ncol(track)] = track[, 4:ncol(track)]
    output
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
