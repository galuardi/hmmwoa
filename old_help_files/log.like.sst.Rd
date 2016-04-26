\name{log.like.sst}
\alias{log.like.sst}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
log.like.sst(ingrid, datax = as.numeric(vec[11]), sigma = as.numeric(vec[12]))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{ingrid}{
%%     ~~Describe \code{ingrid} here~~
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
function (ingrid, datax = as.numeric(vec[11]), sigma = as.numeric(vec[12])) 
{
    if (is.na(sigma)) 
        sigma = 5
    tempL = ingrid * 0
    tdim = dim(ingrid)
    tempL = dlnorm(log(ingrid), log(datax), log(sigma), log = F)
    tempL
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
