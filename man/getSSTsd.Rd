\name{getSSTsd}
\alias{getSSTsd}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
getSSTsd(psat, surf = -11)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{psat}{
%%     ~~Describe \code{psat} here~~
}
  \item{surf}{
%%     ~~Describe \code{surf} here~~
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
function (psat, surf = -11) 
{
    tzidx = psat$Z >= (surf)
    tzidx[is.na(tzidx)] = F
    temp = psat$T
    temp[tzidx == F] = NA
    tsd = apply(temp, 1, sd, na.rm = T)
    tsd[tsd == 0] = NA
    tsd
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
