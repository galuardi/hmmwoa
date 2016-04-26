\name{ind2sub}
\alias{ind2sub}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
ind2sub(siz, ndx)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{siz}{
%%     ~~Describe \code{siz} here~~
}
  \item{ndx}{
%%     ~~Describe \code{ndx} here~~
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
function (siz, ndx) 
{
    nout = 2
    siz = as.double(siz)
    len = end(siz)[1]
    if (length(siz) <= nout) {
        siz = c(siz, numeric(nout - length(siz)))
    }
    else {
        siz = c(siz[1:(nout - 1)], prod(siz[nout:len]))
    }
    n = length(siz)
    k = c(1, cumprod(siz[1:(len - 1)]))
    varargout = numeric(n)
    for (i in rev(1:n)) {
        vi = rem(ndx - 1, k[i]) + 1
        vj = (ndx - vi)/k[i] + 1
        varargout[i] = vj
        ndx = vi
    }
    varargout
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
