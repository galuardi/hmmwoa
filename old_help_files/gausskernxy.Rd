\name{gausskernxy}
\alias{gausskernxy}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
gausskernxy(siz, xsigma, ysigma, muadv = 0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{siz}{
%%     ~~Describe \code{siz} here~~
}
  \item{xsigma}{
%%     ~~Describe \code{xsigma} here~~
}
  \item{ysigma}{
%%     ~~Describe \code{ysigma} here~~
}
  \item{muadv}{
%%     ~~Describe \code{muadv} here~~
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
function (siz, xsigma, ysigma, muadv = 0) 
{
    x = 1:round(siz)
    mu = c(mean(x), mean(x)) + muadv
    fx = (matrix(exp((-0.5 * (x - mu[1])/xsigma)^2))/(sqrt(2 * 
        pi) * xsigma))
    options(digits = 5)
    fx = exp(-0.5 * ((x - mu[1])/xsigma)^2)/sqrt((2 * pi) * xsigma)
    fy = exp(-0.5 * ((x - mu[2])/ysigma)^2)/sqrt((2 * pi) * ysigma)
    fx[!is.finite(fx)] = 0
    fy[!is.finite(fy)] = 0
    kern = (fx \%*\% t(fy))
    kern = kern/(sum(sum(kern, na.rm = T), na.rm = T))
    kern[is.nan(kern)] = 0
    kern
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
