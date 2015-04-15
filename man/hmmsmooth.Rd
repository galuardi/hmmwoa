\name{hmmsmooth}
\alias{hmmsmooth}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
hmmsmooth(allpost, fmat, ks, sigma)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{allpost}{
%%     ~~Describe \code{allpost} here~~
}
  \item{fmat}{
%%     ~~Describe \code{fmat} here~~
}
  \item{ks}{
%%     ~~Describe \code{ks} here~~
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
function (allpost, fmat, ks, sigma) 
{
    time1 = Sys.time()
    allpost[is.nan(allpost)] = 0
    allpost2 = allpost * 0
    post = allpost[, , dim(allpost)[3]]
    allpost2[, , dim(allpost)[3]] = post
    allpost2[, , 1] = allpost[, , 1]
    for (i in rev(2:dim(allpost)[3])) {
        if (fmat$behav[i] == 1) {
            sig = sigma[1]
            kern = ks[1]
        }
        else {
            sig = sigma[2]
            kern = ks[2]
        }
        print(i)
        post[is.nan(post)] = 0
        f = gausskern(kern, sig)
        f = f/sum(f)
        rat = filter2(allpost[, , i - 1], f)
        ratio = (allpost2[, , i]/(rat@.Data + 1e-20))
        S <- filter2(ratio, f)
        S <- S@.Data
        post = S * allpost[, , i - 1]
        allpost2[, , i - 1] = normalise(post)
    }
    print(Sys.time() - time1)
    allpost2
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
