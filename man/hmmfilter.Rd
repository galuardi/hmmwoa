\name{hmmfilter}
\alias{hmmfilter}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
hmmfilter(sstL, fmat, ks = 29, sigma = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{sstL}{
%%     ~~Describe \code{sstL} here~~
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
function (sstL, fmat, ks = 29, sigma = 1) 
{
    allpost = sstL[[3]] * NaN
    lon = L$lon
    lat = L$lat
    xidx = which.min((fmat[1, 8] - (lon))^2)
    yidx = which.min((fmat[1, 9] - lat)^2)
    allpost[xidx, yidx, 1] = 1
    len = nrow(fmat)
    xidx = which.min((fmat[len, 8] - (lon))^2)
    yidx = which.min((fmat[len, 9] - lat)^2)
    allpost[xidx, yidx, len] = 1
    post = allpost[, , 1]
    post[is.nan(post)] = 0
    datediff = c(1, diff(as.numeric(mdy.date(fmat[, 2], fmat[, 
        3], fmat[, 1]))))
    f = gausskern(ks[1], sigma[1])
    f = f/sum(f)
    post = Image(post)
    P <- filter2(post, f)
    post = P@.Data
    post = normalise(post)
    normaliser = numeric(nrow(fmat) - 1) + 1
    time1 = Sys.time()
    for (i in 2:(nrow(fmat) - 1)) {
        print(i)
        if (fmat$behav[i - 1] == 1) {
            sig = sigma[1]
            kern = ks[1]
        }
        else {
            sig = sigma[2]
            kern = ks[2]
        }
        post[is.nan(post)] = 0
        f = gausskern(kern, sig)
        f = f/sum(f)
        P <- filter2(post, f)
        P = P@.Data
        post = P * sstL[[3]][, , i]
        normaliser[i - 1] = sum(post, na.rm = T)
        post = normalise(post)
        allpost[, , i] = post
    }
    print(Sys.time() - time1)
    allpost
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
