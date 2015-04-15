\name{mpt.viterbi}
\alias{mpt.viterbi}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
mpt.viterbi(allpost2, L, fmat, sst, D = c(100, 500), D2s = 0.09)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{allpost2}{
%%     ~~Describe \code{allpost2} here~~
}
  \item{L}{
%%     ~~Describe \code{L} here~~
}
  \item{fmat}{
%%     ~~Describe \code{fmat} here~~
}
  \item{sst}{
%%     ~~Describe \code{sst} here~~
}
  \item{D}{
%%     ~~Describe \code{D} here~~
}
  \item{D2s}{
%%     ~~Describe \code{D2s} here~~
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
function (allpost2, L, fmat, sst, D = c(100, 500), D2s = 0.09) 
{
    land = make.landmask(sst)$mask
    print(sprintf("Number of days: \%i\n", dim(L[[3]][3])))
    s = D * D2s
    rrow = dim(L[[3]])[1]
    ccol = dim(L[[3]])[2]
    icalc = dim(L[[3]])[3]
    numnames = 1
    unc = sqrt(2 * s[1])
    ks = ceiling(unc * 10 + 1)
    ks = ks + mod(ks, 2) + 1
    ks1 = max(15, ks)
    kern1 = gausskern(ks1, unc)
    unc = sqrt(2 * s[2])
    ks = ceiling(unc * 10 + 1)
    ks = ks + mod(ks, 2) + 1
    ks2 = max(15, ks)
    kern2 = gausskern(ks2, unc)
    mpt.maplong = matrix(L$lon, ccol, rrow, byrow = T)
    mpt.maplat = matrix(L$lat, ccol, rrow, byrow = F)
    dlong = diff(L$lon)[1]
    dlat = diff(L$lat)[1]
    R = mapmatrix(L$lat[1], L$lon[1], dlat, dlong)
    smatrix = allpost2
    smatrix[is.nan(smatrix)] = 0
    M = smatrix[, , 1]
    M = log(M)
    subject = (M != -Inf) * 1
    zro = matrix(0, rrow, ccol)
    theend = icalc
    Tprevx = numeric(rrow * ccol * theend)
    dim(Tprevx) = c(rrow, ccol, theend)
    Tprevy = Tprevx
    xidx = which.min((fmat[1, 8] - (L$lon))^2)
    yidx = which.min((fmat[1, 9] - L$lat)^2)
    Tprevx[xidx, yidx, 1] = xidx
    Tprevy[xidx, yidx, 1] = yidx
    Ltotal = smatrix * 0 + 1
    for (j in 1:numnames) {
        Ltotal = Ltotal * smatrix
    }
    Ltotal[is.nan(Ltotal)] = 0
    Ltotal[, , 1] = smatrix[, , 1]
    print("Starting iterations...")
    print(sprintf("Day   1 -  9..."))
    Tx = numeric(rrow * ccol * theend)
    dim(Tx) = c(rrow, ccol, theend)
    Ty = Tx
    for (j in 2:(theend - 1)) {
        if (!mod(j, 10)) {
            print(sprintf("\n done! time = \%4.3f", Sys.time()))
            print(sprintf("Day \%3.0i -\%3.0i...", j, j + 9))
            time1
        }
        Mtemp = log(zro)
        Ttempx = -1 + zro
        Ttempy = Ttempx
        time1 = Sys.time()
        if (fmat$behav[j - 1] == 1) {
            ks = ks1
            kern = kern1
        }
        else {
            ks = ks2
            kern = kern2
        }
        for (xx in 1:ccol) {
            for (yy in 1:rrow) {
                if (as.logical(subject[yy, xx])) {
                  kminlat = 1 + max(ceiling(ks/2) - yy, 0)
                  kmaxlat = min(ks - (yy + floor(ks/2) - rrow), 
                    ks)
                  kminlong = 1 + max(ceiling(ks/2) - xx, 0)
                  kmaxlong = min(ks - (xx + floor(ks/2) - ccol), 
                    ks)
                  klat = kminlat:kmaxlat
                  klong = kminlong:kmaxlong
                  mminlat = max(yy - floor(ks/2), 1)
                  mmaxlat = min(yy + floor(ks/2), rrow)
                  mminlong = max(xx - floor(ks/2), 1)
                  mmaxlong = min(xx + floor(ks/2), ccol)
                  mlat = mminlat:mmaxlat
                  mlong = mminlong:mmaxlong
                  B = log(Ltotal[mlat, mlong, j - 1] * kern[klat, 
                    klong])
                  Msub = B + M[yy, xx]
                  Mupdate = Mtemp[mlat, mlong]
                  Txupdate = Ttempx[mlat, mlong]
                  Tyupdate = Ttempy[mlat, mlong]
                  update = (Mupdate < Msub)
                  update[is.na(update)] = F
                  Mupdate[update] = Msub[update]
                  Txupdate[update] = xx
                  Tyupdate[update] = yy
                  Mtemp[mlat, mlong] = Mupdate
                  Ttempx[mlat, mlong] = Txupdate
                  Ttempy[mlat, mlong] = Tyupdate
                }
            }
        }
        Mtemp[land == 1] = -Inf
        subject = (Mtemp != -Inf) * 1
        for (xx in 1:ccol) {
            for (yy in 1:rrow) {
                if (as.logical(subject[yy, xx])) {
                  Tx[yy, xx, 1:(j - 1)] = Tprevx[Ttempy[yy, xx], 
                    Ttempx[yy, xx], 1:(j - 1)]
                  Ty[yy, xx, 1:(j - 1)] = Tprevy[Ttempy[yy, xx], 
                    Ttempx[yy, xx], 1:(j - 1)]
                  Tx[yy, xx, j] = xx
                  Ty[yy, xx, j] = yy
                }
            }
        }
        print(Sys.time() - time1)
        M = Mtemp
        Tprevx = Tx
        Tprevy = Ty
    }
    Sys.time() - time1
    M[land == 1] = Inf * -1
    val = max(M)
    ind = which.max(M)
    txy = ind2sub(c(rrow, ccol), ind)
    xm = txy[2]
    ym = txy[1]
    mpt.long = Tx[ym, xm, ]
    mpt.long_clean = mpt.long
    mpt.lat = Ty[ym, xm, ]
    mpt.lat_clean = mpt.lat
    mpt.lat = mpt.lat_clean = jitter(mpt.lat)
    mpt.long = mpt.long_clean = jitter(mpt.long)
    mpt = pixtomap(R, mpt.long_clean, mpt.lat_clean)
    mpt[, 1] = mpt[, 1]
    mpt[1, ] = as.numeric(fmat[1, 8:9])
    mpt[nrow(mpt), ] = as.numeric(fmat[nrow(fmat), 8:9])
    mpt
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
