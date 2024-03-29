% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/refinedBoundaries.R
\name{ldSecondaryBoundary}
\alias{ldSecondaryBoundary}
\title{Calculate Refined Secondary Boundary, Error Spending Approach}
\usage{
ldSecondaryBoundary(
  alpha,
  tVec,
  cvec,
  locPeak,
  type = 2,
  initIntvl = c(0.6, 4)
)
}
\arguments{
\item{alpha}{original significance level.}

\item{tVec}{information vector.}

\item{cvec}{primary group sequential boundary.}

\item{locPeak}{location of maximum, a number between 1 and the number of interims.}

\item{type}{type of the test procedure for the secondary endpoint. O'Brien- Fleming (OBF) type error spending funciton is 1, Pocock (POC) type error spending funciton is 2.}

\item{initIntvl}{computing paramter, a pair of numbers containing the end-points of the interval to be searched for the root.}
}
\value{
refined secondary boundaries.
}
\description{
Refined secondary boundaries are calculated by using the error spending approach.
}
\details{
This function calculates the refined secondary boundaries of any Lan-DeMets error spending boundary based on the primary boundaries.
}
\examples{
\dontrun{
require(mvtnorm)
require(ldbounds)
K <- 6;
tVec <- c(140,328,453,578,659,1080)/1080;
alpha = 0.025;
cvec.obf <- ldbounds::ldBounds(t=tVec,iuse=c(1),alpha=c(alpha),sides = 1);
cvec <- cvec.obf$upper.bounds;
secbound <- ldSecondaryBoundary(alpha,tVec,cvec,locPeak=4,type=2,
    initIntvl=c(0.8,8))
}

}
\references{
Lan, K. K. G., and Demets, D. L. (1983). Discrete sequential boundaries for clinical trials. \emph{Biometrika} \bold{70}, 659-663.

Tamhane, A. C., Gou, J., Jennison, C., Mehta, C. R., and Curto, T. (2017+). A gatekeeping procedure to test a primary and a secondary endpoint in a group sequential design with multiple interim looks. \emph{Biometrics}, 74, 40-48.
}
\seealso{
\code{secondaryBoundary}, \code{secondaryBoundaryVecLD}
}
\author{
Jiangtao Gou
}
