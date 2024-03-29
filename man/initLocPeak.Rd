% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/refinedBoundaries.R
\name{initLocPeak}
\alias{initLocPeak}
\title{Find the Location of Maximum, Standard OBF and POC}
\usage{
initLocPeak(alpha, tVec, cvec, type = 2, initIntvl = c(1, 4))
}
\arguments{
\item{alpha}{type I error.}

\item{tVec}{information vector.}

\item{cvec}{primary group sequential boundary.}

\item{type}{type of the test procedure for the secondary endpoint. O'Brien- Fleming (OBF) type error spending funciton is 1, Pocock (POC) type error spending funciton is 2.}

\item{initIntvl}{computing paramter, a pair of numbers containing the end-points of the interval to be searched for the root.}
}
\value{
location of maximum, a number between 1 and the number of interims
}
\description{
Calculate the location of maximal tyep I error of the standard O'Brien-Fleming and Pocock refined secondary boundaries.
}
\details{
This function search the location of the maximal point, in order to calculate the standard (origiinal) O'Brien-Fleming (OBF) and Pocock (POC) refined secondary boundaries.
}
\examples{
require(mvtnorm)
K <- 8
gammaVec <- sqrt((1:K)/K)
tVec <- gammaVec^2
alpha <- 0.025
c <- 2.072274
cvec <- c/gammaVec
loc <- initLocPeak(alpha,tVec,cvec,type=2,initIntvl=c(1,3))

}
\references{
O'Brien, P. C., and Fleming, T. R. (1979). A multiple testing procedure for clinical trials. \emph{Biometrics} \bold{35}, 549-556.

Pocock, S. J. (1977). Group sequential methods in the design and analysis of clinical trials. \emph{Biometrika} \bold{64}, 191-199.

Tamhane, A. C., Gou, J., Jennison, C., Mehta, C. R., and Curto, T. (2018). A gatekeeping procedure to test a primary and a secondary endpoint in a group sequential design with multiple interim looks. \emph{Biometrics}, to appear.
}
\seealso{
\code{SecondaryBoundary}, \code{ldInitLocBeak}
}
\author{
Jiangtao Gou
}
