% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/refinedBoundaries.R
\name{psbTeXtable}
\alias{psbTeXtable}
\title{Summarize Primary and Refined Secondary Boundaries in a TeX table}
\usage{
psbTeXtable(
  alpha,
  tVec,
  pOBF = TRUE,
  sOBF = FALSE,
  LanDeMets = FALSE,
  digits = 2
)
}
\arguments{
\item{alpha}{type I error probability.}

\item{tVec}{vector of relative information levels. The last element in the vector is 1.}

\item{pOBF}{type of primary boundary, \code{TURE} is the O'Brien-Fleming boundary, \code{FALSE} is the Pocock boundary.}

\item{sOBF}{type of secondary boundary, \code{TURE} is the O'Brien-Fleming boundary, \code{FALSE} is the Pocock boundary.}

\item{LanDeMets}{type of boundary, \code{TRUE} is the error spending approach, \code{FALSE} is the original approach.}

\item{digits}{number of digits after decimal point to display in the table.}
}
\value{
a TeX format table including both primary boundary and refined secondary boundary.
}
\description{
Primary boundaries and refined secondary boundaries are listed in a TeX table.
}
\details{
This function gives a TeX format table including both primary boundary and refined secondary boundary.
The number of digits after decimal point can be specified through parameter \code{digits}.
}
\examples{
require(mvtnorm)
require(ldbounds)
require(xtable)
psbTeXtable(alpha=0.025,tVec=c(1/2,3/4,1),pOBF=TRUE,sOBF=FALSE,LanDeMets=FALSE)

}
\references{
Glimm, E., Maurer, W., and Bretz, F. (2010). Hierarchical testing of multiple endpoints in group-sequential trials. \emph{Statistics in Medicine} \bold{29}, 219-228.

Hung, H. M. J., Wang, S.-J., and O'Neill, R. (2007). Statistical considerations for testing multiple endpoints in group sequential or adaptive clinical trials. \emph{Journal of Biopharmaceutical Statistics} \bold{17}, 1201-1210.

Jennison, C. and Turnbull, B. W. (2000). \emph{Group Sequential Methods with Applications to Clinical Trials}. Chapman and Hall/CRC, New York.

Lan, K. K. G., and Demets, D. L. (1983). Discrete sequential boundaries for clinical trials. \emph{Biometrika} \bold{70}, 659-663.

O'Brien, P. C., and Fleming, T. R. (1979). A multiple testing procedure for clinical trials. \emph{Biometrics} \bold{35}, 549-556.

Pocock, S. J. (1977). Group sequential methods in the design and analysis of clinical trials. \emph{Biometrika} \bold{64}, 191-199.

Tamhane, A. C., Mehta, C. R., and Liu, L. (2010). Testing a primary and a secondary endpoint in a group sequential design. \emph{Biometrics} \bold{66}, 1174-1184.

Tamhane, A. C., Gou, J., Jennison, C., Mehta, C. R., and Curto, T. (2018). A gatekeeping procedure to test a primary and a secondary endpoint in a group sequential design with multiple interim looks. \emph{Biometrics}, 74, 40-48.
}
\author{
Jiangtao Gou

Fengqing (Zoe) Zhang
}
