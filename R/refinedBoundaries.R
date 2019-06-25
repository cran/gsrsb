##########
# A gatekeeping test on a primary and a secondary endpoint in a group sequential design with multiple interim looks
# R program for Refined Boundaries
# Version 1.1
# Programmer: Jiangtao Gou
# Version 1.0: 2017/Mar/14-Apr/12
# Version 1.1: 2019/Jun/24
# Reference: Tamhane AC, Gou J, Jennison C, Mehta CR, Curto T.
#   A Gatekeeping Test on a Primary and a Secondary Endpoint in a Group
#   Sequential Design with Multiple Interim Looks.
#   Biometrics (2018).
##########

require(mvtnorm)
require(ldbounds)
require(xtable)

#
# Function 1
# function cdBoundary
# Generate lower and upper bounds for programs calculating the secondary endpoint's type I error when rho=1
# Input:
#  c boundaries cvec (vector),
#  d boundaries dvec (vector),
#  gamma, square root of information gammaVec (vector),
#  effect size of secondary endpoint dlt (number),
#  lower or upper bound upper (binary).
# Output:
#  lower or upper bounds (vector).
# Dependency:
#  None.
#
#' Lower and Upper Bounds Generator
#'
#' Generate lower and upper bounds for programs calculating the secondary endpoint's type I error when the correlation rho between the primary endpoint and the secondary endpoint equals 1.
#'
#' @param cvec primary boundary.
#' @param dvec secondary boundary.
#' @param gammaVec square root of information vector.
#' @param dlt test statistic of the primary endpoint follows a normal distribution with mean \code{dlt} and standard deviation 1.
#' @param upper type of bounds, upper bound is \code{TRUE}, lower bound is \code{FALSE}.
#' @return lower and upper bounds for programs calculating the secondary endpoint's type I error when the correlation rho is 1.
#'
#' @author Jiangtao Gou
#'
#' @details
#' This function generates upper and lower bounds for further computation. For more details, refer to Tamhane et al. (2018, Biometrics), section 4.2.
#'
#' @references
#'  Tamhane, A. C., Gou, J., Jennison, C., Mehta, C. R., and Curto, T. (2018). A gatekeeping procedure to test a primary and a secondary endpoint in a group sequential design with multiple interim looks. \emph{Biometrics}, 74(1), 40-48.
#'
#' @examples
#' cvec <- rep(1.992,3)
#' dvec <- c(1.535*sqrt(3),1.535*sqrt(3/2),1.535)
#' gammaVec <- c(sqrt(1/3),sqrt(2/3),1)
#' dlt <- 2
#' uBoundary <- cdBoundary(cvec, dvec, gammaVec, dlt, upper=TRUE)
#'
#' @export
#' @import mvtnorm
#' @import stats
#' @import ldbounds
#'
cdBoundary <- function (cvec, dvec, gammaVec, dlt, upper=TRUE) {
  K <- length(cvec);
  if (upper == TRUE) {
    uBoundary <- rep(Inf,K);
    uBoundary[1:(K-1)] = cvec[1:(K-1)] - gammaVec[1:(K-1)]*dlt;
    return(uBoundary);
  } else {
    lBoundary <- rep(-Inf,K);
    lBoundary[K] = max(cvec[K] - gammaVec[K]*dlt,dvec[K]);
    return(lBoundary);
  }
}
# End of function cdBoundary
#

#
# Function 2
# function genCorrMat
# Correlation matrix generator.
# Input:
#  gammaVec: a vector which contains gamma_{1}, ..., gamma_{K-1}, gamma_{K},
#  type: For primary endpoint calculation, type is 1,
#   the returned matrix is K by K. For secondary endpoint calculation,
#   type is 2, the returned matrix is (K+1) by (K+1).
#  rhoPS: correlation between primary and secondary endpoints
# Output:
#  Correlation matrix generator, K by K for primary endpoint, (K+1) by (K+1) for secondary endpoint.
# Dependency:
#  None.
#
#' Correlation Matrix Generator
#'
#' Generate correlation matrix between standardized sample mean test statistics for the two endpoint at different looks.
#'
#' @author Jiangtao Gou
#' @author Fengqing (Zoe) Zhang
#'
#' @param gammaVec a vector which contains gamma_(1), ..., gamma_(K-1), gamma_(K), square root of information vector.
#' @param type type of primary or secondary endpoint. For primary endpoint calculation, \code{type} is 1, the returned matrix is K by K. For secondary endpoint calculation, \code{type} is 2, the returned matrix is (K+1) by (K+1).
#' @param rhoPS correlation between primary and secondary endpoints.
#' @return correlation matrix, K by K for primary endpoint, (K+1) by (K+1) for secondary endpoint, where K is the number of interims.
#'
#' @details
#' This function generates correlation matrix between different mean statistics. For more details, refer to Tamhane et al. (2018, Biometrics), section 2.
#' @export
#' @import mvtnorm
#' @import stats
#' @import ldbounds
#'
#' @references
#'  Tamhane, A. C., Gou, J., Jennison, C., Mehta, C. R., and Curto, T. (2018). A gatekeeping procedure to test a primary and a secondary endpoint in a group sequential design with multiple interim looks. \emph{Biometrics}, 74(1), 40-48.
#'
#' @examples
#' corrMat <- genCorrMat(gammaVec=c(sqrt(1/3),sqrt(2/3),1), type=2, rhoPS = 0.3)
#'
genCorrMat <- function (gammaVec, type, rhoPS = 0) {
  K <- length(gammaVec);
  if (type == 1) {
    corrMat <- matrix(rep(1,K*K),nrow=K,ncol=K,byrow=TRUE);
    for (i in 1:K) {
      for (j in 1:K){
        corrMat[i,j] <- min(gammaVec[i],gammaVec[j])/max(gammaVec[i],gammaVec[j]);
      }
    }
    return(corrMat);
  } else if (type == 2) {
    corrMat <- matrix(rep(1,(K+1)*(K+1)),nrow=K+1,ncol=K+1,byrow=TRUE);
    for (i in 1:K) {
      for (j in 1:K){
        corrMat[i,j] <- min(gammaVec[i],gammaVec[j])/max(gammaVec[i],gammaVec[j]);
      }
    }
    for (i in 1:K) {
      corrMat[i,K+1] <- rhoPS*gammaVec[i]/gammaVec[K];
      corrMat[K+1,i] <- rhoPS*gammaVec[i]/gammaVec[K];
    }
    return(corrMat);
  } else {
    print("Type is 1 (for primary endpoint calculation) or 2 (for secondary endpoint calculation).");
    return(type);
  }
}
# End of function genCorrMat
#

#
# Function 3
# function ldPrimaryBoundary
# Primary boundaries calculation of Lan-DeMets OBF and POC
# Input:
#  tVec: a vector of information, gammaVec = sqrt(tVec)
#  alpha: significance level
#  type: type 1 OBF, type 2 POC
#  initIntvl: paramter for function uniroot (two numbers)
# Output:
#  Primary boundaries (vector)
# Dependency:
#  genCorrMat,
#  pmvnorm (package mvtnorm).
#  uniroot (package stats)
#
#' Calculate Primary Boundaries, the Error Spending Approach
#'
#' Primary boundaries calculation of Lan-DeMets OBF and POC.
#'
#' @export
#' @import mvtnorm
#' @import stats
#' @import ldbounds
#'
#' @author Jiangtao Gou
#'
#' @seealso \code{primaryBoundary}
#'
#' @param tVec a vector of information, gammaVec = sqrt(tVec).
#' @param alpha significance level
#' @param type type of sequential procedure. OBF is 1, POC is 2.
#' @param initIntvl paramter for function uniroot (two numbers)
#' @return a vector of primary boundaries.
#'
#' @references
#' Lan, K. K. G., and Demets, D. L. (1983). Discrete sequential boundaries for clinical trials. \emph{Biometrika} \bold{70}, 659-663.
#'
#'  Tamhane, A. C., Gou, J., Jennison, C., Mehta, C. R., and Curto, T. (2018). A gatekeeping procedure to test a primary and a secondary endpoint in a group sequential design with multiple interim looks. \emph{Biometrics}, 74(1), 40-48.
#'
ldPrimaryBoundary <- function (tVec,alpha,type=1,initIntvl=c(0.8,8)) {
  K <- length(tVec);
  cVec <- rep(0,K);
  gammaVec <- sqrt(tVec);
  target <- function (ckk,cVec,gammaVec,alpha) {
    kk <- length(gammaVec);
    lowerB <- rep(-Inf,kk);
    upperB <- c(cVec,ckk);
    meanV <- rep(0,kk);
    corrM <- genCorrMat(gammaVec,1);
    result <- pmvnorm(lowerB, upperB, meanV, corrM, algorithm=Miwa(steps=128));
    return(result - 1 + alpha);
  }
  if (type == 1) {
    # print('It is OBF');
    alphaVec <- 2 - 2*pnorm(qnorm(1-alpha/2)/gammaVec);
    cVec[1] <- qnorm(1-alphaVec[1]);
    for (i in 2:K) {
      result.OF <- uniroot(target, lower = initIntvl[1], upper = initIntvl[2], tol = 2.5e-16, cVec = cVec[1:(i-1)], gammaVec = gammaVec[1:i], alpha = alphaVec[i]);
      cVec[i] = result.OF$root;
    }
    return(cVec);
  } else if (type == 2) {
    # print('It is POC');
    alphaVec <- alpha*log(1 + (exp(1)-1)*tVec);
    cVec[1] <- qnorm(1-alphaVec[1]);
    for (i in 2:K) {
      result.PO <- uniroot(target, lower = initIntvl[1], upper = initIntvl[2], tol = 2.5e-16, cVec = cVec[1:(i-1)], gammaVec = gammaVec[1:i], alpha = alphaVec[i]);
      cVec[i] = result.PO$root;
    }
    return(cVec);
  } else {
    print("Only OBF (1) and POC (2) are currently available.");
    return(K);
  }
}
# End of function ldPrimaryBoundary
#

#
# Function 4
# function ldSecControl
# Secondary boundaries calculation of Lan-DeMets OBF and POC
# Input:
#  ap: significance level for the primary endpoint
#  alpha: targeted significance level for the secondary endpoint
#  cvec: a vector of calculated primary boundaries
#  tVec: a vector of information, gammaVec = sqrt(tVec)
#  ExtrmLoc: an integer between 1 and K, locate the maximum of type I error of secondary endpoint
#  type: Type 1 OBF d, Type 2 POC d.
# Output:
#  Difference between alpha and the calculated error rate
# Dependency:
#  genCorrMat,
#  cdBoundary,
#  pmvnorm (package mvtnorm),
#  pnorm (package stats),
#  bounds (package ldbounds).
#
#'
#' Difference between the Error Rate and Significance Level, the Error Spending Approach
#'
#' Calculate the difference between the error rate and significance level for the secondary endpoint, Lan-DeMets error spending approach.
#'
#' @param ap significance level for the primary endpoint
#' @param alpha targeted significance level for the secondary endpoint
#' @param cvec a vector of calculated primary boundaries
#' @param tVec a vector of information, gammaVec = sqrt(tVec)
#' @param ExtrmLoc an integer between 1 and K, locate the maximum of type I error of secondary endpoint
#' @param type type of sequential procedures. Type 1 OBF d, Type 2 POC d.
#' @return difference between alpha and the calculated error rate.
#'
#' @export
#' @import mvtnorm
#' @import stats
#' @import ldbounds
#'
#' @author Jiangtao Gou
#'
#' @seealso \code{secControl}
#'
#' @references
#' Lan, K. K. G., and Demets, D. L. (1983). Discrete sequential boundaries for clinical trials. \emph{Biometrika} \bold{70}, 659-663.
#'
#'  Tamhane, A. C., Gou, J., Jennison, C., Mehta, C. R., and Curto, T. (2018). A gatekeeping procedure to test a primary and a secondary endpoint in a group sequential design with multiple interim looks. \emph{Biometrics}, 74(1), 40-48.
#'
ldSecControl <- function (ap,alpha,cvec,tVec,ExtrmLoc,type=2) {
  K <- length(tVec);
  gammaVec <- sqrt(tVec);
  if (type == 1) {
    dvec.obf <- bounds(tVec,iuse=c(1),alpha=c(ap));
    dvec <- dvec.obf$upper.bounds;
  } else if (type == 2) {
    dvec.poc <- bounds(tVec,iuse=c(2),alpha=c(ap));
    dvec <- dvec.poc$upper.bounds;
  } else {
    print("Type 1 OBF d, Type 2 POC d.");
    return(K);
  }
  dlt <- (cvec[ExtrmLoc] - dvec[ExtrmLoc])/gammaVec[ExtrmLoc];
  ErrRate <- 1 - pnorm(max(cvec[1]-gammaVec[1]*dlt, dvec[1]));
  for (i in 2:K) {
    lowerB <- cdBoundary(cvec[1:i], dvec[1:i], gammaVec[1:i], dlt, upper=FALSE);
    upperB <- cdBoundary(cvec[1:i], dvec[1:i], gammaVec[1:i], dlt, upper=TRUE);
    meanV <- rep(0,i);
    corrM <- genCorrMat(gammaVec[1:i],1);
    ErrRate = ErrRate + pmvnorm(lowerB, upperB, meanV, corrM, algorithm=Miwa(steps=128));
  }
  return(ErrRate - alpha);
}
# End of function ldSecControl
#

#
# Function 5
# function primaryBoundary
# Original OBF and POC boundaries for primary endpoints
# Input:
#  gammaVec: square root of information (vector)
#  alpha: significance level for the primary endpoint
#  type: type 1: OBF, type 2: POC.
#  initIntvl: parameters for uniroot function.
# Output:
#  Original OBF and POC boundaries (primary endpoints) (a number, cK)
# Dependency:
#  genCorrMat,
#  cdBoundary,
#  pmvnorm (package mvtnorm),
#  uniroot (package stats).
#
#' Calculate Primary Boundaries, Standard Approach
#'
#' Primary boundaries calculation of standard (original) OBF and POC.
#'
#' @export
#' @import mvtnorm
#' @import stats
#' @import ldbounds
#'
#' @author Jiangtao Gou
#'
#' @seealso \code{ldPrimaryBoundary}
#'
#' @param gammaVec a vector of square root of information.
#' @param alpha significance level
#' @param type type of sequential procedure. OBF is 1, POC is 2.
#' @param initIntvl paramter for function uniroot (two numbers)
#' @return original OBF and POC boundaries (primary endpoints) (a number, c_(K)).
#'
#' @references
#'  Tamhane, A. C., Gou, J., Jennison, C., Mehta, C. R., and Curto, T. (2017+). A gatekeeping procedure to test a primary and a secondary endpoint in a group sequential design with multiple interim looks. \emph{Biometrics}, to appear.
#'
primaryBoundary <- function (gammaVec,alpha,type=1,initIntvl=c(1,4)) {
  K <- length(gammaVec);
  if (type == 1) {
    target <- function (c,gammaVec,alpha) {
      lowerB <- rep(-Inf,K);
      upperB <- c/gammaVec;
      meanV <- rep(0,K);
      corrM <- genCorrMat(gammaVec,1);
      result <- pmvnorm(lowerB, upperB, meanV, corrM, algorithm=Miwa(steps=128));
      return(result - 1 + alpha);
    }
    result.OF <- uniroot(target, lower = initIntvl[1], upper = initIntvl[2], tol = 2.5e-16, gammaVec = gammaVec, alpha = alpha);
    return(result.OF$root);
  } else if (type == 2) {
    target <- function (c,gammaVec,alpha) {
      lowerB <- rep(-Inf,K);
      upperB <- c;
      meanV <- rep(0,K);
      corrM <- genCorrMat(gammaVec,1);
      result <- pmvnorm(lowerB, upperB, meanV, corrM, algorithm=Miwa(steps=128));
      return(result - 1 + alpha);
    }
    result.OF <- uniroot(target, lower = initIntvl[1], upper = initIntvl[2], tol = 2.5e-16, gammaVec = gammaVec, alpha = alpha);
    return(result.OF$root);
  } else {
    print("Only OBF (1) and POC (2) are currently available.");
    return(K);
  }
}
# End of function primaryBoundary
#

#
# Function 6
# function secControl
# Secondary boundaries calculation of original OBF and POC
# Input:
#  d: boundary of secondary endpoint (a number, dK)
#  gammaVec: square root of information (vector)
#  alpha: targeted significance level for the secondary endpoint
#  cvec: boundaries of primary endpoints
#  gammaVec: square root of information
#  ExtrmLoc: the location of the maximum of type I error of the secondary endpoint
#  type: Type 1 OBF d, Type 2 POC d.
# Output:
#  Original OBF and POC boundaries (primary endpoints) (a number, cK)
# Dependency:
#  genCorrMat,
#  cdBoundary,
#  pmvnorm (package mvtnorm).
#'
#' Difference between the Error Rate and Significance Level, Standard Approach
#'
#' Calculate the difference between the error rate and significance level for the secondary endpoint, standard (original) approach.
#'
#' @param d boundary of secondary endpoint at the final look (a number, d_(K))
#' @param alpha targeted significance level for the secondary endpoint
#' @param cvec a vector of calculated primary boundaries
#' @param gammaVec square root of information
#' @param ExtrmLoc an integer between 1 and K, locate the maximum of type I error of secondary endpoint
#' @param type type of sequential procedures. Type 1 OBF d, Type 2 POC d.
#' @return difference between alpha and the calculated error rate.
#'
#' @export
#' @import mvtnorm
#' @import stats
#' @import ldbounds
#'
#' @author Jiangtao Gou
#'
#' @seealso \code{ldSecControl}
#'
#' @references
#'  Tamhane, A. C., Gou, J., Jennison, C., Mehta, C. R., and Curto, T. (2018). A gatekeeping procedure to test a primary and a secondary endpoint in a group sequential design with multiple interim looks. \emph{Biometrics}, 74(1), 40-48.
#'
secControl <- function (d,alpha,cvec,gammaVec,ExtrmLoc,type=2) {
  K <- length(gammaVec);
  if (type == 1) {
    dvec = d/gammaVec;
  } else if (type == 2) {
    dvec = d/rep(1,K);
  } else {
    print("Type 1 OBF d, Type 2 POC d.");
    return(K);
  }
  dlt <- (cvec[ExtrmLoc] - dvec[ExtrmLoc])/gammaVec[ExtrmLoc];
  ErrRate <- 1 - pnorm(max(cvec[1]-gammaVec[1]*dlt, dvec[1]));
  for (i in 2:K) {
    lowerB <- cdBoundary(cvec[1:i], dvec[1:i], gammaVec[1:i], dlt, upper=FALSE);
    upperB <- cdBoundary(cvec[1:i], dvec[1:i], gammaVec[1:i], dlt, upper=TRUE);
    meanV <- rep(0,i);
    corrM <- genCorrMat(gammaVec[1:i],1);
    ErrRate = ErrRate + pmvnorm(lowerB, upperB, meanV, corrM, algorithm=Miwa(steps=128));
  }
  return(ErrRate - alpha);
}
# End of function secControl
#

#
# Function 7
# function primaryBoundaryVec
# Calculate the primary boundaries
# Input:
#  alpha: significance level for the primary endpoint
#  tVec: information (vector)
#  initIntvl: parameter for function uniroot (two numbers)
#     for function primaryBoundary or function ldPrimaryBoundary
#  OBF: TRUE for OBF, FALSE for POC
#  LanDeMets: TRUE for Lan-Demets type boundaries, FALSE for original boundaries
#  digits: number of digits for output,
#  printOut: TRUE for printing the boundaries.
# Output:
#  OBF and POC boundaries (primary endpoints) (vector)
# Dependency:
#  genCorrMat,
#  cdBoundary,
#  primaryBoundary,
#  ldPrimaryBoundary,
#  pmvnorm (package mvtnorm).
#
#'
#' Calculate the Primary Boundaries
#'
#' Primary boundaries are calculated, including the standard approach and the error spending approach.
#'
#' @param alpha significance level for the primary endpoint.
#' @param tVec information (vector).
#' @param initIntvl parameter for function uniroot (two numbers) for function primaryBoundary or function ldPrimaryBoundary
#' @param OBF type of procedures. \code{TRUE} for OBF, \code{FALSE} for POC.
#' @param LanDeMets type of procedures. \code{TRUE} for Lan-Demets type boundaries, \code{FALSE} for original boundaries.
#' @param digits number of digits for output,
#' @param printOut \code{TRUE} for printing the boundaries.
#' @return OBF and POC boundaries (primary endpoints) (vector).
#'
#' @author Jiangtao Gou
#'
#' @examples
#' #require(mvtnorm)
#' #K = 4
#' #alpha = 0.025
#' #tVec = (1:K)/K
#' #boundaryVector <- primaryBoundaryVec(alpha,tVec,initIntvl=c(1,4),
#' #   OBF=TRUE,LanDeMets=FALSE,digits=3,printOut=TRUE)
#' #boundaryVector <- primaryBoundaryVec(alpha,tVec,initIntvl=c(1,4),
#' #   OBF=FALSE,LanDeMets=FALSE,digits=3,printOut=TRUE)
#' #boundaryVector <- primaryBoundaryVec(alpha,tVec,initIntvl=c(1,8),
#' #   OBF=TRUE,LanDeMets=TRUE,digits=3,printOut=TRUE)
#' #boundaryVector <- primaryBoundaryVec(alpha,tVec,initIntvl=c(1,4),
#' #   OBF=FALSE,LanDeMets=TRUE,digits=3,printOut=TRUE)
#'
#' @references
#'  Jennison, C. and Turnbull, B. W. (2000). \emph{Group Sequential Methods with Applications to Clinical Trials}. Chapman and Hall/CRC, New York.
#'
#'  Lan, K. K. G., and Demets, D. L. (1983). Discrete sequential boundaries for clinical trials. \emph{Biometrika} \bold{70}, 659-663.
#'
#'  O'Brien, P. C., and Fleming, T. R. (1979). A multiple testing procedure for clinical trials. \emph{Biometrics} \bold{35}, 549-556.
#'
#'  Pocock, S. J. (1977). Group sequential methods in the design and analysis of clinical trials. \emph{Biometrika} \bold{64}, 191-199.
#'
#'  Tamhane, A. C., Gou, J., Jennison, C., Mehta, C. R., and Curto, T. (2018). A gatekeeping procedure to test a primary and a secondary endpoint in a group sequential design with multiple interim looks. \emph{Biometrics}, 74(1), 40-48.
#'
#' @export
#' @import mvtnorm
#' @import stats
#' @import ldbounds
#'
primaryBoundaryVec  <- function (alpha,tVec,OBF=TRUE,LanDeMets=FALSE,digits=2,printOut=TRUE,initIntvl=c(1,8)) {
  K <- length(tVec);
  gammaVec <- sqrt(tVec);
  if (OBF == TRUE) {
    typeOP = 1;
  } else {
    typeOP = 2;
  }
  if (LanDeMets == FALSE) {
    bndyRep <- rep(0,1);
    for (i in 1:1) {
      bndyRep[i] <- primaryBoundary(gammaVec,alpha,typeOP,initIntvl);
    }
    cK <- median(bndyRep);
    if (OBF == TRUE) {
      bndyVec <- rep(cK,K)/gammaVec;
    } else {
      bndyVec <- rep(cK,K);
    }
    rBndyVec <- round(bndyVec, digits);
    if (printOut == TRUE) {
      print(rBndyVec);
    }
    return(bndyVec);
  } else {  # Lan-DeMets
    bndyRep <- matrix(rep(0,1*K),nrow=1);
    for (i in 1:1) {
      bndyRep[i,] <- ldPrimaryBoundary(tVec,alpha,typeOP,initIntvl);
    }
    bndyVec <- rep(0,K);
    for (j in 1:K) {
      bndyVec[j] <- median(bndyRep[,j]);
    }
    rBndyVec <- round(bndyVec, digits);
    if (printOut == TRUE) {
      print(rBndyVec);
    }
    return(bndyVec);
  }
}
#
# End of function primaryBoundaryVec
#

#
# Function 8
# function initLocPeak
# Return the initial location of maximum, with originical d-vec
# Input:
#  alpha: significance level for the primary endpoint
#  tVec: information (vector)
#  cVec: primary boundaries are given
#  type: OBF 1, POC 2
#  initIntvl: parameters for uniroot
# Output:
#  Inital location of maximum (a number between 1 and K)
# Dependency:
#  secControl
#  uniroot (package stats).
#
#' Find the Location of Maximum, Standard OBF and POC
#'
#' Calculate the location of maximal tyep I error of the standard O'Brien-Fleming and Pocock refined secondary boundaries.
#'
#' @author Jiangtao Gou
#' @param alpha type I error.
#' @param tVec information vector.
#' @param cvec primary group sequential boundary.
#' @param type type of the test procedure for the secondary endpoint. O'Brien- Fleming (OBF) type error spending funciton is 1, Pocock (POC) type error spending funciton is 2.
#' @param initIntvl computing paramter, a pair of numbers containing the end-points of the interval to be searched for the root.
#' @return location of maximum, a number between 1 and the number of interims
#'
#' @details
#' This function search the location of the maximal point, in order to calculate the standard (origiinal) O'Brien-Fleming (OBF) and Pocock (POC) refined secondary boundaries.
#'
#' @export
#' @import mvtnorm
#' @import stats
#' @import ldbounds
#'
#' @seealso \code{SecondaryBoundary}, \code{ldInitLocBeak}
#' @examples
#' #require(mvtnorm)
#' #K <- 8
#' #gammaVec <- sqrt((1:K)/K)
#' #tVec <- gammaVec^2
#' #alpha = 0.025
#' #c <- 2.072274
#' #cvec <- c/gammaVec
#' #loc <- initLocPeak(alpha,tVec,cvec,type=2,initIntvl=c(1,3))
#'
#' @references
#'  O'Brien, P. C., and Fleming, T. R. (1979). A multiple testing procedure for clinical trials. \emph{Biometrics} \bold{35}, 549-556.
#'
#'  Pocock, S. J. (1977). Group sequential methods in the design and analysis of clinical trials. \emph{Biometrika} \bold{64}, 191-199.
#'
#'  Tamhane, A. C., Gou, J., Jennison, C., Mehta, C. R., and Curto, T. (2018). A gatekeeping procedure to test a primary and a secondary endpoint in a group sequential design with multiple interim looks. \emph{Biometrics}, to appear.
#'
initLocPeak <- function (alpha,tVec,cvec,type=2,initIntvl=c(1,4)) {
  K <- length(tVec);
  gammaVec <- sqrt(tVec);
  locVec <- rep(0,1);
  for (j in 1:1) {
    resultVec <- rep(0,K); # Need to be predefined, otherwise the program does not work
    for (i in 1:K) {
      result <- uniroot(secControl,lower=initIntvl[1],upper=initIntvl[2],tol=2.5e-16,alpha=alpha,cvec=cvec,gammaVec=gammaVec,ExtrmLoc=i,type=type);
      resultVec[i] <- result$root;
      # print(resultVec);
    }
    locVec[j] <- which.max(resultVec);
  }
  # print(locVec);
  locPeak <- as.numeric(names(sort(-table(locVec)))[1]) # Find the mode
  # Reference: <http://stackoverflow.com/questions/2547402/is-there-a-built-in-function-for-finding-the-mode>
  return(locPeak);
}
# End of function initLocPeak
#

#
# Function 9
# function secondaryBoundary
# Calculate the secondary boundaries (original OBF POC)
# Input:
#  alpha: significance level for the primary endpoint
#  tVec: information (vector)
#  cVec: primary boundaries are given
#  locPeak: location of maximum, returned by function initLocPeak
#  type: OBF 1, POC 2
#  initIntvl: parameters for uniroot
# Output:
#  boundary vector
# Dependency:
#  secControl
#  uniroot (package stats).
#
#
#' Calculate the Refined Secondary Boundaries, Standard OBF and POC
#'
#' Calculate the standard O'Brien-Fleming and Pocock refined secondary boundaries
#'
#' @author Jiangtao Gou
#' @param alpha type I error.
#' @param tVec information vector.
#' @param cvec primary group sequential boundary.
#' @param locPeak location of maximum, a number between 1 and the number of interims.
#' @param type type of the test procedure for the secondary endpoint. O'Brien- Fleming (OBF) type error spending funciton is 1, Pocock (POC) type error spending funciton is 2.
#' @param initIntvl computing paramter, a pair of numbers containing the end-points of the interval to be searched for the root.
#' @return standard O'Brien-Fleming and Pocock refined secondary boundaries.
#'
#' @details
#' This function calculates the standard (origiinal) O'Brien-Fleming (OBF) and Pocock (POC) refined secondary boundaries.
#'
#' @export
#' @import mvtnorm
#' @import stats
#' @import ldbounds
#'
#' @seealso \code{ldSecondaryBoundary}, \code{initLocBeak}
#' @examples
#' #require(mvtnorm)
#' #K <- 8
#' #gammaVec <- sqrt((1:K)/K)
#' #tVec <- gammaVec^2
#' #alpha = 0.025
#' #c <- 2.072274
#' #cvec <- c/gammaVec
#' #loc <- initLocPeak(alpha,tVec,cvec,type=2,initIntvl=c(1,4))
#' #sbvec <- secondaryBoundary(alpha,tVec,cvec,loc,type=2,
#' #       initIntvl=c(1,8))
#'
#' @references
#'  O'Brien, P. C., and Fleming, T. R. (1979). A multiple testing procedure for clinical trials. \emph{Biometrics} \bold{35}, 549-556.
#'
#'  Pocock, S. J. (1977). Group sequential methods in the design and analysis of clinical trials. \emph{Biometrika} \bold{64}, 191-199.
#'
#'  Tamhane, A. C., Gou, J., Jennison, C., Mehta, C. R., and Curto, T. (2017+). A gatekeeping procedure to test a primary and a secondary endpoint in a group sequential design with multiple interim looks. \emph{Biometrics}, 74, 40-48.
#'
secondaryBoundary <- function (alpha,tVec,cvec,locPeak,type=2,initIntvl=c(1,4)) {
  K <- length(tVec);
  gammaVec <- sqrt(tVec);
  resultVec <- rep(0,1);
  for (i in 1:1) {
    result <- uniroot(secControl, lower = initIntvl[1], upper = initIntvl[2], tol = 2.5e-16, alpha=alpha,cvec=cvec,gammaVec=gammaVec,ExtrmLoc=locPeak,type=type);
    resultVec[i] <- result$root;
  }
  dK <- median(resultVec) #
  if (type == 1) {
    bndyVec <- rep(dK,K)/gammaVec;
  } else if (type == 2) {
    bndyVec <- rep(dK,K);
  } else {
    print("Type 1 OBF, Type 2 POC.");
  }
  return(bndyVec);
}
# End of function secondaryBoundary
#

#
# Function 10
# function ldInitLocPeak
# Return the initial location of maximum, with originical d-vec, for Lan-DeMets
# Input:
#  alpha: significance level for the primary endpoint
#  tVec: information (vector)
#  cVec: primary boundaries are given
#  type: OBF 1, POC 2
#  initIntvl: parameters for uniroot
# Output:
#  Inital location of maximum (a number between 1 and K)
# Dependency:
#  ldSecControl
#  uniroot (package stats).
#
#
#' Find the Location of Maximum, Error Spending Approach
#'
#' Calculate the location of maximal type I error of secondary endpoint.
#' @author Jiangtao Gou
#' @param alpha type I error.
#' @param tVec information vector.
#' @param cvec primary group sequential boundary.
#' @param type type of the test procedure for the secondary endpoint. O'Brien- Fleming (OBF) type error spending funciton is 1, Pocock (POC) type error spending funciton is 2.
#' @param initIntvl computing paramter, a pair of numbers containing the end-points of the interval to be searched for the root.
#' @return location of maximum, a number between 1 and the number of interims.
#'
#' @details
#' This function searches the location of maximal type I error of secondary endpoint by using the error spending approach.
#'
#' @export
#' @import mvtnorm
#' @import stats
#' @import ldbounds
#'
#' @seealso \code{ldSecondaryBoundary}, \code{initLocBeak}
#' @examples
#' #require(mvtnorm)
#' #require(ldbounds)
#' #K <- 6;
#' #tVec <- c(140,328,453,578,659,1080)/1080;
#' #alpha = 0.025;
#' #cvec.obf <- bounds(tVec,iuse=c(1),alpha=c(alpha));
#' #cvec <- cvec.obf$upper.bounds;
#' #loc <- ldInitLocPeak(alpha,tVec,cvec,type=2,initIntvl=c(0.9,4))
#'
#' @references
#'  Lan, K. K. G., and Demets, D. L. (1983). Discrete sequential boundaries for clinical trials. \emph{Biometrika} \bold{70}, 659-663.
#'
#'  Tamhane, A. C., Gou, J., Jennison, C., Mehta, C. R., and Curto, T. (2018). A gatekeeping procedure to test a primary and a secondary endpoint in a group sequential design with multiple interim looks. \emph{Biometrics}, 74, 40-48.
#'
#'
ldInitLocPeak <- function (alpha,tVec,cvec,type=2,initIntvl=c(0.8,4)) {
  K <- length(tVec);
  gammaVec <- sqrt(tVec);
  locVec <- rep(0,1);
  for (j in 1:1) {
    resultVec <- rep(0,K); # Need to be predefined, otherwise the program does not work
    for (i in 1:K) {
      result <- uniroot(ldSecControl,lower=initIntvl[1]*alpha,upper=initIntvl[2]*alpha,tol=2.5e-16,alpha=alpha,cvec=cvec,tVec=tVec,ExtrmLoc=i,type=type);
      resultVec[i] <- result$root;
      # print(resultVec);
    }
    locVec[j] <- which.min(resultVec);
  }
  # print(locVec);
  locPeak <- as.numeric(names(sort(-table(locVec)))[1]) # Find the mode
  # Reference: <http://stackoverflow.com/questions/2547402/is-there-a-built-in-function-for-finding-the-mode>
  return(locPeak);
}
# End of function ldInitLocPeak
#

#
# Function 11
# function ldSecondaryBoundary
# Input:
#  alpha: significance level for the primary endpoint
#  tVec: information (vector)
#  cVec: primary boundaries are given
#  locPeak: location of maximum, returned by function initLocPeak
#  type: OBF 1, POC 2
#  initIntvl: parameters for uniroot
# Output:
#  boundary vector
# Dependency:
#  ldSecControl
#  uniroot (package stats).
#
#
#' Calculate Refined Secondary Boundary, Error Spending Approach
#'
#' Refined secondary boundaries are calculated by using the error spending approach.
#'
#' @author Jiangtao Gou
#' @param alpha original significance level.
#' @param tVec information vector.
#' @param cvec primary group sequential boundary.
#' @param locPeak location of maximum, a number between 1 and the number of interims.
#' @param type type of the test procedure for the secondary endpoint. O'Brien- Fleming (OBF) type error spending funciton is 1, Pocock (POC) type error spending funciton is 2.
#' @param initIntvl computing paramter, a pair of numbers containing the end-points of the interval to be searched for the root.
#' @return refined secondary boundaries.
#'
#' @details
#' This function calculates the refined secondary boundaries of any Lan-DeMets error spending boundary based on the primary boundaries.
#'
#' @export
#' @import mvtnorm
#' @import stats
#' @import ldbounds
#'
#' @seealso \code{secondaryBoundary}, \code{secondaryBoundaryVecLD}
#' @examples
#' #require(mvtnorm)
#' #require(ldbounds)
#' #K <- 6;
#' #tVec <- c(140,328,453,578,659,1080)/1080;
#' #alpha = 0.025;
#' #cvec.obf <- bounds(tVec,iuse=c(1),alpha=c(alpha));
#' #cvec <- cvec.obf$upper.bounds;
#' #secbound <- ldSecondaryBoundary(alpha,tVec,cvec,locPeak=4,type=2,
#' #    initIntvl=c(0.8,8))
#'
#' @references
#'  Lan, K. K. G., and Demets, D. L. (1983). Discrete sequential boundaries for clinical trials. \emph{Biometrika} \bold{70}, 659-663.
#'
#'  Tamhane, A. C., Gou, J., Jennison, C., Mehta, C. R., and Curto, T. (2017+). A gatekeeping procedure to test a primary and a secondary endpoint in a group sequential design with multiple interim looks. \emph{Biometrics}, 74, 40-48.
#'
ldSecondaryBoundary <- function (alpha,tVec,cvec,locPeak,type=2,initIntvl=c(0.6,4)) {
  K <- length(tVec);
  gammaVec <- sqrt(tVec);
  resultVec <- rep(0,1);
  for (i in 1:1) {
    result <- uniroot(ldSecControl, lower = initIntvl[1]*alpha, upper = initIntvl[2]*alpha, tol = 2.5e-16, alpha=alpha,cvec=cvec,tVec=tVec,ExtrmLoc=locPeak,type=type);
    resultVec[i] <- result$root;
  }
  aprime <- median(resultVec) #
  if (type == 1) {
    dvec.obf <- bounds(tVec,iuse=c(1),alpha=c(aprime));
    ldBndyVec <- dvec.obf$upper.bounds;
  } else if (type == 2) {
    dvec.poc <- bounds(tVec,iuse=c(2),alpha=c(aprime));
    ldBndyVec <- dvec.poc$upper.bounds;
  } else {
    print("Type 1 OBF, Type 2 POC.");
  }
  return(ldBndyVec);
}
# End of function ldSecondaryBoundary
#

#
# Function 12
# nominal significance
# Input:
#  alpha: significance level for the primary endpoint
#  tVec: information (vector)
#  cVec: primary boundaries are given
#  locPeak: location of maximum, returned by function initLocPeak
#  type: OBF 1, POC 2
#  initIntvl: parameters for uniroot
# Output:
#  nominal significance
# Dependency:
#  ldSecControl
#  uniroot (package stats).
#
#' Calculate Nominal Significance, Error Spending Approach
#'
#' Nominal significance for the secondary endpoint are calculated by using the error spending approach.
#'
#' @param alpha original significance level.
#' @param tVec information vector.
#' @param cvec primary group sequential boundary.
#' @param locPeak location of maximum, a number between 1 and the number of interims.
#' @param type O'Brien- Fleming (OBF) type error spending funciton is 1, Pocock (POC) type error spending funciton is 2.
#' @param initIntvl computing paramter, a pair of numbers containing the end-points of the interval to be searched for the root.
#' @return nominal significance of the secondary group sequential boundary.
#'
#' @details
#' This function calculates the nominal significance level of any Lan-DeMets error spending boundary.
#' The original significance level is used to choose the initial searching range of the nominal significance.
#'
#' @export
#' @import mvtnorm
#' @import stats
#' @import ldbounds
#'
#' @author Jiangtao Gou
#'
#' @seealso \code{nominalSig}, \code{secondaryBoundaryVecLD}
#' @examples
#' #require(mvtnorm)
#' #require(ldbounds)
#' #K <- 6;
#' #tVec <- c(140,328,453,578,659,1080)/1080;
#' #alpha = 0.025;
#' #cvec.obf <- bounds(tVec,iuse=c(1),alpha=c(alpha));
#' #cvec <- cvec.obf$upper.bounds;
#' #alphaprime <- ldNominalSig(alpha,tVec,cvec,locPeak=4,type=2,
#' #      initIntvl=c(1,4))
#'
#' @references
#'  Lan, K. K. G., and Demets, D. L. (1983). Discrete sequential boundaries for clinical trials. \emph{Biometrika} \bold{70}, 659-663.
#'
#'  Tamhane, A. C., Gou, J., Jennison, C., Mehta, C. R., and Curto, T. (2018). A gatekeeping procedure to test a primary and a secondary endpoint in a group sequential design with multiple interim looks. \emph{Biometrics}, 74, 40-48.
#'
ldNominalSig <- function (alpha,tVec,cvec,locPeak,type=2,initIntvl=c(1,4)) {
  K <- length(tVec);
  gammaVec <- sqrt(tVec);
  resultVec <- rep(0,1);
  for (i in 1:1) {
    result <- uniroot(ldSecControl, lower = initIntvl[1]*alpha, upper = initIntvl[2]*alpha, tol = 2.5e-16, alpha=alpha,cvec=cvec,tVec=tVec,ExtrmLoc=locPeak,type=type);
    resultVec[i] <- result$root;
  }
  aprime <- median(resultVec) #
  return(aprime)
}
# End of function ldNominalSig
#

#
# Function 13
# function nominalSig
# Find the nominal Significance level for Original OBF and POC boundaries
# Input:
#  gammaVec: square root of information (vector)
#  alpha: significance level for the primary endpoint
#  cvec: the boundaries
# Output:
#  Nominal significance level
# Dependency:
#  genCorrMat,
#  cdBoundary,
#  pmvnorm (package mvtnorm)
#
#' Calculate Nominal Significance, Standard Approach
#'
#' Nominal significance for the secondary endpoint are calculated by using the standard (original) approach.
#'
#' @param gammaVec square root of information.
#' @param cvec group sequential boundary.
#' @return nominal significance
#'
#' @details
#' This function calculates he nominal significance level of any given boundary.
#'
#' @export
#' @import mvtnorm
#' @import stats
#' @import ldbounds
#'
#' @seealso \code{ldNominalSig}, \code{secondaryBoundaryVecOrig}
#' @examples
#' #require(mvtnorm)
#' #require(ldbounds)
#' #nSig <- nominalSig(gammaVec=c(sqrt(1/3),1),cvec=c(2.2,1.8))
#'
#' @author Jiangtao Gou
#'
#' @references
#'  O'Brien, P. C., and Fleming, T. R. (1979). A multiple testing procedure for clinical trials. \emph{Biometrics} \bold{35}, 549-556.
#'
#'  Pocock, S. J. (1977). Group sequential methods in the design and analysis of clinical trials. \emph{Biometrika} \bold{64}, 191-199.
#'
#'  Tamhane, A. C., Gou, J., Jennison, C., Mehta, C. R., and Curto, T. (2018). A gatekeeping procedure to test a primary and a secondary endpoint in a group sequential design with multiple interim looks. \emph{Biometrics}, 74, 40-48.
#'
nominalSig <- function (gammaVec,cvec) {
  K <- length(gammaVec);
  lowerB <- rep(-Inf,K); #print(lowerB)
  upperB <- cvec; #print(upperB)
  meanV <- rep(0,K);
  corrM <- genCorrMat(gammaVec,1); #print(corrM)
  result <- pmvnorm(lowerB, upperB, meanV, corrM, algorithm=Miwa(steps=128));
  return(1 - result[1]);
}
# End of function primaryBoundary
#

#
# Function 14
# function secondaryBoundaryVecOrig
# Calculate the secondary boundaries
# Input:
#  alpha: significance level for the secondary endpoint
#  tVec: information (vector)
#  initIntvl: parameter for function uniroot (two numbers)
#     for function primaryBoundary or function ldPrimaryBoundary
#  primaryOBF: TRUE for OBF, FALSE for POC
#  secondaryOBF: TRUE for OBF, FALSE for POC
#  digits: number of digits for output,
#  printOut: TRUE for printing the boundaries.
# Output:
#  OBF or POC boundaries (secondary endpoints) (vector)
# Dependency:
#  genCorrMat,
#  cdBoundary,
#  primaryBoundary,
#  ldPrimaryBoundary,
#  pmvnorm (package mvtnorm).
#<http://stackoverflow.com/questions/8936099/returning-multiple-objects-in-an-r-function>
#
#
#' Calculate Refined Secondary Boundaries and Nominal Significance, Standard Approach
#'
#' Standard refined secondary boundaries, and nominal significance for the secondary endpoint are calculated by using the standard (original) approach.
#'
#' @author Jiangtao Gou
#' @param alpha type I error probability.
#' @param tVec vector of relative information levels. The last element in the vector is 1.
#' @param primaryOBF type of primary boundary, \code{TURE} is the O'Brien-Fleming boundary, \code{FALSE} is the Pocock boundary.
#' @param secondaryOBF type of secondary boundary, \code{TURE} is the O'Brien-Fleming boundary, \code{FALSE} is the Pocock boundary.
#' @param initIntvl computing paramter, a pair of numbers containing the end-points of the interval to be searched for the root.
#' @return a result list including standard refined secondary boundary and the nominal significance for the secondary endpoint.
#'
#' @details
#' This function uses the standard approach (O'Brien and Fleming 1979, Pocock 1977),
#' and gives a list including refined secondary boundary and the nominal significance for the secondary endpoint.
#' There is a computing parameter \code{initIntvl}. 
#' Parameter \code{initIntvl} contains the end-points of the interval to be searched for the root.
#' The lower end point should choose a number around 1,
#' and the upper end point should choose a number between 4 and 10.
#'
#'
#' @export
#' @import mvtnorm
#' @import stats
#' @import ldbounds
#'
#' @seealso \code{secondaryBoundaryVec}, \code{secondaryBoundaryVecLD}
#' @examples
#' #require(mvtnorm)
#' #require(ldbounds)
#' #result <- secondaryBoundaryVecOrig(alpha=0.025,tVec=c(1/2,1),primaryOBF=TRUE,
#' #        secondaryOBF=FALSE, initIntvl=c(1,4))
#' #result$secondaryBoundary
#' #result$nomialSignificance
#'
#' @references
#'  Glimm, E., Maurer, W., and Bretz, F. (2010). Hierarchical testing of multiple endpoints in group-sequential trials. \emph{Statistics in Medicine} \bold{29}, 219-228.
#'
#'  Hung, H. M. J., Wang, S.-J., and O'Neill, R. (2007). Statistical considerations for testing multiple endpoints in group sequential or adaptive clinical trials. \emph{Journal of Biopharmaceutical Statistics} \bold{17}, 1201-1210.
#'
#'  Jennison, C. and Turnbull, B. W. (2000). \emph{Group Sequential Methods with Applications to Clinical Trials}. Chapman and Hall/CRC, New York.
#'
#'  O'Brien, P. C., and Fleming, T. R. (1979). A multiple testing procedure for clinical trials. \emph{Biometrics} \bold{35}, 549-556.
#'
#'  Pocock, S. J. (1977). Group sequential methods in the design and analysis of clinical trials. \emph{Biometrika} \bold{64}, 191-199.
#'
#'  Tamhane, A. C., Mehta, C. R., and Liu, L. (2010). Testing a primary and a secondary endpoint in a group sequential design. \emph{Biometrics} \bold{66}, 1174-1184.
#'
#'  Tamhane, A. C., Gou, J., Jennison, C., Mehta, C. R., and Curto, T. (2018). A gatekeeping procedure to test a primary and a secondary endpoint in a group sequential design with multiple interim looks. \emph{Biometrics}, 74, 40-48. 
#'
secondaryBoundaryVecOrig  <- function (alpha,tVec,primaryOBF=TRUE,secondaryOBF=FALSE,initIntvl=c(1,8)) {
  #
  K <- length(tVec);
  gammaVec <- sqrt(tVec);
  if (primaryOBF == secondaryOBF) {
    bndyVec <- primaryBoundaryVec(alpha,tVec,OBF=primaryOBF,LanDeMets=FALSE,digits=10,printOut=FALSE,initIntvl=initIntvl);
    nomialSigLvl <- alpha;
  } else {
    primaryBndyVec <- primaryBoundaryVec(alpha=alpha,tVec=tVec,OBF=primaryOBF,LanDeMets=FALSE,digits=10,printOut=FALSE,initIntvl=initIntvl);
    if (secondaryOBF == TRUE) {
      sType = 1;
    } else {
      sType = 2;
    }
    peakLocation <- initLocPeak(alpha=alpha,tVec=tVec,cvec=primaryBndyVec,type=sType,initIntvl=initIntvl);
    bndyVec <- secondaryBoundary(alpha,tVec,cvec=primaryBndyVec,locPeak=peakLocation,type=sType,initIntvl=initIntvl);
    #print(bndyVec)
    #
    nomialSigLvl <- nominalSig(gammaVec=gammaVec,cvec=bndyVec);
  }
  resultlist <- list("secondaryBoundary" = bndyVec, "nomialSignificance" = nomialSigLvl)
  return(resultlist)
}
#
# End of function secondaryBoundaryVecOrig
#

#
# Function 15
# function secondaryBoundaryVecLD
#
#
#' Calculate Refined Secondary Boundaries and Nominal Significance, the Error Spending Approach
#'
#' Lan-DeMets refined secondary boundaries, and nominal significance for the secondary endpoint are calculated by using the error spending approach.
#'
#' @author Jiangtao Gou
#' @param alpha type I error probability.
#' @param tVec vector of relative information levels. The last element in the vector is 1.
#' @param primaryOBF type of primary boundary, \code{TURE} is the O'Brien-Fleming boundary, \code{FALSE} is the Pocock boundary.
#' @param secondaryOBF type of secondary boundary, \code{TURE} is the O'Brien-Fleming boundary, \code{FALSE} is the Pocock boundary.
#' @param initIntvl computing paramter, a pair of numbers containing the end-points of the interval to be searched for the root.
#' @return a result list including Lan-DeMets refined secondary boundary and the nominal significance for the secondary endpoint.
#'
#' @details
#' This function uses the Lan-DeMets error spending approach,
#' and gives a list including refined secondary boundary and the nominal significance for the secondary endpoint.
#' There is a computing parameter \code{initIntvl}. 
#' Parameter \code{initIntvl} contains the end-points of the interval to be searched for the root.
#' For Lan-DeMets error spending approach, the lower end point should choose a number slightly less than 1,
#' and the upper end point should choose a number between 4 and 10.
#'
#'
#' @export
#' @import mvtnorm
#' @import stats
#' @import ldbounds
#'
#' @seealso \code{secondaryBoundaryVec}, \code{secondaryBoundaryVecOrig}
#' @examples
#' #require(mvtnorm)
#' #require(ldbounds)
#' #result <- secondaryBoundaryVecLD(alpha=0.025,tVec=c(1/2,1),primaryOBF=TRUE,
#' #        secondaryOBF=FALSE,initIntvl=c(0.8,6))
#' #result$secondaryBoundary
#' #result$nomialSignificance
#'
#' @references
#'  Glimm, E., Maurer, W., and Bretz, F. (2010). Hierarchical testing of multiple endpoints in group-sequential trials. \emph{Statistics in Medicine} \bold{29}, 219-228.
#'
#'  Hung, H. M. J., Wang, S.-J., and O'Neill, R. (2007). Statistical considerations for testing multiple endpoints in group sequential or adaptive clinical trials. \emph{Journal of Biopharmaceutical Statistics} \bold{17}, 1201-1210.
#'
#'  Jennison, C. and Turnbull, B. W. (2000). \emph{Group Sequential Methods with Applications to Clinical Trials}. Chapman and Hall/CRC, New York.
#'
#'  Lan, K. K. G., and Demets, D. L. (1983). Discrete sequential boundaries for clinical trials. \emph{Biometrika} \bold{70}, 659-663.
#'
#'  O'Brien, P. C., and Fleming, T. R. (1979). A multiple testing procedure for clinical trials. \emph{Biometrics} \bold{35}, 549-556.
#'
#'  Pocock, S. J. (1977). Group sequential methods in the design and analysis of clinical trials. \emph{Biometrika} \bold{64}, 191-199.
#'
#'  Tamhane, A. C., Mehta, C. R., and Liu, L. (2010). Testing a primary and a secondary endpoint in a group sequential design. \emph{Biometrics} \bold{66}, 1174-1184.
#'
#'  Tamhane, A. C., Gou, J., Jennison, C., Mehta, C. R., and Curto, T. (2018). A gatekeeping procedure to test a primary and a secondary endpoint in a group sequential design with multiple interim looks. \emph{Biometrics}, 74, 40-48.
#'
secondaryBoundaryVecLD <- function (alpha,tVec,primaryOBF=TRUE,secondaryOBF=FALSE,initIntvl=c(0.8,8)) {
#
  K <- length(tVec);
  gammaVec <- sqrt(tVec);
  if (primaryOBF == secondaryOBF) {
    bndyVec <- primaryBoundaryVec(alpha,tVec,OBF=primaryOBF,LanDeMets=TRUE,digits=10,printOut=FALSE,initIntvl=initIntvl);
    nomialSigLvl <- alpha;
  } else {
    if (primaryOBF == TRUE) {
      pType = 1;
    } else {
      pType = 2;
    }
    if (secondaryOBF == TRUE) {
      sType = 1;
    } else {
      sType = 2;
    }
    cvec.bounds <- bounds(tVec,iuse=c(pType),alpha=c(alpha));
    cvec <- cvec.bounds$upper.bounds;
    peakLocation <- ldInitLocPeak(alpha=alpha,tVec=tVec,cvec=cvec,type=sType,initIntvl=initIntvl);
    bndyVec <- ldSecondaryBoundary(alpha=alpha,tVec=tVec,cvec=cvec,locPeak=peakLocation,type=sType,initIntvl=initIntvl);
    #
    nomialSigLvl <- ldNominalSig(alpha=alpha,tVec=tVec,cvec=cvec,locPeak=peakLocation,type=sType,initIntvl=initIntvl);
  }
  resultlist <- list("secondaryBoundary" = bndyVec, "nomialSignificance" = nomialSigLvl)
  return(resultlist)
}
# End of function secondaryBoundaryVecLD
#

#
# Function 16
#
# function secondaryBoundaryVec
#
#
#' Calculate Refined Secondary Boundaries and Nominal Significance
#'
#' Refined secondary boundaries, and nominal significance for the secondary endpoint are calculated.
#'
#' @author Jiangtao Gou
#' @param alpha type I error probability.
#' @param tVec vector of relative information levels. The last element in the vector is 1.
#' @param pOBF type of primary boundary, \code{TURE} is the O'Brien-Fleming boundary, \code{FALSE} is the Pocock boundary.
#' @param sOBF type of secondary boundary, \code{TURE} is the O'Brien-Fleming boundary, \code{FALSE} is the Pocock boundary.
#' @param LanDeMets type of boundary, \code{TRUE} is the error spending approach, \code{FALSE} is the original approach.
#' @param initIntvl computing paramter, a pair of numbers containing the end-points of the interval to be searched for the root.
#' @return a result list including refined secondary boundary and the nominal significance for the secondary endpoint.
#'
#' @details
#' This function gives a list including refined secondary boundary and the nominal significance for the secondary endpoint.
#' There are a computing parameter \code{initIntvl}. 
#' Parameter \code{initIntvl} contains the end-points of the interval to be searched for the root.
#' For Lan-DeMets error spending approach, the lower end point should choose a number slightly less than 1,
#' and the upper end point should choose a number between 4 and 10.
#'
#'
#' @export
#' @import mvtnorm
#' @import stats
#' @import ldbounds
#'
#' @seealso \code{secondaryBoundaryVecLD}, \code{secondaryBoundaryVecOrig}
#' @examples
#' #require(mvtnorm)
#' #require(ldbounds)
#' #result <- secondaryBoundaryVec(alpha=0.025,tVec=c(1/2,1),pOBF=TRUE,sOBF=FALSE,
#' #       LanDeMets=FALSE,initIntvl=c(0.8,5))
#' #result$secondaryBoundary
#' #result$nomialSignificance
#'
#' @references
#'  Glimm, E., Maurer, W., and Bretz, F. (2010). Hierarchical testing of multiple endpoints in group-sequential trials. \emph{Statistics in Medicine} \bold{29}, 219-228.
#'
#'  Hung, H. M. J., Wang, S.-J., and O'Neill, R. (2007). Statistical considerations for testing multiple endpoints in group sequential or adaptive clinical trials. \emph{Journal of Biopharmaceutical Statistics} \bold{17}, 1201-1210.
#'
#'  Jennison, C. and Turnbull, B. W. (2000). \emph{Group Sequential Methods with Applications to Clinical Trials}. Chapman and Hall/CRC, New York.
#'
#'  Lan, K. K. G., and Demets, D. L. (1983). Discrete sequential boundaries for clinical trials. \emph{Biometrika} \bold{70}, 659-663.
#'
#'  O'Brien, P. C., and Fleming, T. R. (1979). A multiple testing procedure for clinical trials. \emph{Biometrics} \bold{35}, 549-556.
#'
#'  Pocock, S. J. (1977). Group sequential methods in the design and analysis of clinical trials. \emph{Biometrika} \bold{64}, 191-199.
#'
#'  Tamhane, A. C., Mehta, C. R., and Liu, L. (2010). Testing a primary and a secondary endpoint in a group sequential design. \emph{Biometrics} \bold{66}, 1174-1184.
#'
#'  Tamhane, A. C., Gou, J., Jennison, C., Mehta, C. R., and Curto, T. (2018). A gatekeeping procedure to test a primary and a secondary endpoint in a group sequential design with multiple interim looks. \emph{Biometrics}, 74, 40-48. 
#'
secondaryBoundaryVec <- function (alpha,tVec,pOBF=TRUE,sOBF=FALSE,LanDeMets=FALSE,initIntvl=c(0.8,8)) {
  #
  if (LanDeMets == FALSE) {
    bdns <- secondaryBoundaryVecOrig(alpha=alpha,tVec=tVec,primaryOBF=pOBF,secondaryOBF=sOBF,initIntvl=initIntvl);
    # bdns: boundaries and nominal significance
  } else {
    bdns <- secondaryBoundaryVecLD(alpha=alpha,tVec=tVec,primaryOBF=pOBF,secondaryOBF=sOBF,initIntvl=initIntvl);
  }
  return(bdns);
}
# End of function secondaryBoundaryVec
#

#
# Function 17
#
# function refinedBoundary
#
#
#' Summarize Primary and Refined Secondary Boundaries, Nominal Significance
#'
#' Primary boundaries, refined secondary boundaries, and nominal significance for the secondary endpoint are listed.
#'
#' @author Jiangtao Gou
#' @param alpha type I error probability.
#' @param tVec vector of relative information levels. The last element in the vector is 1.
#' @param pOBF type of primary boundary, \code{TURE} is the O'Brien-Fleming boundary, \code{FALSE} is the Pocock boundary.
#' @param sOBF type of secondary boundary, \code{TURE} is the O'Brien-Fleming boundary, \code{FALSE} is the Pocock boundary.
#' @param LanDeMets type of boundary, \code{TRUE} is the error spending approach, \code{FALSE} is the original approach.
#' @param digits number of digits after decimal point for primary and secondary boundaries.
#' @return a result list including primary boundary, refined secondary boundary, and the nominal significance for the secondary endpoint.
#'
#' @details
#' This function gives a list including primary boundary, refined secondary boundary, and the nominal significance for the secondary endpoint.
#' The number of digits for the nominal significance depends on parameter \code{alpha}.
#'
#' @export
#' @import mvtnorm
#' @import stats
#' @import ldbounds
#'
#' @examples
#' require(mvtnorm)
#' require(ldbounds)
#' result <- refinedBoundary(alpha=0.05,tVec=c(0.2,0.6,1))
#' result$primaryBoundary
#' result$secondaryBoundary
#' result$nomialSignificance
#'
#' @references
#'  Glimm, E., Maurer, W., and Bretz, F. (2010). Hierarchical testing of multiple endpoints in group-sequential trials. \emph{Statistics in Medicine} \bold{29}, 219-228.
#'
#'  Hung, H. M. J., Wang, S.-J., and O'Neill, R. (2007). Statistical considerations for testing multiple endpoints in group sequential or adaptive clinical trials. \emph{Journal of Biopharmaceutical Statistics} \bold{17}, 1201-1210.
#'
#'  Jennison, C. and Turnbull, B. W. (2000). \emph{Group Sequential Methods with Applications to Clinical Trials}. Chapman and Hall/CRC, New York.
#'
#'  Lan, K. K. G., and Demets, D. L. (1983). Discrete sequential boundaries for clinical trials. \emph{Biometrika} \bold{70}, 659-663.
#'
#'  O'Brien, P. C., and Fleming, T. R. (1979). A multiple testing procedure for clinical trials. \emph{Biometrics} \bold{35}, 549-556.
#'
#'  Pocock, S. J. (1977). Group sequential methods in the design and analysis of clinical trials. \emph{Biometrika} \bold{64}, 191-199.
#'
#'  Tamhane, A. C., Mehta, C. R., and Liu, L. (2010). Testing a primary and a secondary endpoint in a group sequential design. \emph{Biometrics} \bold{66}, 1174-1184.
#'
#'  Tamhane, A. C., Gou, J., Jennison, C., Mehta, C. R., and Curto, T. (2018). A gatekeeping procedure to test a primary and a secondary endpoint in a group sequential design with multiple interim looks. \emph{Biometrics}, 74, 40-48
#'
refinedBoundary <- function (alpha,tVec,pOBF=TRUE,sOBF=FALSE,LanDeMets=FALSE,digits=2) {
  #
  initIntvl=c(0.8,8);
  pd <- primaryBoundaryVec(alpha=alpha,tVec=tVec,OBF=pOBF,LanDeMets=LanDeMets,digits=10,printOut=FALSE,initIntvl=initIntvl);
  bdns <- secondaryBoundaryVec(alpha=alpha,tVec=tVec,pOBF=pOBF,sOBF=sOBF,LanDeMets=LanDeMets,initIntvl=initIntvl);
  #
  rawlist <- list("primaryBoundary" = pd, "secondaryBoundary" = bdns$secondaryBoundary, "nomialSignificance" = bdns$nomialSignificance)
  dgt = digits;
  rpd <- round(pd, dgt);
  rsd <- round(bdns$secondaryBoundary, dgt);
  nsdgt <- 3+round(log10(1/alpha));
  nsdgt <- max(dgt, nsdgt); 
  rns <- round(bdns$nomialSignificance, nsdgt);
  resultlist <- list("primaryBoundary" = rpd, "secondaryBoundary" = rsd, "nomialSignificance" = rns)
  return(resultlist)
}
# End of function refinedBoundary
#

#
# Function 18
#
# function psbTeXtable
#
#' Summarize Primary and Refined Secondary Boundaries in a TeX table
#'
#' Primary boundaries and refined secondary boundaries are listed in a TeX table.
#'
#' @author Jiangtao Gou
#' @author Fengqing (Zoe) Zhang
#' @param alpha type I error probability.
#' @param tVec vector of relative information levels. The last element in the vector is 1.
#' @param pOBF type of primary boundary, \code{TURE} is the O'Brien-Fleming boundary, \code{FALSE} is the Pocock boundary.
#' @param sOBF type of secondary boundary, \code{TURE} is the O'Brien-Fleming boundary, \code{FALSE} is the Pocock boundary.
#' @param LanDeMets type of boundary, \code{TRUE} is the error spending approach, \code{FALSE} is the original approach.
#' @param digits number of digits after decimal point to display in the table.
#' @return a TeX format table including both primary boundary and refined secondary boundary.
#'
#' @details
#' This function gives a TeX format table including both primary boundary and refined secondary boundary. 
#' The number of digits after decimal point can be specified through parameter \code{digits}.
#'
#' @export
#' @import mvtnorm
#' @import stats
#' @import ldbounds
#' @import xtable
#'
#' @examples
#' #require(mvtnorm)
#' #require(ldbounds)
#' #require(xtable)
#' #psbTeXtable(alpha=0.025,tVec=c(1/2,3/4,1),pOBF=TRUE,sOBF=FALSE,LanDeMets=FALSE)
#' @references
#'  Glimm, E., Maurer, W., and Bretz, F. (2010). Hierarchical testing of multiple endpoints in group-sequential trials. \emph{Statistics in Medicine} \bold{29}, 219-228.
#'
#'  Hung, H. M. J., Wang, S.-J., and O'Neill, R. (2007). Statistical considerations for testing multiple endpoints in group sequential or adaptive clinical trials. \emph{Journal of Biopharmaceutical Statistics} \bold{17}, 1201-1210.
#'
#'  Jennison, C. and Turnbull, B. W. (2000). \emph{Group Sequential Methods with Applications to Clinical Trials}. Chapman and Hall/CRC, New York.
#'
#'  Lan, K. K. G., and Demets, D. L. (1983). Discrete sequential boundaries for clinical trials. \emph{Biometrika} \bold{70}, 659-663.
#'
#'  O'Brien, P. C., and Fleming, T. R. (1979). A multiple testing procedure for clinical trials. \emph{Biometrics} \bold{35}, 549-556.
#'
#'  Pocock, S. J. (1977). Group sequential methods in the design and analysis of clinical trials. \emph{Biometrika} \bold{64}, 191-199.
#'
#'  Tamhane, A. C., Mehta, C. R., and Liu, L. (2010). Testing a primary and a secondary endpoint in a group sequential design. \emph{Biometrics} \bold{66}, 1174-1184.
#'
#'  Tamhane, A. C., Gou, J., Jennison, C., Mehta, C. R., and Curto, T. (2018). A gatekeeping procedure to test a primary and a secondary endpoint in a group sequential design with multiple interim looks. \emph{Biometrics}, 74, 40-48.
#'
psbTeXtable <- function (alpha,tVec,pOBF=TRUE,sOBF=FALSE,LanDeMets=FALSE,digits=2) {
  #
  result <- refinedBoundary(alpha=alpha,tVec=tVec,pOBF=pOBF,sOBF=sOBF,LanDeMets=LanDeMets,digits=digits);
  psb <- data.frame(result$primaryBoundary, result$secondaryBoundary);
  colnames(psb) <- c("Primary Boundary", "Refined Secondary Boundary");
  print(xtable(psb,label = 'tab:psb',
               caption = 'Primary Boundary and Refined Secondary Boundary', digits=digits));
}
# End of function psbTeXtable
#

##########
