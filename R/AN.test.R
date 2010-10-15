## Copyright 2010 Laurent Jacob, Pierre Neuvial and Sandrine Dudoit.

## This file is part of DEGraph.

## DEGraph is free software: you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.

## DEGraph is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.

## You should have received a copy of the GNU General Public License
## along with DEGraph.  If not, see <http://www.gnu.org/licenses/>.

#########################################################################/**
## @RdocFunction AN.test
##
## @title "Performs the Adaptive Neyman test of Fan and Lin (1998)"
##
## \description{
##  @get "title".
## }
##
## @synopsis
##
## \arguments{
##   \item{X1}{A n1 x p @matrix, observed data for class 1: p variables, n1
##     observations.}
##   \item{X2}{A n2 x p @matrix, observed data for class 2: p variables, n2
##     observations.}
##   \item{candK}{A @vector, candidate values for the true number of Fourier
##     components.}
##   \item{na.rm}{A @logical value indicating whether variables with @NA in
##     at least one of the n1 + n2 observations should be discarder before
##     the test is performed.}
## }
##
## \value{
##  A @list with class "htest" containing the following components:
##  \describe{
##    \item{statistic}{A @numeric value, the test statistic.}
##    \item{p.value}{A @numeric value, the corresponding p-value.}
##    \item{kstar}{A @numeric value, the estimated true number of Fourier
##      components.}
##   }
## }
##
## @author
##
## \seealso{
##   @see "BS.test"
##   @see "graph.T2.test"
##   @see "hyper.test"
## }
##
## @examples "../incl/tests.Rex"
##
##*/########################################################################

AN.test <- function(X1, X2, candK=1:ncol(X1), na.rm=FALSE) {
  if (na.rm) {
    na1 <- apply(X1, 2, FUN=function(x) sum(is.na(x)))
    na2 <- apply(X2, 2, FUN=function(x) sum(is.na(x)))
    idxs <- which((na1==0) & (na2==0))
    X1 <- X1[, idxs]
    X2 <- X2[, idxs]    
  }
  p <- ncol(X1)
  n1 <- nrow(X1)
  n2 <- nrow(X2)
  Tstar <- -Inf
  kstar <- NA
  for (kk in candK) {
    var1 <- mean(diag(var(X1[, 1:kk, drop=FALSE])))
    var2 <- mean(diag(var(X2[, 1:kk, drop=FALSE])))
    X <- colMeans(X1[, 1:kk, drop=FALSE]) - colMeans(X2[, 1:kk, drop=FALSE])
    X <- X/sqrt(var1/n1 + var2/n2)
    tmp <- sum(X^2 - 1)/sqrt(2*kk)
    if (tmp > Tstar) {
      Tstar <- tmp
      kstar <- kk
    }
  }
  T <- sqrt(2*log(log(p)))*Tstar - (2*log(log(p)) + log(log(log(p)))/2 - log(4*pi)/2)
  p <- 1-exp(-exp(-T)) ## Gumbel cdf
  res <- list(statistic=T,  p.value=p, kstar=kstar)
  class(res) <- "htest"
  res
}

############################################################################
## HISTORY
## 2010-09-23
## o Added an example.
## o Added an option to remove NA:s.
## 2010-09-22
## o Now returns an object of class "htest".
## o Gumbel cdf is now calculated manually to avoid dependency 'evd'.
############################################################################
