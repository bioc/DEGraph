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
## @RdocFunction BS.test
##
## @title "Performs the test of Bai and Saranadasa (1996)"
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
##   }
## }
##
## @author
##
## \seealso{
##   @see "AN.test"
##   @see "graph.T2.test"
##   @see "hyper.test"
## }
##
## @examples "../incl/tests.Rex"
##
##*/########################################################################

BS.test <- function(X1, X2, na.rm=FALSE) {
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
  n <- n1 + n2 - 2
  X1m <- colMeans(X1)
  X2m <- colMeans(X2)
  dm <- X1m - X2m
  
  dm2 <- t(dm) %*% dm
  
  X1c <- X1 - t(matrix(rep(colMeans(X1),n1),p,n1))
  X2c <- X2 - t(matrix(rep(colMeans(X2),n2),p,n2))  
  Sn <- ((t(X1c) %*% X1c) + (t(X2c) %*% X2c))/n
  
  Bn <- sqrt((n^2/((n+2)*(n-1))) * (sum(diag(Sn%*%Sn)) - (1/n)*(sum(diag(Sn)))^2))
  
  z <- ((n1*n2/(n1+n2))*dm2 - sum(diag(Sn)))/(sqrt((2*(n+1))/n)*Bn)
  p <- 1-pnorm(z)
  
  res <- list(statistic=z,  p.value=p)
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
############################################################################
