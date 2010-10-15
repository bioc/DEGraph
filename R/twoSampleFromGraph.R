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
## @RdocFunction twoSampleFromGraph
##
## @title "Given a basis (typically the eigenvectors of a graph Laplacian), builds two multivariate normal samples with mean shift located in the first elements of the basis"
##
## \description{
##  @get "title".
## }
##
## @synopsis
##
## \arguments{
##   \item{n1}{An @integer value specifying the number of points in the first sample.}
##   \item{n2}{An @integer value specifying the number of points in the second sample.}
##   \item{shiftM2}{A @numeric value giving the desired squared Mahalanobis norm of the mean shift between the two samples.}
##   \item{sigma}{A matrix giving the covariance structure of each sample.}
##   \item{U}{A matrix giving the desired basis.}
##   \item{k}{An @integer value giving the number of basis elements in which the mean shift must be located.}
## }
##
## \value{ A @list with named elements:
##  \describe{
##    \item{X1}{The first sample in the original basis (before transformation by U).}
##    \item{X2}{The second sample in the original basis (before transformation by U).}
##    \item{X1}{The first sample in the specified basis (after transformation by U).}
##    \item{X2}{The second sample in the specified basis (after transformation by U).}
##    \item{mu1}{The population mean of F1}
##    \item{mu2}{The population mean of F2}
##    \item{diff}{mu1 - mu2}
##   }
##  }
##
## @author
##
## @examples "../incl/randomWAMGraph.Rex"
##
##*/########################################################################

twoSampleFromGraph <-
function(n1=20, n2=n1, shiftM2=0, sigma, U, k=ceiling(ncol(U)/3)){
  p = nrow(U)

  ## Generate two samples in Fourier space

  diff <- rep(0,p)

  diff[1:k] <- rnorm(k)

  ## tmp <- rnorm(k)
  ## diff[1:k] <- shift*tmp/sqrt((t(tmp)%*%tmp))
  rawShiftNorm <- t(diff)%*%solve(sigma)%*%diff
  diff[1:k] <- sqrt(shiftM2)*diff[1:k]/sqrt(rawShiftNorm)  ## so that the norm of this shift is shiftM2 !

  mu1 <- rnorm(p)
  mu2 <- mu1 + diff

  F1 <- rmvnorm(n=n1, mean=mu1, sigma=sigma, method="chol")
  F2 <- rmvnorm(n=n2, mean=mu2, sigma=sigma, method="chol")

  ## Build the inverse Fourier transform from A
  ## iDs <- diag(1/sqrt(rowSums(abs(A))))
  ## L <- diag(rep(1,p)) - (iDs %*% A %*% iDs) # Normalized version
  ## D <- diag(rowSums(abs(A)))
  ## L <- D - A # Unormalized version
  ## edL <- eigen(L)
  ## U <- edL$vectors

  ## Return the data in the graph space

  X1 <- F1 %*% t(U)
  X2<- F2 %*% t(U)

  list(X1=X1,X2=X2,F1=F1,F2=F2,mu1=mu1,mu2=mu2,diff=diff)

}

