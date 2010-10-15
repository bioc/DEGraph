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
## @RdocFunction graph.T2.test
##
## @title "Performs the Hotelling T2 test in Fourier space"
##
## \description{
##  @get "title".
## }
##
## @synopsis
##
## \arguments{
##   \item{X1}{A n1 x p @numeric @matrix, observed data for class 1: p variables, n1
##     observations.}
##   \item{X2}{A n2 x p @numeric @matrix, observed data for class 2: p variables, n2
##     observations.}
##   \item{G}{An object of class \code{\link[=graphAM-class]{graphAM}} or
##     \code{\link[=graphNEL-class]{graphNEL}}, the graph to be used in the two-sample test.}
##   \item{lfA}{A list returned by @see "laplacianFromA",
##     containing the Laplacian eigen vectors and eigen values}
##   \item{...}{Further arguments to be passed to @see "laplacianFromA".}
##   \item{k}{A @numeric value, number of Fourier components retained for the
##     test.}
## }
##
## \value{
##  A @list with class "htest", as returned by @see "T2.test".
## }
##
## @author
##
## \seealso{
##   @see "T2.test"
##   \code{\link[=graphAM-class]{graphAM}}
## }
##
## @examples "../incl/graph.T2.test.Rex"
##
##*/########################################################################
graph.T2.test <- function(X1, X2, G=NULL, lfA=NULL, ..., k=ncol(X1))
{
  tol <- 1e-8
  p <- ncol(X1)
  
  if (is.null(lfA))
    {
      if(is.null(G))
        A <- diag(rep(1,p))
      else
        A <- as(G, "graphAM")@adjMat
      
      ## Build basis based on subgraph Laplacian
      lfA <- laplacianFromA(A, ..., k=k)            
    }    

  U <- lfA$U
  egVal <- lfA$l
  kIdx <- (egVal <= max(egVal[k],tol)) #lfA$kIdx
  rk <- max(which(kIdx)) ## "round" the number of kept eigenvectors to take into account eigenvalue multiplicity
  ##print(egVal)
  ##print(U[,1:rk, drop=FALSE])
  ##print(X1[1:3,])
  ##print((X1%*%U[, 1:rk, drop=FALSE]))
  ## Perform the test in the new basis
  ut <- T2.test(X1%*%U[, 1:rk, drop=FALSE], X2%*%U[, 1:rk, drop=FALSE])
}

############################################################################
## HISTORY
## 2010-09-23
## o Added an example.
############################################################################
