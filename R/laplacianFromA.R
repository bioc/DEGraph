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
## @RdocFunction laplacianFromA
##
## @title "Calculates the Laplacian associated to an adjacency matrix"
##
## \description{
##  @get "title".
## }
##
## @synopsis
##
## \arguments{
##   \item{A}{The adjacency matrix of the graph.}
##   \item{k}{...}
##   \item{ltype}{A @character value specifying the type of Laplacian to be
##     calculated.  Defaults to meanInfluence.}
## }
##
## \value{
##  A @list containing the following components:
##  \describe{
##    \item{U}{Eigenvectors of the graph Laplacian.}
##    \item{l}{Eigenvalues of the graph Laplacian}
##    \item{kIdx}{Multiplicity of '0' as eigenvalue.}
##   }
## }
##
## @author
##
## @examples "../incl/randomWAMGraph.Rex"
##
##*/########################################################################

## Types: unnormalized, normalized, meanInfluence, totalInfluence

laplacianFromA <- function(A, k=1, ltype=c("meanInfluence", "normalized", "unnormalized", "totalInfluence")) {
  ltype <- match.arg(ltype)
  tol <- 1e-8
  rownames(A) <- NULL
  colnames(A) <- NULL

  ltype <- match.arg(ltype)

  ##if (!isSymmetric.matrix(A)) {
  ##  print(A)
  ##  throw("Argument 'A' must be a symmetric matrix")
  ##}
  
  p <- nrow(A)
  I <- diag(rep(1,p))

  if(ltype %in% c("normalized","unnormalized")) # A must be made symmetric
    {
      tsIdx <- ((A == 0) & (t(A) != 0))
      A[tsIdx] <- t(A)[tsIdx]
    }
  
  if(ltype == "normalized")
    {
      iDs <- diag(1/sqrt(rowSums(abs(A))))
      L <- I - (iDs %*% A %*% iDs)
    }

  if(ltype == "unnormalized")
    {
      D <- diag(rowSums(abs(A)))
      L <- D - A
    }
              
  if(ltype == "meanInfluence")
    {
      ImA <- diag(as.integer(rowSums(abs(t(A))) != 0)) - diag(1/pmax(1,rowSums(abs(t(A)))))%*%t(A)
      L <- t(ImA)%*%ImA
    }

  if(ltype == "totalInfluence")
    {
      ImA <- diag(as.integer(rowSums(abs(t(A))) != 0)) - t(A)
      L <- t(ImA)%*%ImA
    }  

  edL <- eigen(L, symmetric=TRUE)
  egVal <- rev(edL$values)
  kIdx <- (egVal <= max(egVal[k], tol))
  return(list(U=edL$vectors[,ncol(edL$vectors):1], l=egVal, kIdx=kIdx))
}
