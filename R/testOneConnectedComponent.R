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
## @RdocFunction testOneConnectedComponent
##
## @title "Applies a series of two-sample tests to a connected graph using various statistics"
##
## \description{
##  @get "title".
## }
##
## @synopsis
##
## \arguments{
##   \item{graph}{A \code{\link[=graph-class]{graph}} object.}
##   \item{data}{A '@numeric @matrix (size: number 'p' of genes x number 'n' of samples) of gene expression.}
##   \item{classes}{A @character @vector (length: 'n') of class assignments.}
##   \item{...}{Further arguments to be passed to @see "laplacianFromA".}
##   \item{prop}{A @numeric value, percentage of components retained for Fourier and PCA.}
##   \item{verbose}{If @TRUE, extra information is output.}
## }
##
## \value{
##  A structured @list containing the p-values of the tests, the 
##  \code{\link[=graph-class]{graph}} object of the connected component and the number 
##  of retained Fourier dimensions.
## }
##
## @author
##
## \details{
##   This function performs the test, assuming that all genes
##   in the graph are represented in the expression data set,
##   in order not to have to modify the graph topology.
##
##   Interaction signs are used if available in the graph
##   ('getSignedGraph' is not called here, in order not to
##   have to modify the graph topology.).
##
##   The graph given as input has to have only one
##   connex component.  It can be retrieved from the output of
##   @see "getConnectedComponentList".
## }
##
## \seealso{
##   @see "testOneGraph"
##   @see "getConnectedComponentList"
## }
##
## @examples "../incl/graph.T2.test.Rex"
##
##*/########################################################################

testOneConnectedComponent <- function(graph, data, classes, ..., prop=0.2, verbose=FALSE) {
  ## - - - - - - - - - - - - - - - - - - -
  ## Validate arguments
  ## - - - - - - - - - - - - - - - - - - -
  ## Argument 'classes'
  classes <- Arguments$getNumerics(classes)
  ll <- length(unique(classes))
  if (ll != 2) {
    throw("Number of classes should be 2, not ", ll)
  }
  rm(ll)

  ## Argument 'data'
  cls <- class(data)[1]
  if (cls != "matrix") {
    throw("Arguments 'data' should be a 'matrix', not a ", cls)
  }
  n <- ncol(data)
  nc <- length(classes)
  if (n != nc) {
    throw("Number of samples in 'data' (", n, ") should match number of samples in 'classes' (", nc, ")")
  }
  rm(nc)

  ## Argument 'graph'
  if (!validGraph(graph)) {
    throw("Argument 'graph' is not a valid graph object")
  }
  ## Check that there is only one connex component
  cc <- connectedComp(graph)
  if (length(cc)>1) {
    throw("More than one connex component in graph")
  }
  rm(cc)

  ## Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose)
  if (verbose) {
    cat <- R.utils::cat
    pushState(verbose)
    on.exit(popState(verbose))
  } 

  ## check consistency between gene and node names
  verbose && enter(verbose, "Keeping genes in the graph *and* the expression data set")
  dataGN <- rownames(data)
  graphGN <- nodes(graph)
  mm <- match(graphGN, dataGN)
  ww <- which(is.na(mm))
  ll <- length(ww)
  if (ll) {
    throw(ll, " genes of the graph were not found in the expression data set: ")
  }
  ## sorting (and subsetting if necessary)
  data <- data[mm, ]

  ## sanity check (ordering should be the same now)
  dataGN <- rownames(data)
  stopifnot(all.equal(dataGN, graphGN))
  rm(graphGN)
  rm(dataGN)
  verbose && exit(verbose)

  p <- numNodes(graph)
  geneNames <- rownames(data)

  cls <- sort(unique(classes))
  X1 <- t(data[, classes==cls[1]])
  X2 <- t(data[, classes==cls[2]])

  ## get sign information (if any) and infer from its presence the type of graph
  #signMat <- attr(graph, 'signMat')
  signMat <- graph@graphData$signMat
  if (is.null(signMat)) {
    verbose && cat(verbose, "Unsigned graph")
    ## use adjacency matrix of the corresponding undirected graph
    ##ugraph <- ugraph(graph)
    ##graphAM <- as(ugraph, "graphAM")
    ##adjMat <- graphAM@adjMat
    adjMat <- as(graph,"graphAM")@adjMat
    ## sanity check: symmetry
    ##stopifnot(sum(sum(t(adjMat)!=adjMat))==0)
    mat <- adjMat
  } else {
    verbose && cat(verbose, "Signed graph")
    mat <- signMat
  }

  ##print(mat)
  
  ## Argument 'prop'
  if (is.null(prop)) {
    prop <- c(0.05, 0.1, 0.2, 0.4)
    ndim <- unique(c(1:5, ceiling(p*prop)))
    ndim <- sort(ndim[ndim<=p])
  } else {
    prop <- Arguments$getNumerics(prop, range=c(0, 1))
    ndim <- unique(ceiling(p*prop))
  }

  dg <- rowSums(abs(mat))
  verbose && cat(verbose, "Degrees")
  verbose && str(verbose, dg)
  ## TODO: take multiplicity of eigen value into account using parameter 'k'
  lfA <- laplacianFromA(mat, ..., ltype="meanInfluence", k=1)
  U <- lfA$U

  res <- NULL
  resNames <- NULL

  ## T2 in original space
  verbose && enter(verbose, "Testing")
  verbose && cat(verbose, "Raw T2 (Hotelling)")
  T2H <- try(T2.test(X1, X2))
  if (class(T2H)=="try-error") {
    pH <- NA
  } else {
    pH <- T2H$p.value
    verbose & str(verbose, pH)
  }
  res <- c(res, pH)
  resNames <- c(resNames, "T2")

  ## T2 in Fourier space
  verbose && cat(verbose, "T2 on Fourier components")
  pU <- rep(NA, length=length(ndim))
  names(pU) <- ndim
  for (kk in seq(along=ndim)) {
    kR <- ndim[kk]
    T2U <- graph.T2.test(X1, X2, G=NULL, lfA=lfA, k=kR)
    ##T2U <- T2.test(X1%*%U[, 1:kR], X2%*%U[, 1:kR])
    pU[kk] <- T2U$p.value
  }
  verbose & str(verbose, pU)
  res <- c(res, pU)
  resNames <- c(resNames, paste("T2 (", ndim, " Fourier components)", sep=""))

  ## T2 with PCA, other test statistics
  if (FALSE) {
  verbose && cat(verbose, "T2 on PCA components")
  edX <- svd(rbind(X1, X2))
  E <- edX$v
  pPCA <- rep(NA, length=length(ndim))
  names(pPCA) <- ndim
  for (kk in seq(along=ndim)) {
    kR <- ndim[kk]
    T2PCA <- T2.test(X1%*%E[, 1:kR], X2%*%E[, 1:kR])
    pPCA[kk] <- T2PCA$p.value
  }
  res <- c(res, pPCA)
  resNames <- c(resames, paste("T2 (", ndim, " PCA components)", sep=""))

  verbose & str(verbose, pPCA)

  verbose && cat(verbose, "Adaptive Neyman (Fan)")
  T2AN <- AN.test(X1%*%U, X2%*%U, p)
  pAN <- T2AN$p.value
  verbose & str(verbose, pAN)

  res <- c(res, pAN)
  resNames <- c(resames, "Adaptive Neyman")
  } ## if (FALSE)

  verbose && exit(verbose)

  names(res) <- resNames
  ret <- list(p.value=res, graph=graph, k=ndim)

  ret
}

############################################################################
## HISTORY:
## 2010-10-08
## o Now validating argument 'verbose'.
## 2010-10-01
## o Now manipulates A instead of t(A) (transposition made within
##   laplacianFromA).
## 2010-08-05
## o CLEAN-UP: Now uses 'laplacianFromA'.
## 2010-07-15
## o ROBUSTIFICATION: Now use 'cls <- sort(unique(classes))' to make 'cls'
##   independent of the input order of the label vector.
## 2010-05
## o Created.
############################################################################
