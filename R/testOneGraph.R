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
# @RdocFunction testOneGraph
##
## @title "Applies a serie of two-sample tests to each connected component
##   of a graph using various statistics"
##
## \description{
##  @get "title".
## }
##
## @synopsis
##
## \arguments{
##   \item{graph}{A \code{\link[=graph-class]{graph}} object.}
##   \item{data}{A 'matrix' (size: number 'p' of genes x number 'n' of samples) of gene expression.}
##   \item{classes}{A 'vector' (length: 'n') of class assignments.}
##   \item{useInteractionSigns}{A @logical value indicating whether the sign of interaction should be taken into account.}
##   \item{...}{Further arguments to be passed to testOneConnectedComponent.}
##   \item{verbose}{If @TRUE, extra information is output.}
## }
##
## \value{
##  A structured @list containing the p-values of the tests, the 
##  \code{\link[=graph-class]{graph}} object of the connected component and the number of 
##  retained Fourier dimensions.
## }
##
## @author
##
## \seealso{
##   @see "testOneConnectedComponent"
## }
##
## @examples "../incl/testOneGraph.Rex"
##
##*/########################################################################

testOneGraph <- function(graph, data, classes, useInteractionSigns=TRUE, ..., verbose=FALSE) {
  ## 'graph': an object of class 'graph' or 'graphAM' or 'graphNEL'.
  ## 'data': a 'matrix' (size: number 'p' of genes x number 'n' of samples) of gene expression
  ## 'classes': a 'vector' (length: 'n') of class assignments
  ## 'useInteractionSigns': should the sign of interaction be taken into account ?
  ## '...': further arguments to be passed to testOneConnectedComponent

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

  ## Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose)
  if (verbose) {
    cat <- R.utils::cat
    pushState(verbose)
    on.exit(popState(verbose))
  } 
  
  ## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ## 1. Keep only expression data from nodes in the graph
  ## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Keeping genes in the graph *and* the expression data set")
  dataGN <- rownames(data)
  if(is.NCIgraph(graph))
    {
      graphGN <- translateNCI2GeneID(graph)
      names(graphGN) <- NULL
    }
  else
    graphGN <- translateKEGGID2GeneID(nodes(graph))
  commonGN <- intersect(dataGN, graphGN)

  mm <- match(graphGN, commonGN)
  ww <- which(is.na(mm))
  ll <- length(ww)
  if (ll) {
    verbose && cat(verbose, ll, " genes of the graph were not found in the expression data set: ")
    verbose && str(verbose, graphGN[ww])
  }
  ## sorting (and subsetting if necessary)
  mm <- match(commonGN, graphGN)
  graph <- subGraph(nodes(graph)[mm], graph)
  if (!validGraph(graph)) {
    throw("Argument 'graph' is not a valid graph object")
  }

  mm <- match(dataGN, commonGN)
  ww <- which(is.na(mm))
  ll <- length(ww)
  if (ll) {
    verbose && cat(verbose, ll, " genes of the expression data set are absent from the graph: ")
    verbose && str(verbose, dataGN[ww])
  }
  mm <- match(commonGN, dataGN)
  ## sorting (and subsetting if necessary)
  data <- data[mm, ]

  ## sanity check (ordering should be the same now)
  if(is.NCIgraph(graph))
    {
      graphGN <- translateNCI2GeneID(graph)
      names(graphGN) <- NULL
    }
  else
    graphGN <- translateKEGGID2GeneID(nodes(graph))
  dataGN <- rownames(data)
  stopifnot(all.equal(dataGN, graphGN))

  ## TODO: Rename graph nodes (so that testOneConnectedComponent is not KEGGGraph-specific)
  ##   nodes(graph) <- graphGN
  ## currently cannot be done because we're using getSubtypeDisplay to get the interaction info
  ## and it is not updated when node name change.
  ## OK, we'll do it the other way around:
  rownames(data) <- nodes(graph)

  rm(graphGN)
  rm(dataGN)
  verbose && exit(verbose)

  ## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ## 2. Retrieve interaction signs (if required)
  ## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if (useInteractionSigns) {
    tt <- try(graph <- getSignedGraph(graph, verbose=verbose))
    if (class(tt)=="try-error") {
      warning("Failed to retrieve interaction sign information: returning NULL")
      return(NULL)
    }
  }

  ## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ## After 1. and 2. the graph may have several connected components
  ## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  graphList <- getConnectedComponentList(graph, verbose=verbose)
  verbose && cat(verbose, "Found ", length(graphList), " connected component(s)")
  gsizes <- sapply(graphList, numNodes)
  w1 <- which(gsizes==1)
  if (length(w1)) {
    verbose && cat(verbose, "Excluding ", length(w1), " components with only one node from further analysis")
    graphList <- graphList[-w1]
  }

  ## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ## sign information has been lost in subsetting...
  ## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if (useInteractionSigns) {
    graphList <- lapply(graphList, getSignedGraph, verbose=verbose)
  }

  res <- lapply(graphList, testOneConnectedComponent, data, classes, ..., verbose=verbose)
  res
}

############################################################################
## HISTORY:
## 2011-03-06
## o Dealing with NCI graphs
## 2010-10-08
## o Now validating argument 'verbose'.
############################################################################
