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
## @RdocFunction getSignedGraph
##
## @title "Given a graph, builds a signed version of the adjacency matrix
##   taking into account the type of interaction (e.g., activation or
##   inhibition)"
##
## \description{
##  @get "title".
## }
##
## @synopsis
##
## \arguments{
##   \item{graph}{A \code{\link[=graph-class]{graph}} object.}
##   \item{positiveInteractionLabels}{A @character @vector specifying which
##     interaction labels correspond to positive interactions. Defaults to
##     'c("activation", "expression")'.}
##   \item{negativeInteractionLabels}{A @character @vector specifying which
##     interaction labels correspond to negative interactions. Defaults to
##     'c("inhibition", "repression")'.}
##   \item{verbose}{If @TRUE, extra information is output.}
## }
##
## \value{
##   This function returns a squared matrix whose (i,j) entry is:
##   \describe{
##     \item{0}{if edges i and j are not connected}
##     \item{1}{if edges i and j are connected by a positive interaction}
##     \item{-1}{if edges i and j are connected by a negative interaction.}
##   }
##   By construction, the absolute value of this matrix is the adjacency
##   matrix of the graph. Edges which cannot interpreted as corresponding
##   to a positive or a negative interaction are marked as not connected.
## }
##
## @author
##
## @examples "../incl/getSignedGraph.Rex"
##
##*/########################################################################

## NOTES:
## adjMat cannot be signed in KEGGgraph objects.
## Long-term fix: implement a 'signedGraph' class that allows that
## Short-term fix: add an attribute 'signMat' to the graph object
##   and use it manually afterwards (it is not properly subsetted
##   when removing nodes from the graph for example.  Actually it
##   is even lost)

getSignedGraph <- function(graph, positiveInteractionLabels=c("activation", "expression"), negativeInteractionLabels=c("inhibition", "repression"), verbose=FALSE) {
  ## - - - - - - - - - - - - - - - - - - -
  ## Validate arguments
  ## - - - - - - - - - - - - - - - - - - -
  ## Argument 'graph'
  if (!validGraph(graph)) {
    throw("Argument 'graph' is not a valid graph object")
  }
  if (!isDirected(graph)) {
    throw("Undirected graphs are not supported yet")
  }
  ## Argument 'positiveInteractionLabels'
  posIntLabs <- Arguments$getCharacters(positiveInteractionLabels)

  ## Argument 'negativeInteractionLabels'
  negIntLabs <- Arguments$getCharacters(negativeInteractionLabels)
  
  ## Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose)
  if (verbose) {
    cat <- R.utils::cat
    pushState(verbose)
    on.exit(popState(verbose))
  } 

  
  ## retrieve edge types
  st <- getSubtype(graph)  ## FIXME: does not work for undirected graphs
  edgeNames <- names(st)
  intNames <- sapply(st, FUN=function(x) {
    x$subtype@name
  })

  verbose && cat(verbose, "Interactions:")
  verbose && str(verbose, intNames)
  verbose && print(verbose, table(intNames))

  lap <- lapply(posIntLabs, grep, intNames)
  names(lap) <- posIntLabs
  verbose && cat(verbose, "Positive interactions:")
  verbose && str(verbose, lap)
  idxPos <- unlist(lap) 
  
  lan <- lapply(negIntLabs, grep, intNames)
  names(lan) <- negIntLabs
  verbose && cat(verbose, "Negative interactions:")
  verbose && str(verbose, lan)
  idxNeg <- unlist(lan)

  posNames <- edgeNames[idxPos]
  negNames <- edgeNames[idxNeg]

  idx0 <- intersect(idxPos, idxNeg)
  if (length(idx0)) {
      warning("Probable annotation error: positive and negative edges. These edges will be removed")
      idxPos <- idxPos[-match(idx0, idxPos)]
      idxNeg <- idxNeg[-match(idx0, idxNeg)]
  }

  ## Remove other edges
  edgesToRemove <- intNames[-c(idxPos, idxNeg)]
  if (length(edgesToRemove)) {
    verbose && enter(verbose, "Removing unsigned edges from graph")
    verbose && cat(verbose, "Edges to remove:")
    verbose && str(verbose, edgesToRemove)
    verbose && print(verbose, table(edgesToRemove))

    etrNames <- edgeNames[-c(idxPos, idxNeg)]
    ftNames <- sapply(etrNames, FUN=function(edge) {
      unlist(strsplit(edge, "~"))
    })

    oldGraph <- graph
    graph <- removeEdge(from=ftNames[1, ], ftNames[2, ], graph)
    verbose && print(verbose, graph)

    rm(ftNames)
    verbose && exit(verbose)
  }

  verbose && enter(verbose, "Creating sign matrix")
  ## use adjacency matrix of the corresponding undirected graph
  ## ugraph <- ugraph(graph)
  ## graphAM <- as(ugraph, "graphAM")
  ## sanity check: symmetry
  ## stopifnot(sum(sum(t(adjMat)!=adjMat))==0)
  
  graphAM <- as(graph, "graphAM") # Now use directed version as some energy functions require edge directions
  
  adjMat <- graphAM@adjMat

  signMat <- adjMat

  if (length(negNames)) {
    ## set negative edges to -1
    ftNames <- sapply(negNames, FUN=function(edge) {
      unlist(strsplit(edge, "~"))
    })
    nodeNames <- colnames(signMat)
    for (ii in seq(length=ncol(ftNames))) {
      ftName <- ftNames[, ii]
      idx <- graph:::getIndices(nodeNames, ftName[1], ftName[2])
      adj <- adjMat[idx$from, idx$to]
      if (adj != 1) { ## Sanity check
        throw("Inconsistent data in adjacency matrix for edge ", paste(ftName, collapse="~"))
      } else {
        signMat[idx$from, idx$to] <- -1
        ##signMat[idx$to, idx$from] <- -1
      }
    }
  }

  if (length(posNames)) {
    ## check that positive edges are 1
    ftNames <- sapply(posNames, FUN=function(edge) {
      unlist(strsplit(edge, "~"))
    })
    nodeNames <- colnames(signMat)
    for (ii in seq(length=ncol(ftNames))) {
      ftName <- ftNames[, ii]
      idx <- graph:::getIndices(nodeNames, ftName[1], ftName[2])
      adj <- adjMat[idx$from, idx$to]
      if (adj != 1) { ## Sanity check
        throw("Inconsistent data in adjacency matrix for edge", ftName)
      }
    }
  }
  
  attr(graph, 'signMat') <- signMat
  verbose && exit(verbose)

  graph
}

############################################################################
## HISTORY:
## 2010-10-08
## o Now validating argument 'verbose'.
## o Updated arguments.
## 2010-09-17
## o Now returns an directed graph (non symmetric adjacency matrix that
##   takes edge direction into account).
## 2010-07-16
## o Now returns an undirected graph (so that adjacency matrix is symmetric).
##   
## 2010-05
## o Created.
############################################################################
