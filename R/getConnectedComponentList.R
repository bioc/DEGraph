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
## @RdocFunction getConnectedComponentList
##
## @title "Given a graph, returns a list of its connected components
## (which are also graph objects), ordered by decreasing number of
## nodes"
##
## \description{
##  @get "title".
## }
##
## @synopsis
##
## \arguments{
##   \item{graph}{A \code{\link[=graph-class]{graph}} object.}
##   \item{verbose}{If @TRUE, extra information is output.}
## }
##
## \value{ A @list containing a \code{\link[=graph-class]{graph}} object for each connected
##  component of the input graph, ordered by decreasing number of nodes
##  }
##
## @author
##
## \seealso{@see "connectedComp".}
##
## @examples "../incl/getSignedGraph.Rex"
##
##*/#########################################################################

getConnectedComponentList <- function(graph, verbose=FALSE) {
  ## Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose)
  if (verbose) {
    pushState(verbose)
    on.exit(popState(verbose))
  } 

  nodeList <- connectedComp(graph)
  sizes <- sapply(nodeList, length)
  oo <- order(sizes, decreasing=TRUE)
  graphList <- lapply(nodeList[oo], subGraph, graph)
  verbose && print(verbose, graphList)
  graphList
}

############################################################################
## HISTORY:
## 2010-10-08
## o Now validating argument 'verbose'.
## 2010-05-01
## o Created.
############################################################################
