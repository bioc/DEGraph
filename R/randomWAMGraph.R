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
## @RdocFunction randomWAMGraph
##
## @title "Generates a random graph"
##
## \description{
##  @get "title".
## }
##
## @synopsis
##
## \arguments{
##   \item{nnodes}{A @numeric value, the desired number of nodes.}
##   \item{nedges}{A @numeric value, the desired number of edges.}
##   \item{verbose}{If @TRUE, extra information is output.}
## }
##
## \value{
##  An object of class \code{\link[=graphAM-class]{graphAM}}.
## }
##
## @author
##
## \seealso{
##   \code{\link[=graphAM-class]{graphAM}}.
## }
##
## @examples "../incl/randomWAMGraph.Rex"
##
##*/########################################################################

randomWAMGraph <- function(nnodes=5, nedges=nnodes, verbose=FALSE){
  ## Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose)
  if (verbose) {
    cat <- R.utils::cat
    pushState(verbose)
    on.exit(popState(verbose))
  } 

  i <- 0

  repeat{
    g <- randomEGraph(as.character(1:nnodes),edges=nedges)
    if(isConnected(g))
      break
    i <- i+1
    verbose && print(verbose, i)
  }

  g <- as(g,"graphAM")
}


############################################################################
## HISTORY
## 2010-10-08
## o Now validating argument 'verbose'.
############################################################################
