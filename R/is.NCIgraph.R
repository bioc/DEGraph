## Copyright 2011 Laurent Jacob

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
# @RdocFunction is.NCIgraph
##
## @title "Assess whether a graph is a NCI graph"
##
##
## \description{
##  @get "title".
## }
##
## @synopsis
##
## \arguments{
##   \item{gr}{A \code{\link[=graph-class]{graph}} object.}
## }
##
## \value{A @logical, @TRUE if the graph is a NCI graph, @FALSE otherwise.}
##
## @author
##
## \seealso{
##   @see "parseNCInetwork"
## }
##
##*/########################################################################

is.NCIgraph <- function(gr)
  {
    ##(length(gr@nodeData@data) > 0) && (length(grep('biopax',names(gr@nodeData@data[[1]]))) > 0)
    isTRUE(class(gr) == 'NCIgraph')
  }
