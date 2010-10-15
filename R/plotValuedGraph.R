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
## @RdocFunction plotValuedGraph
##
## @title "Plots a graph with nodes colored according to a quantitative variable"
##
## \description{
##  @get "title".
## }
##
## @synopsis
##
## \arguments{
##   \item{graph}{A \code{\link[=graph-class]{graph}} object.}
##   \item{values}{A named @vector of @numeric values according to which the
##     graph nodes should be colored.}
##   \item{nodeLabels}{A @character @vector of the same length and in the
##     same order as 'nodes(graph)': node labels to be displayed.  Defaults
##     to 'nodes(graph)'.}
##   \item{qMax}{A @numeric value, fraction of the data to be truncated in order
##     to avoid outliers.}
##   \item{colorPalette}{A @character vector, the set of colors to be used.}
##   \item{adjustColorRange}{A @logical value.  If @TRUE, the color range is
##     adjusted to the range of values of nodes actually present in the graph.
##     Defaults to @FALSE, i.e. the color range spans range(values) regardless
##     of which nodes are present in the graph.}
##   \item{symmetrizeArrows}{A @logical value.  If @TRUE, arrow tails are
##     drawn as the corresponding arrow heads.  Defaults to @FALSE.}
##   \item{height}{A @numeric value, the (common) size of nodes.}
##   \item{lwd}{A @numeric value, the (common) width of edges.}
##   \item{cex}{A @numeric value, the relative size of the text for gene names.}
##   \item{...}{Further arguments to be passed to 'edgeRenderInfo' and
##     'nodeRenderInfo'.}
##   \item{verbose}{If @TRUE, extra information is output.}
## }
##
## \value{
##   A @list containing the following components:
##   \describe{
##     \item{graph}{The 'graph' object as plotted.}
##     \item{breaks}{The break points in the supplied values (can be used for
##       plotting a legend).}
##   }
## }
##
## @author
##
## \seealso{
##   @see "plotKEGGgraph"
##   @see "plot"
## }
##
## @examples "../incl/testOneGraph.Rex"
##
##*/########################################################################

plotValuedGraph <- function(graph, values=NULL, nodeLabels=nodes(graph), qMax=0.95, colorPalette=heat.colors(10), adjustColorRange=FALSE, symmetrizeArrows=FALSE, height=1, lwd=1, cex=1, ..., verbose=FALSE){
  
  ##par(oma=c( 0,0,0,4))

  ## Validate arguments

  ## Argument 'graph'
  if (!inherits(graph, "graph")) {
    throw("Argument 'graph' should derive from class 'graph'")
  }
  gnodes <- nodes(graph)
  nnodes <- length(gnodes)

  ## Argument 'values'
  values <- Arguments$getNumerics(values)
  vnodes <- names(values)
  if (length(values) && is.null(vnodes)) {
    throw("Names of argument 'values' should be non NULL")
  }

  ## Argument 'nodeLabels'
  nodeLabels <- Arguments$getCharacters(nodeLabels)
  if (length(nodeLabels) != nnodes) {
    throw("Length of argument 'nodeLabels' should match the number of nodes in the graph")
  }
  nLabels <- nodeLabels
  names(nLabels) <- gnodes

  ## Argument 'qMax'
  qMax <- Arguments$getNumeric(qMax)

  ## Argument 'colorPalette'
  colorPalette <- Arguments$getCharacters(colorPalette)

  ## Argument 'adjustColorRange'
  adjustColorRange <- Arguments$getLogical(adjustColorRange)

  ## Argument 'symmetrizeArrows'
  symmetrizeArrows <- Arguments$getLogical(symmetrizeArrows)

  ## Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose)
  if (verbose) {
    cat <- R.utils::cat
    pushState(verbose)
    on.exit(popState(verbose))
  } 

  verbose && cat(verbose, "Nodes and their labels")
  verbose && str(verbose, nLabels)

  ## Associate values to the corresponding nodes on the graph
  graphValues <- rep(NA, nnodes)
  names(graphValues) <- gnodes

  cnodes <- intersect(vnodes, gnodes)
  graphValues[match(cnodes, gnodes)] <- values[cnodes]

  if (length(values)) {
    if (adjustColorRange) {
      cs <- graphValues
    } else {  ## color range from all the values
      cs <- values
    }
    cs <- abs(cs)  ## enforce color scale symmetry
    MM <- quantile(cs, qMax, na.rm=TRUE)  ## try to avoid outliers
    
    ## truncate outliers
    graphValues[graphValues< -MM] <- -MM
    graphValues[graphValues> MM] <- MM
    
    breaks <- seq(from=-MM, to=MM, length=length(colorPalette))
    verbose && cat(verbose, "Color scale breaks")
    verbose && str(verbose, breaks)
    
    nodeCols <- level.colors(graphValues, at=breaks, col.regions=colorPalette)
    names(nodeCols) <- names(graphValues)
    verbose && cat(verbose, "Node colors")
    verbose && str(verbose, nodeCols)
  } else {
    breaks=NULL
  }

  ##par(mfrow=c(1,2))
  ##par(mar = c(5, 0, 0, 5))
  ##image(cbind(1L:length(pal)), col = pal, axes = FALSE)
  ##par(mar = c(0, 0, 0, 0))
  
  ed <- edgeData(graph)
  ke <- ed[[1]]$KEGGEdge
  if (!is.null(ke)) {  ## only way to know if graph is KEGGgraph-compliant
    ## BEGIN code borrowed from plotKEGGgraph
    subdisplay <- subtypeDisplay(graph)
    eLabel <- subdisplay["label", ]
    eCol <- subdisplay["color", ]
    eTextCol <- subdisplay["fontcolor", ]
    eLty <- subdisplay["style", ]
    eArrowhead <- subdisplay["arrowhead", ]
    eArrowhead[eArrowhead=="normal"] <- "normArrow"  
    if (ncol(subdisplay) == 1) {
      tmp <- colnames(subdisplay)[1]
      names(eLabel) <- names(eCol) <- names(eTextCol) <- tmp
      names(eLty) <- names(eArrowhead) <- tmp
    }
    edgeRenderInfo(graph) <- list(lty=eLty, col=eCol, textCol=eTextCol, 
                                  label=eLabel, arrowhead=eArrowhead, label=eLabel)
    if (symmetrizeArrows) {
      edgeRenderInfo(graph) <- list(arrowtail=eArrowhead)
    }
  } else {
    ## graphAM
    if (inherits(graph, "graphAM")) {
     ahd <- rep("normal", length(ed))
     names(ahd) <- names(ed)
     adjMat <- graph@adjMat
     ## not useful: adjacency matrix is not signed...
   }
  }
  graph <- layoutGraph(graph)
  nodeRenderInfo(graph) <- list(label=nLabels, height=height, cex=cex, ...)
  if (length(values)) {
    nodeRenderInfo(graph) <- list(fill=nodeCols,textCol="black")
  }
  edgeRenderInfo(graph) <- list(lwd=lwd, ...)
  renderGraph(graph)  
  ## END code borrowed from plotKEGGgraph

  ##par(oma=c( 0,0,15,1))# reset margin to be much smaller.
  ##image.plot(legend.only=TRUE, zlim=range(breaks), col=colorPalette, legend.shrink=0.3, legend.width=0.8, legend.lab="t-scores", legend.mar=5) 
  ##set.panel() # reset plotting device

  invisible(list(graph=graph, breaks=breaks))
}

############################################################################
## HISTORY
## 2010-10-08
## o Now validating argument 'verbose'.
## 2010-09-17
## o Fixed node labels.
## 2010-09-16
## o Removed dependency on KEGGgraph.
## o Added option 'symmetrizeArrows'.
## 2010-09-07
## o Added option 'translateGeneIDs'.
## o Added Rdoc.
## o removed brewer.pal to avoid depending on RColorBrewer.
## o renamed 'shift' into 'values' for more general applicability.
## 2010-08-05
## o Color scale is now symmetric / 0.
## o Legend is not drawn within plotDERes anymore.
## 2010-08-04
## o BUG FIX: colors were wrong (due to a factor/character problem).
############################################################################
