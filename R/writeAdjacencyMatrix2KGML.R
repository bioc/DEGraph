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
## @RdocFunction writeAdjacencyMatrix2KGML
##
## @title "Writes an adjacency matrix into an XML file"
##
## \description{
##  @get "title".
## }
##
## @synopsis
##
## \arguments{
##   \item{mat}{A @matrix, interpreted of the adjacency matrix of a graph.}
##   \item{pathname}{The full path name of the XML file to be written.}
##   \item{nodePrefix}{A @character value giving the prefix to which the node
##     index in 'mat' will be appended.}
##   \item{overwrite}{If @TRUE and file already exists, overwrite it.}
##   \item{...}{Further arguments to be passed to plotKEGGgraph.}
##   \item{verbose}{If @TRUE, extra information is output.}
## }
##
## \value{
##   None.
## }
##
## @author
##
## \seealso{
##   @see "parseKGML2Graph"
## }
##
## @examples "../incl/randomWAMGraph.Rex"
##
##*/########################################################################

writeAdjacencyMatrix2KGML <- function(mat, pathname, nodePrefix="n", overwrite=FALSE, ..., verbose=FALSE){
  ## Validate arguments

  ## Argument 'mat'
  if (!is.element("matrix", class(mat))) {
    throw("Argument 'mat' should be of class 'matrix'")
  }
  nc <- ncol(mat)
  nr <- nrow(mat)
  if (nr != nc) {
    throw("Argument 'mat' should be a *symmetric* matrix, but has ", nr, " rows and ", nc, " columns")
  }
  nnodes <- nr
  rm(nr, nc)
  
  ## Checking that 'mat' only has '-1', '0', and '1' values.
  nmat <- mat
  dim(nmat) <- NULL
  nmat <- Arguments$getIntegers(nmat, range=c(-1, 1))
  rm(nmat)
  
  ## Argument 'pathname'
  pathname <- Arguments$getCharacter(pathname)
  path <- dirname(pathname)
  path <- Arguments$getWritablePath(path)

  ## Argument 'nodePrefix'
  nodePrefix <- Arguments$getCharacter(nodePrefix)

  ## Argument 'overwrite'
  overwrite <- Arguments$getLogical(overwrite)
  if (!overwrite && file.exists(pathname)) {
    throw("File exists and will not be overwritten: ", pathname)
  }

  ## Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose)
  if (verbose) {
    cat <- R.utils::cat
    pushState(verbose)
    on.exit(popState(verbose))
  } 

  ## local 'cat' function:
  locCat <- function(..., append=TRUE) {
    cat(..., "\n", file=pathname, append=append)
  }

  ## nodes
  nodes <- seq(length=nnodes)
  nodeNames <- sprintf("%s%s", nodePrefix, nodes)

  edgeMap <- data.frame(key=c(-1, 1),
                        type=c("inhibition", "activation"),
                        label= c("-|", "-&gt;"))
  
  edgeList <- lapply(nodes, FUN=function(node) {
    mm <- match(mat[node, ], edgeMap[["key"]])
    idxs <- which(!is.na(mm))
    mm <- na.omit(mm)

    res <- nodes[idxs]
    attr(res, "type") <- edgeMap[mm, "type"]
    attr(res, "label") <- edgeMap[mm, "label"]
    res
  })
  names(edgeList) <- nodes
  
  ## write XML file header
  txt <- "<?xml version=\"1.0\"?>"
  locCat(txt, append=FALSE)
  txt <- "<!DOCTYPE pathway SYSTEM \"http://www.genome.jp/kegg/xml/KGML_v0.7.1_.dtd\">"
  locCat(txt)

  ## write pathway
  txt <- "<pathway name=\"path:dummy\" org=\"dummy\" number=\"1\">"
  locCat(txt)

  ## nodes
  entries <- seq(along=nodes)
  for (nn in entries) {
    txt <- sprintf("  <entry id=\"%s\" name=\"%s\" type=\"gene\">", nn, nodeNames[nn])
    locCat(txt)
    txt <- "  </entry>"
    locCat(txt)
  }

  ## edges
  edgesNames <- names(edgeList)
  for (ee1 in seq(along=edgeList)) {
    entry1 <- edgesNames[ee1]
    entries <- edgeList[[ee1]]
    types <- attr(entries, "type")
    labels <- attr(entries, "label")
    for (ee2 in seq(along=entries)) {
      entry2 <- entries[ee2]
      type <- types[ee2]
      label <- labels[ee2]
      txt <- sprintf("  <relation entry1=\"%s\" entry2=\"%s\" type=\"%s\">", entry1, entry2, "NA")
      locCat(txt)
      txt <- sprintf("    <subtype name=\"%s\" value=\"-%s\"/>", type, label)
      locCat(txt)
      txt <- "  </relation>"
      locCat(txt)
    }
  }
  txt <- "</pathway>"
  locCat(txt)
}
  
############################################################################
## HISTORY
## 2010-10-08
## o Now validating argument 'verbose'.
############################################################################
