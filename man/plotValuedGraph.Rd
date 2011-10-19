%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  ./plotValuedGraph.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{plotValuedGraph}
\alias{plotValuedGraph}


\title{Plots a graph with nodes colored according to a quantitative variable}

\description{
 Plots a graph with nodes colored according to a quantitative variable.
}

\usage{plotValuedGraph(graph, values=NULL, nodeLabels=nodes(graph), qMax=0.95, colorPalette=heat.colors(10), adjustColorRange=FALSE, symmetrizeArrows=FALSE, height=1, lwd=1, cex=1, ..., verbose=FALSE)}

\arguments{
  \item{graph}{A \code{\link[=graph-class]{graph}} object.}
  \item{values}{A named \code{\link[base]{vector}} of \code{\link[base]{numeric}} values according to which the
    graph nodes should be colored.}
  \item{nodeLabels}{A \code{\link[base]{character}} \code{\link[base]{vector}} of the same length and in the
    same order as 'nodes(graph)': node labels to be displayed.  Defaults
    to 'nodes(graph)'.}
  \item{qMax}{A \code{\link[base]{numeric}} value, fraction of the data to be truncated in order
    to avoid outliers.}
  \item{colorPalette}{A \code{\link[base]{character}} vector, the set of colors to be used.}
  \item{adjustColorRange}{A \code{\link[base]{logical}} value.  If \code{\link[base:logical]{TRUE}}, the color range is
    adjusted to the range of values of nodes actually present in the graph.
    Defaults to \code{\link[base:logical]{FALSE}}, i.e. the color range spans range(values) regardless
    of which nodes are present in the graph.}
  \item{symmetrizeArrows}{A \code{\link[base]{logical}} value.  If \code{\link[base:logical]{TRUE}}, arrow tails are
    drawn as the corresponding arrow heads.  Defaults to \code{\link[base:logical]{FALSE}}.}
  \item{height}{A \code{\link[base]{numeric}} value, the (common) size of nodes.}
  \item{lwd}{A \code{\link[base]{numeric}} value, the (common) width of edges.}
  \item{cex}{A \code{\link[base]{numeric}} value, the relative size of the text for gene names.}
  \item{...}{Further arguments to be passed to 'edgeRenderInfo' and
    'nodeRenderInfo'.}
  \item{verbose}{If \code{\link[base:logical]{TRUE}}, extra information is output.}
}

\value{
  A \code{\link[base]{list}} containing the following components:
  \describe{
    \item{graph}{The 'graph' object as plotted.}
    \item{breaks}{The break points in the supplied values (can be used for
      plotting a legend).}
  }
}

\author{Laurent Jacob, Pierre Neuvial and Sandrine Dudoit}

\seealso{
  \code{\link{plotKEGGgraph}}
  \code{\link{plot}}()
}

\examples{
library("Rgraphviz")
library("KEGGgraph")
library("NCIgraph")

data("Loi2008_DEGraphVignette")
exprData <- exprLoi2008
classData <- classLoi2008
annData <- annLoi2008

rn <- rownames(exprData)

## Retrieve expression levels data for genes from one KEGG pathway
graph <- grListKEGG[[1]]
pname <- attr(graph, "label")
print(pname)

## DEGraph T2 test
resList <- testOneGraph(graph, exprData, classData, verbose=TRUE, prop=0.2)

## Largest connected component
res <- resList[[1]]
gr <- res$graph

## individual t statistics
shift <- apply(exprData, 1, FUN=function(x) {
  tt <- t.test(x[classData==0], x[classData==1])
  tt$statistic
})
names(shift) <- translateGeneID2KEGGID(names(shift))

## color palette
if (require(marray)) {
  pal <- maPalette(low="red", high="green", mid="black", k=100)
} else {
  pal <- heat.colors(100)
}

## plot results
dn <- getDisplayName(gr, shortLabel=TRUE)
mm <- match(translateKEGGID2GeneID(nodes(gr)), rownames(annData))
dn <- annData[mm, "NCBI.gene.symbol"]
 
pvg <- plotValuedGraph(gr, values=shift, nodeLabels=dn, qMax=0.95, colorPalette=pal, height=40, lwd=1, verbose=TRUE, cex=0.5)
title(pname)

txt1 <- sprintf("p(T2)=\%s", signif(res$p.value[1], 2))
txt2 <- sprintf("p(T2F[\%s])=\%s", res$k, signif(res$p.value[2]))
txt <- paste(txt1, txt2, sep="\n")
stext(side=3, pos=1, txt)
if (require(fields)) {
  image.plot(legend.only=TRUE, zlim=range(pvg$breaks), col=pal, legend.shrink=0.3, legend.width=0.8, legend.lab="t-scores", legend.mar=3.3)
}
}

