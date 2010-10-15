%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  ./getConnectedComponentList.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{getConnectedComponentList}
\alias{getConnectedComponentList}


\title{Given a graph, returns a list of its connected components
(which are also graph objects), ordered by decreasing number of
nodes}

\description{
 Given a graph, returns a list of its connected components
(which are also graph objects), ordered by decreasing number of
nodes.
}

\usage{getConnectedComponentList(graph, verbose=FALSE)}

\arguments{
  \item{graph}{A \code{\link[=graph-class]{graph}} object.}
  \item{verbose}{If \code{\link[base:logical]{TRUE}}, extra information is output.}
}

\value{ A \code{\link[base]{list}} containing a \code{\link[=graph-class]{graph}} object for each connected
 component of the input graph, ordered by decreasing number of nodes
 }

\author{Laurent Jacob, Pierre Neuvial and Sandrine Dudoit}

\seealso{\code{\link{connectedComp}}.}

\examples{
data("Loi2008_DEGraphVignette")
exprData <- exprLoi2008
rn <- rownames(exprData)

## Retrieve expression levels data for genes from one KEGG pathway
graph <- grListKEGG[[1]]
pname <- attr(graph, "label")
cat(verbose, "Pathway name: ", pname)

sgraph <- getSignedGraph(graph, verbose=TRUE)
print(sgraph)

graphList <- getConnectedComponentList(graph, verbose=TRUE)
print(graphList)

}

