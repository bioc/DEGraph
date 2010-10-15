%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  ./getSignedGraph.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{getSignedGraph}
\alias{getSignedGraph}


\title{Given a graph, builds a signed version of the adjacency matrix
  taking into account the type of interaction (e.g., activation or
  inhibition)}

\description{
 Given a graph, builds a signed version of the adjacency matrix
  taking into account the type of interaction (e.g., activation or
  inhibition).
}

\usage{getSignedGraph(graph, positiveInteractionLabels=c("activation", "expression"), negativeInteractionLabels=c("inhibition", "repression"), verbose=FALSE)}

\arguments{
  \item{graph}{A \code{\link[=graph-class]{graph}} object.}
  \item{positiveInteractionLabels}{A \code{\link[base]{character}} \code{\link[base]{vector}} specifying which
    interaction labels correspond to positive interactions. Defaults to
    'c("activation", "expression")'.}
  \item{negativeInteractionLabels}{A \code{\link[base]{character}} \code{\link[base]{vector}} specifying which
    interaction labels correspond to negative interactions. Defaults to
    'c("inhibition", "repression")'.}
  \item{verbose}{If \code{\link[base:logical]{TRUE}}, extra information is output.}
}

\value{
  This function returns a squared matrix whose (i,j) entry is:
  \describe{
    \item{0}{if edges i and j are not connected}
    \item{1}{if edges i and j are connected by a positive interaction}
    \item{-1}{if edges i and j are connected by a negative interaction.}
  }
  By construction, the absolute value of this matrix is the adjacency
  matrix of the graph. Edges which cannot interpreted as corresponding
  to a positive or a negative interaction are marked as not connected.
}

\author{Laurent Jacob, Pierre Neuvial and Sandrine Dudoit}

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

