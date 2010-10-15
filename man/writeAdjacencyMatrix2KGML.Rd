%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  ./writeAdjacencyMatrix2KGML.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{writeAdjacencyMatrix2KGML}
\alias{writeAdjacencyMatrix2KGML}


\title{Writes an adjacency matrix into an XML file}

\description{
 Writes an adjacency matrix into an XML file.
}

\usage{writeAdjacencyMatrix2KGML(mat, pathname, nodePrefix="n", overwrite=FALSE, ..., verbose=FALSE)}

\arguments{
  \item{mat}{A \code{\link[base]{matrix}}, interpreted of the adjacency matrix of a graph.}
  \item{pathname}{The full path name of the XML file to be written.}
  \item{nodePrefix}{A \code{\link[base]{character}} value giving the prefix to which the node
    index in 'mat' will be appended.}
  \item{overwrite}{If \code{\link[base:logical]{TRUE}} and file already exists, overwrite it.}
  \item{...}{Further arguments to be passed to plotKEGGgraph.}
  \item{verbose}{If \code{\link[base:logical]{TRUE}}, extra information is output.}
}

\value{
  None.
}

\author{Laurent Jacob, Pierre Neuvial and Sandrine Dudoit}

\seealso{
  \code{\link{parseKGML2Graph}}
}

\examples{
library("KEGGgraph")
library("rrcov")

## Create a random graph
graph <- randomWAMGraph(nnodes=5, nedges=7, verbose=TRUE)
plot(graph)

## Retrieve its adjacency matrix
A <- graph@adjMat

## write it to KGML file
grPathname <- "randomWAMGraph.xml"
writeAdjacencyMatrix2KGML(A, pathname=grPathname, verbose=TRUE, overwrite=TRUE)

## read it from file
gr <- parseKGML2Graph(grPathname)

## Two examples of Laplacians from the same graph
lapMI <- laplacianFromA(A, ltype="meanInfluence")
print(lapMI)

lapN <- laplacianFromA(A, ltype="normalized")
print(lapN)

U <- lapN$U
p <- nrow(A)
sigma <- diag(p)/sqrt(p)

X <- twoSampleFromGraph(100, 120, shiftM2=1, sigma, U=U, k=3)

## T2
t <- T2.test(X$X1,X$X2)
str(t)

tu <- graph.T2.test(X$X1, X$X2, lfA=lapMI, k=3)
str(tu)
}

