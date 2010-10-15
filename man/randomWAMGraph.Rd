%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  ./randomWAMGraph.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{randomWAMGraph}
\alias{randomWAMGraph}


\title{Generates a random graph}

\description{
 Generates a random graph.
}

\usage{randomWAMGraph(nnodes=5, nedges=nnodes, verbose=FALSE)}

\arguments{
  \item{nnodes}{A \code{\link[base]{numeric}} value, the desired number of nodes.}
  \item{nedges}{A \code{\link[base]{numeric}} value, the desired number of edges.}
  \item{verbose}{If \code{\link[base:logical]{TRUE}}, extra information is output.}
}

\value{
 An object of class \code{\link[=graphAM-class]{graphAM}}.
}

\author{Laurent Jacob, Pierre Neuvial and Sandrine Dudoit}

\seealso{
  \code{\link[=graphAM-class]{graphAM}}.
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

