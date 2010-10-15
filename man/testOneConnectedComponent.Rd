%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  ./testOneConnectedComponent.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{testOneConnectedComponent}
\alias{testOneConnectedComponent}


\title{Applies a series of two-sample tests to a connected graph using various statistics}

\description{
 Applies a series of two-sample tests to a connected graph using various statistics.
}

\usage{testOneConnectedComponent(graph, data, classes, ..., prop=0.2, verbose=FALSE)}

\arguments{
  \item{graph}{A \code{\link[=graph-class]{graph}} object.}
  \item{data}{A '\code{\link[base]{numeric}} \code{\link[base]{matrix}} (size: number 'p' of genes x number 'n' of samples) of gene expression.}
  \item{classes}{A \code{\link[base]{character}} \code{\link[base]{vector}} (length: 'n') of class assignments.}
  \item{...}{Further arguments to be passed to \code{\link{laplacianFromA}}().}
  \item{prop}{A \code{\link[base]{numeric}} value, percentage of components retained for Fourier and PCA.}
  \item{verbose}{If \code{\link[base:logical]{TRUE}}, extra information is output.}
}

\value{
 A structured \code{\link[base]{list}} containing the p-values of the tests, the
 \code{\link[=graph-class]{graph}} object of the connected component and the number
 of retained Fourier dimensions.
}

\author{Laurent Jacob, Pierre Neuvial and Sandrine Dudoit}

\details{
  This function performs the test, assuming that all genes
  in the graph are represented in the expression data set,
  in order not to have to modify the graph topology.

  Interaction signs are used if available in the graph
  ('getSignedGraph' is not called here, in order not to
  have to modify the graph topology.).

  The graph given as input has to have only one
  connex component.  It can be retrieved from the output of
  \code{\link{getConnectedComponentList}}().
}

\seealso{
  \code{\link{testOneGraph}}()
  \code{\link{getConnectedComponentList}}()
}

\examples{
library("rrcov")

## Some parameters
n1 <- n2 <- 20
nnodes <- nedges <- 20
k <- 3
ncp <- 0.5
sigma <- diag(nnodes)/sqrt(nnodes)


## Build graph, decompose laplacian
G <- randomWAMGraph(nnodes=nnodes,nedges=nedges)
A <- G@adjMat
lfA <- laplacianFromA(A,ltype="unnormalized")
U <- lfA$U
l <- lfA$l

## Build two samples with smooth mean shift
X <- twoSampleFromGraph(n1,n2,shiftM2=ncp,sigma,U=U,k=k)

## Do hypothesis testing
t <- T2.test(X$X1,X$X2) # Raw T-square
print(t$p.value)
tu <- graph.T2.test(X$X1,X$X2,lfA=lfA,k=k) # Filtered T-squares
print(tu$p.value)
}

