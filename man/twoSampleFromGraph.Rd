%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  ./twoSampleFromGraph.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{twoSampleFromGraph}
\alias{twoSampleFromGraph}


\title{Given a basis (typically the eigenvectors of a graph Laplacian), builds two multivariate normal samples with mean shift located in the first elements of the basis}

\description{
 Given a basis (typically the eigenvectors of a graph Laplacian), builds two multivariate normal samples with mean shift located in the first elements of the basis.
}

\usage{twoSampleFromGraph(n1=20, n2=n1, shiftM2=0, sigma, U, k=ceiling(ncol(U)/3))}

\arguments{
  \item{n1}{An \code{\link[base]{integer}} value specifying the number of points in the first sample.}
  \item{n2}{An \code{\link[base]{integer}} value specifying the number of points in the second sample.}
  \item{shiftM2}{A \code{\link[base]{numeric}} value giving the desired squared Mahalanobis norm of the mean shift between the two samples.}
  \item{sigma}{A matrix giving the covariance structure of each sample.}
  \item{U}{A matrix giving the desired basis.}
  \item{k}{An \code{\link[base]{integer}} value giving the number of basis elements in which the mean shift must be located.}
}

\value{ A \code{\link[base]{list}} with named elements:
 \describe{
   \item{X1}{The first sample in the original basis (before transformation by U).}
   \item{X2}{The second sample in the original basis (before transformation by U).}
   \item{X1}{The first sample in the specified basis (after transformation by U).}
   \item{X2}{The second sample in the specified basis (after transformation by U).}
   \item{mu1}{The population mean of F1}
   \item{mu2}{The population mean of F2}
   \item{diff}{mu1 - mu2}
  }
 }

\author{Laurent Jacob, Pierre Neuvial and Sandrine Dudoit}

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

