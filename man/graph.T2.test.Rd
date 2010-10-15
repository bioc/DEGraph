%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  ./graph.T2.test.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{graph.T2.test}
\alias{graph.T2.test}


\title{Performs the Hotelling T2 test in Fourier space}

\description{
 Performs the Hotelling T2 test in Fourier space.
}

\usage{graph.T2.test(X1, X2, G=NULL, lfA=NULL, ..., k=ncol(X1))}

\arguments{
  \item{X1}{A n1 x p \code{\link[base]{numeric}} \code{\link[base]{matrix}}, observed data for class 1: p variables, n1
    observations.}
  \item{X2}{A n2 x p \code{\link[base]{numeric}} \code{\link[base]{matrix}}, observed data for class 2: p variables, n2
    observations.}
  \item{G}{An object of class \code{\link[=graphAM-class]{graphAM}} or
    \code{\link[=graphNEL-class]{graphNEL}}, the graph to be used in the two-sample test.}
  \item{lfA}{A list returned by \code{\link{laplacianFromA}}(),
    containing the Laplacian eigen vectors and eigen values}
  \item{...}{Further arguments to be passed to \code{\link{laplacianFromA}}().}
  \item{k}{A \code{\link[base]{numeric}} value, number of Fourier components retained for the
    test.}
}

\value{
 A \code{\link[base]{list}} with class "htest", as returned by \code{\link{T2.test}}.
}

\author{Laurent Jacob, Pierre Neuvial and Sandrine Dudoit}

\seealso{
  \code{\link{T2.test}}
  \code{\link[=graphAM-class]{graphAM}}
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

