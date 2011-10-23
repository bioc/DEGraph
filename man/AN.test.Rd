%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  ./AN.test.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{AN.test}
\alias{AN.test}


\title{Performs the Adaptive Neyman test of Fan and Lin (1998)}

\description{
 Performs the Adaptive Neyman test of Fan and Lin (1998).
}

\usage{AN.test(X1, X2, candK=1:ncol(X1), na.rm=FALSE)}

\arguments{
  \item{X1}{A n1 x p \code{\link[base]{matrix}}, observed data for class 1: p variables, n1
    observations.}
  \item{X2}{A n2 x p \code{\link[base]{matrix}}, observed data for class 2: p variables, n2
    observations.}
  \item{candK}{A \code{\link[base]{vector}}, candidate values for the true number of Fourier
    components.}
  \item{na.rm}{A \code{\link[base]{logical}} value indicating whether variables with \code{\link[base]{NA}} in
    at least one of the n1 + n2 observations should be discarder before
    the test is performed.}
}

\value{
 A \code{\link[base]{list}} with class "htest" containing the following components:
 \describe{
   \item{statistic}{A \code{\link[base]{numeric}} value, the test statistic.}
   \item{p.value}{A \code{\link[base]{numeric}} value, the corresponding p-value.}
   \item{kstar}{A \code{\link[base]{numeric}} value, the estimated true number of Fourier
     components.}
  }
}

\author{Laurent Jacob, Pierre Neuvial and Sandrine Dudoit}

\seealso{
  \code{\link{BS.test}}()
  \code{\link{graph.T2.test}}()
  \code{\link{hyper.test}}()
}

\examples{
library("KEGGgraph")
## library("NCIgraph")
library("rrcov")

data("Loi2008_DEGraphVignette")
exprData <- exprLoi2008
classData <- classLoi2008
rn <- rownames(exprData)

## Retrieve expression levels data for genes from one KEGG pathway
gr <- grListKEGG[[1]]
gids <- translateKEGGID2GeneID(nodes(gr))
mm <- match(gids, rownames(exprData))

## Keep genes from the graph that are present in the expression data set
idxs <- which(!is.na(mm))
gr <- subGraph(nodes(gr)[idxs], gr)

idxs <- which(is.na(mm))
if(length(idxs)) {
  print("Gene ID not found in expression data: ")
  str(gids[idxs])
}
dat <- exprData[na.omit(mm), ]
str(dat)

X1 <- t(dat[, classData==0])
X2 <- t(dat[, classData==1])

## DEGraph T2 test
res <- testOneGraph(gr, exprData, classData, verbose=TRUE, prop=0.2)

## T2 test (Hotelling)
rT2 <- T2.test(X1, X2)
str(rT2)

## Adaptive Neyman test
rAN <- AN.test(X1, X2, na.rm=TRUE)
str(rAN)

## Adaptive Neyman test from Fan and Lin (1998)
rAN <- AN.test(X1, X2, na.rm=TRUE)
str(rAN)

## Test from Bai and Saranadasa (1996)
rBS <- BS.test(X1, X2, na.rm=TRUE)
str(rBS)

## Hypergeometric test
pValues <- apply(exprData, 1, FUN=function(x) {
  tt <- t.test(x[classData==0], x[classData==1])
  tt$p.value
})
str(pValues)
names(pValues) <- rownames(exprData)
rHyper <- hyper.test(pValues, gids, thr=0.01)
str(rHyper)

}

