%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  ./hyper.test.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{hyper.test}
\alias{hyper.test}


\title{Performs an hypergeometric test of enrichment of a set of
  hypotheses in significant elements}

\description{
 Performs an hypergeometric test of enrichment of a set of
  hypotheses in significant elements.
}

\usage{hyper.test(p.values, testSet, thr=0.001, universe=length(p.values), verbose=FALSE)}

\arguments{
  \item{p.values}{A named \code{\link[base]{numeric}} vector giving the p-values of all
    tested elements.}
  \item{testSet}{A \code{\link[base]{character}} vector giving the ids of the elements in the
    tested set. Elements of 'testSet' must have a match in 'names(p.values)'.}
  \item{thr}{A \code{\link[base]{numeric}} value between 0 and 1 giving the threshold on
    p-values at which an element is declared to be significant.}
  \item{universe}{An \code{\link[base]{integer}} value giving the number of elelments in the
    considered universe.  Defaults to 'length(p.values)'.}
  \item{verbose}{If \code{\link[base:logical]{TRUE}}, extra information is output.}
}

\value{
 A \code{\link[base]{list}} with class "htest" containing the following components:
 \describe{
   \item{statistic}{A \code{\link[base]{numeric}} value, the test statistic.}
   \item{p.value}{A \code{\link[base]{numeric}} value, the corresponding p-value.}
  }
}

\author{Laurent Jacob, Pierre Neuvial and Sandrine Dudoit}

\seealso{
  \code{\link{AN.test}}()
  \code{\link{BS.test}}()
  \code{\link{graph.T2.test}}()
}

\examples{
library("KEGGgraph")
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

