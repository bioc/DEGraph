library("KEGGgraph")
library("NCIgraph")
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

