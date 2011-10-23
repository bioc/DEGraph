library("Rgraphviz")
library("KEGGgraph")
## library("NCIgraph")

data("Loi2008_DEGraphVignette")
exprData <- exprLoi2008
classData <- classLoi2008
annData <- annLoi2008

rn <- rownames(exprData)

## Retrieve expression levels data for genes from one KEGG pathway
graph <- grListKEGG[[1]]
pname <- attr(graph, "label")
print(pname)

## DEGraph T2 test
resList <- testOneGraph(graph, exprData, classData, verbose=TRUE, prop=0.2)

## Largest connected component
res <- resList[[1]]
gr <- res$graph

## individual t statistics
shift <- apply(exprData, 1, FUN=function(x) {
  tt <- t.test(x[classData==0], x[classData==1])
  tt$statistic
})
names(shift) <- translateGeneID2KEGGID(names(shift))

## color palette
if (require(marray)) {
  pal <- maPalette(low="red", high="green", mid="black", k=100)
} else {
  pal <- heat.colors(100)
}

## plot results
dn <- getDisplayName(gr, shortLabel=TRUE)
mm <- match(translateKEGGID2GeneID(nodes(gr)), rownames(annData))
dn <- annData[mm, "NCBI.gene.symbol"]
  
pvg <- plotValuedGraph(gr, values=shift, nodeLabels=dn, qMax=0.95, colorPalette=pal, height=40, lwd=1, verbose=TRUE, cex=0.5)
title(pname)

txt1 <- sprintf("p(T2)=%s", signif(res$p.value[1], 2))
txt2 <- sprintf("p(T2F[%s])=%s", res$k, signif(res$p.value[2]))
txt <- paste(txt1, txt2, sep="\n")
stext(side=3, pos=1, txt)
if (require(fields)) {
  image.plot(legend.only=TRUE, zlim=range(pvg$breaks), col=pal, legend.shrink=0.3, legend.width=0.8, legend.lab="t-scores", legend.mar=3.3)
}
