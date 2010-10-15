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

