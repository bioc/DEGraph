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
