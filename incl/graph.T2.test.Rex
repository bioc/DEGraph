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