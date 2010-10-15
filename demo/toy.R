library("R.utils")
library("DEGraph")
pathname <- system.file("demoScripts", "001.setup.R", package="DEGraph")
source(pathname)

## Parameters

NP <- 500
NN <- 500
k <- 3  ## true number of Fourier components
kRs <- 2:4  ## number of Fourier components retained
n1 <- 20
n2 <- 20
# shift <- 0.5 ## Amplitude of the shift
ncp <- 0.5

print("[toy] Building data")

## Get graph
nnodes <- 20
nedges <- 20
G <- randomWAMGraph(nnodes=nnodes,nedges=nedges)
A <- G@adjMat

#A <- rbind(c(0,1,0,0,1,0),
#           c(1,0,1,0,1,0),
#           c(0,1,0,1,0,0),
#           c(0,0,1,0,1,1),
#           c(1,1,0,1,0,0),
#           c(0,0,0,1,0,0))

p <- nrow(A)
##D <- diag(rowSums(abs(A)))
##L <- D - A # Unormalized version
iDs <- diag(1/sqrt(rowSums(abs(A))))
L <- diag(rep(1,p)) - (iDs %*% A %*% iDs)
##edL <- eigen(L,symmetric=TRUE)
##U <- fliplr(edL$vectors)
##l <- fliplr(edL$values)

lfA <- laplacianFromA(A,ltype="unnormalized")
U <- lfA$U
l <- lfA$l

sigma <- diag(p)/sqrt(p)
sigma[1:k,1:k] <- sigma[1:k,1:k] + 0.5/sqrt(p)
diag(sigma)[1:k] <- diag(sigma)[1:k] - 0.6/sqrt(p)

## Build and test examples

scores = rep(NA,NN+NP)
uscores = matrix(NA, NN+NP, length(kRs))
pcascores = rep(NA,NN+NP)
anscores = rep(NA,NN+NP)
ant = rep(NA,NN+NP)
ank = rep(NA,NN+NP)
cqscores = rep(NA,NN+NP)
cqq = rep(NA,NN+NP)
bsscores = rep(NA,NN+NP)
bsz = rep(NA,NN+NP)

print("[toy] Starting tests")

for(ii in 1:(NP+NN))
{
  ##X <- twoSmoothSampleFromGraph(n1,n2,shiftM2=ncp*as(ii<NP,"integer"),sigma,U=U,l=l)
  X <- twoSampleFromGraph(n1,n2,shiftM2=ncp*as(ii<NP,"integer"),sigma,U=U,k=k)
  ## Raw T-square
  t <- T2.test(X$X1,X$X2)
  scores[ii] <- t$p.value
  ## Filtered T-squares
  for (kk in seq(along=kRs)) {
    kR <- kRs[kk]
    tu <- graph.T2.test(X$X1,X$X2,lfA=lfA,k=kR)
    ##tu <- T2.test(X$X1%*%U[,1:kR],X$X2%*%U[,1:kR])
    uscores[ii, kk] <- tu$p.value
  }
  ## PCA T-square
  edX <- svd(rbind(X$X1,X$X2))
  E <- edX$v
  tpca <- T2.test(X$X1%*%E[,1:k],X$X2%*%E[,1:k])
  pcascores[ii] <- tpca$p.value
  ## Adaptive Neyman (Fan)
  tan <- AN.test(X$X1%*%U,X$X2%*%U,p)
  anscores[ii] <- tan$p.value
  ant[ii] <- tan$statistic
  ank[ii] <- tan$kstar
  ## Chen & Qin
  ##tcq <- CQ.Test(X$X1,X$X2)
  ##cqscores[ii] <- tcq$p.value
  ##cqq[ii] <- tcq$q
  ## Chen & Qin
  tbs <- BS.test(X$X1,X$X2)
  bsscores[ii] <- tbs$p.value
  bsz[ii] <- tbs$statistic
}

print("[toy] Tests done")

## ROC


# Empirical
thr = seq(0,1,0.001) # c(0.00001,0.0001,0.001,0.01,0.1)
nT <- length(thr)
tpower = rep(NA, nT)
tlevel = rep(NA, nT)
utpower = matrix(NA, nT, length(kRs))
utlevel = matrix(NA, nT, length(kRs))
pcatpower = rep(NA, nT)
pcatlevel = rep(NA, nT)
anpower = rep(NA, nT)
anlevel = rep(NA, nT)
cqpower = rep(NA, nT)
cqlevel = rep(NA, nT)
bspower = rep(NA, nT)
bslevel = rep(NA, nT)

for(cc in seq(along=thr))
{
  i <- thr[cc]
  tpower[cc] <- mean(scores[1:NP]<i)
  tlevel[cc] <- mean(scores[NP+1:NN]<i)
  for (kk in seq(along=kRs)) {
    utpower[cc, kk] <- mean(uscores[1:NP, kk]<i)
    utlevel[cc, kk] <- mean(uscores[NP+1:NN, kk]<i)
  }
  pcatpower[cc] <- mean(pcascores[1:NP]<i)
  pcatlevel[cc] <- mean(pcascores[NP+1:NN]<i)
  anpower[cc] <- mean(anscores[1:NP]<i)
  anlevel[cc] <- mean(anscores[NP+1:NN]<i)
  cqpower[cc] <- mean(cqscores[1:NP]<i)
  cqlevel[cc] <- mean(cqscores[NP+1:NN]<i)
  bspower[cc] <- mean(bsscores[1:NP]<i)
  bslevel[cc] <- mean(bsscores[NP+1:NN]<i)
}


# Theoretical
kR <- k
#ncp <- (n1*n2/(n1+n2)) * shift^2*sqrt(p)
ncp <- (n1*n2/(n1+n2)) * ncp
df1 <- p
df2 <- n1+n2-p-1
nullquant <- qf(1-thr, df1, df2, 0)
thtpower <- 1-pf(nullquant, p, df2, ncp)

## thutpower <- sapply(kRs, FUN=function(kR) {
##   df1 <- kR
##   df2 <- n1+n2-kR-1
##   unullquant <- qf(1-thr, df1, df2, 0)
##   df1 <- k
##   df2 <- n1+n2-k-1
##   1-pf(unullquant, df1, df2, ncp)
## })

thutpower <- sapply(kRs, FUN=function(kR) {
  if (kR==k) {
    df1 <- k
    df2 <- n1+n2-k-1
    unullquant <- qf(1-thr, df1, df2, 0)
    res <- 1-pf(unullquant, df1, df2, ncp)
  } else if (kR>k) {
    df1 <- kR
    df2 <- n1+n2-kR-1
    unullquant <- qf(1-thr, df1, df2, 0)
    res <- 1-pf(unullquant, df1, df2, ncp)
  } else {
    res <- rep(NA, length(thr))
  }
  res
})

## plot parameters
lim <- c(0,1)

txt <- paste("Graph of size ", nnodes, " Shift in the first ", k, " Fourier components", sep="")

lgd <- paste("T2, k=", kRs, sep="")
lgd <- c("T2", lgd)
lgd <- c(paste("theoretical", lgd), paste("simulated", lgd))
lgd <- c(lgd, paste("simulated", k, "first PC"), "simulated AN")

nn <- length(kRs)+1
cols <- seq(nn+2)
ltys <- c(rep(1:2, each=nn), 3, 3)

## Plot: power vs threshold
if(FALSE){
x11()
plot(thr,tpower,type='s', col=cols[1], xlab="Threshold", ylab="Test power", xlim=lim, ylim=lim, lty=2)
lines(thr,thtpower, col=cols[1])
for (kk in seq(along=kRs)) {
  lines(thr, utpower[, kk], col=cols[1+kk], lty=2, type='s')
  lines(thr,thutpower[, kk], col=cols[1+kk])
}
lines(thr, pcatpower, col=cols[2+kk], lty=3, type='s')
lines(thr, anpower, col=cols[3+kk], lty=3, type='s')
legend("bottomright", lgd, col=cols, lty=ltys, ncol=2)
mtext(txt, side=3, adj=1)
}

if (FALSE) {
## Plot: level vs threshold
x11()
plot(thr,tlevel,type='s', col=cols[1], xlab="Threshold", ylab="Test level", xlim=lim, ylim=lim)
for (kk in seq(along=kRs)) {
  lines(thr, utlevel[, kk], col=cols[1+kk], lty=1, type='s')
}
lines(thr, pcatlevel, col=cols[1+kk]+1, lty=3, type='s')
lines(thr, anlevel, col=cols[1+kk]+2, lty=3, type='s')
lines(thr, bslevel, col=cols[1+kk]+3, lty=3, type='s')
legend("bottomright", lgd[nn+(1:(nn+2))], col=cols, lty=c(ltys[1:nn],3,3))
mtext(txt, side=3, adj=1)
}

if(TRUE){
## Plot: power vs level
x11()
plot(tlevel,tpower,type='s', col=cols[1], xlab="Test level", ylab="Test power", xlim=lim, ylim=lim)
for (kk in seq(along=kRs)) {
  lines(utlevel[, kk], utpower[, kk], col=cols[1+kk], lty=1, type='s')
}
lines(pcatlevel, pcatpower, col=nn+1, lty=3, type='s')
lines(anlevel, anpower, col=nn+2, lty=3, type='s')
##lines(cqlevel, cqpower, col=nn+2, lty=4, type='s')
lines(bslevel, bspower, col=nn+2, lty=4, type='s')
legend("bottomright", lgd[nn+(1:(nn+2))], col=c(cols,nn+2), lty=c(ltys[1:nn],3,3))
mtext(txt, side=3, adj=1)
}

## Plot: theoretical power vs threshold
x11()
plot(thr,thtpower,type='l', col=cols[1], xlab="Threshold", ylab="Theoretical test power", xlim=lim, ylim=lim)
for (kk in seq(along=kRs)) {
  lines(thr,thutpower[, kk], col=cols[1+kk])
}
legend("bottomright", lgd[1:nn], col=cols, lty=1)
mtext(txt, side=3, adj=1)
