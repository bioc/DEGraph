## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
## Test all (locally stored) pathways from KEGG on Loi's data
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

library('R.utils')
library('xtable')

verbose <- Arguments$getVerbose(-8, timestamp=TRUE)

## force <- TRUE

path <- system.file("downloadScripts", package="DEGraph")
sourceDirectory(path)

path <- system.file("demoScripts", package="DEGraph")
sourceDirectory(path)

## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
## get all NCI pathways
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

##library('NCIgraph')
##path <- system.file("downloadScripts", package="NCIgraph")
##sourceDirectory(path)
##load(file.path('rawNCINetworks','NCI-cyList.RData'))

##grList <- getNCIPathways(cyList=NCI.cyList, verbose=verbose)$pList

## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
## get all KEGG pathways
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

grList <- getKEGGPathways(organism="hsa", metaTag="non-metabolic", verbose=verbose)

## grList <- getKEGGPathways(organism="hsa", metaTag="non-metabolic", patt="04060", verbose=verbose)

## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
## test all KEGG pathways
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

## Individual t-test p-values
res <- apply(exprData, 1, FUN=function(x) {
 tt <- t.test(x[classData==0], x[classData==1])  ## resistant (0) vs sensitive (1)
 c(statistic=tt$statistic[[1]], p.value=tt$p.value)
})
ttpv <- res["p.value", ]
tts <- res["statistic", ]

if(!is.NCIgraph(grList[[1]]))
  {
    names(ttpv) <- translateGeneID2KEGGID(names(ttpv))
    names(tts) <- translateGeneID2KEGGID(names(tts))
  }

prop <- 0.2 ## proportion of the spectrum to be retained for T2-Fourier 

## Multivariate tests
resList <- NULL
## for (ii in 1:5) {
for (ii in seq(along=grList)) {
  verbose && cat(verbose, ii)
  gr <- grList[[ii]]
  if(min(length(nodes(gr)),length(gr@edgeData@data))>0) {
    res <- testOneGraph(gr, exprData, classData, verbose=verbose, prop=prop)
  } else {
    res <- NULL
  }
  if(!is.null(res)) {
    res <- lapply(res, FUN=function(x) {
      if(!is.null(x)) {
          if(is.NCIgraph(x$graph)) {
            geneIDs <- translateNCI2GeneID(x$graph)
          } else {
            geneIDs <- translateKEGGID2GeneID(nodes(x$graph))
          }
          ht <- hyper.test(ttpv, geneIDs, thr=0.001)
          x$p.value <- c(x$p.value, ht$p.value)
          names(x$p.value)[length(names(x$p.value))] <- "Hypergeometric test"
          return(x)
        } else {
          return(NULL)
        }
    })
  }
  resList <- c(resList, list(res))
}
resNames <- names(grList)

if(is.NCIgraph(grList[[1]])) {
  pLabels <- names(grList)
} else {
  pLabels <- sapply(grList, attr, "label")
}

## get rid of NULL results (no connected component of size > 1)
isNULL <- sapply(resList, is.null)
if (sum(isNULL)) {
  grList[isNULL]
  resList <- resList[!isNULL]
  resNames <- names(grList)[!isNULL]
  pLabels <- pLabels[!isNULL]
}

resL <- sapply(resList, length)
graphNames <- rep(resNames, times=resL)
pathwayNames <- rep(pLabels, times=resL)
verbose && str(verbose, graphNames)
verbose && print(verbose, table(graphNames))

graphList <- NULL
for (res in resList) {
  grl <- lapply(res, FUN=function(x) {
    x$graph
  })
  graphList <- c(graphList, as.list(grl))
}

ndims <- NULL
for (res in resList) {
  ndim <- sapply(res, FUN=function(x) {
    x$k
  })
  ndims <- c(ndims, ndim)
}

## Compare results
pKEGG <- NULL
for (res in resList) {
  pp <- sapply(res, FUN=function(x) {
    x$p.value
  })
  if(length(pp))
    pKEGG <- cbind(pKEGG, pp)
}
colnames(pKEGG) <- graphNames
rn <- rownames(pKEGG)
rownames(pKEGG)[grep("Fourier", rn)] <- paste("T2 (", round(100*prop), "% Fourier components)", sep="")

pH <- pKEGG[1,]
pF <- pKEGG[2,]
pHG <- pKEGG[3,]

## Before multiple testing correction
plot(sort(pKEGG[1, ]), t='l', xlab="Rank", ylab="p-value")
lines(sort(pKEGG[2, ]), col=2)
lines(sort(pKEGG[3, ]), col=3)
legend("topleft", rownames(pKEGG), lty=1, col=1:3) 
     
## multiple testing correction(s)
BHpKEGG <- t(apply(pKEGG, 1, p.adjust, method="BH"))
BYpKEGG <- t(apply(pKEGG, 1, p.adjust, method="BY"))  ## more conservative, but OK w. any dependency

thr <- c(0.05, 0.1, 0.2)[2]

## number of pathways passing the threshold
apply(BHpKEGG, 1, FUN=function(p) sum(p<thr))
apply(BYpKEGG, 1, FUN=function(p) sum(p<thr))

wHF <- which((pH<thr) & (pF>thr))
wFH <- which((pF<thr) & (pH>thr))
wFHG <- which((pF<thr) & (pHG>thr))
wHGF <- which((pHG<thr) & (pF>thr))
wHHG <- which((pH<thr) & (pHG>thr))
wHGH <- which((pHG<thr) & (pH>thr))

dev.new()
pKEGG[is.na(pKEGG)] <- 1
c3 <- t(pKEGG < 0.01)
a <- vennCounts(c3)
vennDiagram(a,names=c("Hotelling","Fourier","Hypergeometric"))
title(sprintf("Number of significant genes\n with the three uncorrected tests at level %3f",0.01))

dev.new()
BHpKEGG[is.na(BHpKEGG)] <- 1
c3 <- t(BHpKEGG < thr)
a <- vennCounts(c3)
vennDiagram(a,names=c("Hotelling","Fourier","Hypergeometric"))
title(sprintf("Number of significant genes\n with the three BH-corrected tests at level %3f",thr))

dev.new()
BYpKEGG[is.na(BYpKEGG)] <- 1
c3 <- t(BYpKEGG < thr)
a <- vennCounts(c3)
vennDiagram(a,names=c("Hotelling","Fourier","Hypergeometric"))
title(sprintf("Number of significant genes\n with the three BY-corrected tests at level %3f",thr))

## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
## Visualization of a pathway
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

##pathwayNames[graphNames[order(pHG)[[1]]]]

oo <- order(pF/pH)[1:25]

gIdx <- oo[1]

gr <- graphList[[gIdx]]
pathwayNames[graphNames[gIdx]]
if (require(marray)) {
  pal <- maPalette(low="red", high="green", mid="black", k=100)
} else {
  pal <- heat.colors(100)
}
shift <- tts # Plot t-statistics
##shift <- -getMeanShift(exprData, classData) # Plot mean shifts

## dn <- getDisplayName(gr, shortLabel=TRUE)

if(is.NCIgraph(gr)) {
  mm <- match(translateNCI2GeneID(gr), rownames(annData))
  nodes(gr) <- translateNCI2GeneID(gr)
} else {
  mm <- match(translateKEGGID2GeneID(nodes(gr)), rownames(annData))
}
dn <- annData[mm, "NCBI.gene.symbol"]

##pdf("Loi2008-leukocyte.pdf",width=10,height=10)
res <- plotValuedGraph(gr, values=shift, nodeLabels=dn, qMax=0.95, colorPalette=pal, height=40, lwd=1, cex=0.7, verbose=verbose)
##devDone()
title(paste(pathwayNames[graphNames[gIdx]],"pF=",as.character(signif(pF[gIdx],4)),"pH=",as.character(signif(pH[gIdx],4))))

## mm <- match(translateKEGGID2GeneID(nodes(gr)), names(shift))
## plotValuedGraph(gr, values=shift[mm], qMax=0.95, colorPalette=pal, height=1, lwd=1, verbose=verbose)

## Example: look at the pathways which are significant with Fourier
## after correction but would have large p-values without the Fourier
## step.

fSignif <- which(BHpKEGG[2,]<0.05)
fSignif <- fSignif[order(BHpKEGG[1,fSignif],decreasing=TRUE)]
##fSignif <- fSignif[order(BHpKEGG[2,fSignif],decreasing=FALSE)] # Order by increasing pF


## Plot all graphs
##dev.new()
##printFileName <- "Tech-MI-smoothGraphs-loi-dec"
printFileName <- "Loi-NCI"
printStarted <- FALSE
pdf(paste(printFileName,".pdf",sep=""))
for (gIdx in fSignif) {
  ## Get graph
  gr <- graphList[[gIdx]]

  if(is.NCIgraph(gr)) {
    mm <- match(translateNCI2GeneID(gr), rownames(annData))
    nodes(gr) <- translateNCI2GeneID(gr)
  } else {
    mm <- match(translateKEGGID2GeneID(nodes(gr)), rownames(annData))
  }
  dn <- annData[mm, "NCBI.gene.symbol"]

  ## Try to plot
  tt <- try(res <- plotValuedGraph(gr, values=shift, nodeLabels=dn, qMax=0.95, colorPalette=pal, height=40, lwd=1, cex=0.3, verbose=verbose))
  ## Add legend
  if (class(tt)!="try-error") {
    stext(side=3, pos=0, pathwayNames[gIdx])
    ps <- signif(pKEGG[, gIdx],2)
    txt1 <- paste("p(T2)=", ps[1], sep="")
    txt2 <- paste("p(T2F[", ndims[gIdx], "])=", ps[2], sep="")
    txt <- paste(txt1, txt2, sep="\n")
    stext(side=3, pos=1, txt)
    image.plot(legend.only=TRUE, zlim=range(res$breaks), col=pal, legend.shrink=0.3, legend.width=0.8, legend.lab="t-scores", legend.mar=3.3)
    ## Now print gene table in a tex file
    cIdx <- match(nodes(gr),names(shift)) #which(names(shift)%in% nodes(gr))
    gTable <- cbind(dn,signif(shift[cIdx],2),signif(ttpv[cIdx],2))
    colnames(gTable) <- c("Gene symbol","t-stat","p-value")
    xt <- xtable(gTable,caption=sprintf("%s, p(Hotelling)=%g, p(netHotelling)=%g",pathwayNames[gIdx], ps[1], ps[2]))
    print(xt,file=paste(printFileName,"-Tables.tex",sep=""),tabular.environment='longtable',floating=FALSE, append=printStarted)
    printStarted <- TRUE
  } else {
    warning("check gIdx=", gIdx)
  }
}
devDone()

## Save results
ts <- format(Sys.time(), "%Y-%m-%d,%X")
filename <- sprintf("pKEGG,signed,normalized,%s.rda", ts)
save(resList, pKEGG, graphList, pathwayNames, ndims, file=filename)
