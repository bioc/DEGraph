library(R.utils)
verbose <- Arguments$getVerbose(-8, timestamp=TRUE)

## setup
srcRootPath <- "http://www.stat.berkeley.edu/~pierre/share"

dataSet <- "Loi08"
chipType <- "HG-U133A+B"

ttTags <- c("tam", "untreated")
outTags <- c("metastasis", "relapse")

if (!exists("force", mode="character")) {
  force <- FALSE
}

verbose && enter(verbose, "Downloading expression data")
path <- file.path("exprData", dataSet, chipType)

srcPath <- file.path(srcRootPath, path)
tgtPath <- path
tgtPath <- Arguments$getWritablePath(tgtPath)

for (ttTag in ttTags) {
  filename <- sprintf("%s,%s,geneLevel,median.xdr", dataSet, ttTag)
  srcPathname <- file.path(srcPath, filename)
  tgtPathname <- file.path(tgtPath, filename)
  
  pn <- downloadFile(srcPathname, tgtPathname, overwrite=force, verbose=verbose)
}
verbose && exit(verbose)

## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
## Annotation data (classes)
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

verbose && enter(verbose, "Downloading annotation data (classes of patients)")
path <- file.path("annotationData/samples", dataSet)

srcPath <- file.path(srcRootPath, path)
tgtPath <- path
tgtPath <- Arguments$getWritablePath(tgtPath)

for (ttTag in ttTags) {
  for (outTag in outTags) {
    filename <- sprintf("%s,%s,%s.xdr", dataSet, outTag, ttTag)
    srcPathname <- file.path(srcPath, filename)
    tgtPathname <- file.path(tgtPath, filename)
  
    pn <- downloadFile(srcPathname, tgtPathname, skip=!force, overwrite=force, verbose=verbose)
  }
}
verbose && exit(verbose)

## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
## Annotation data (classes)
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

verbose && enter(verbose, "Downloading gene and probe annotation data")
path <- file.path("annotationData/chipTypes", chipType)

srcPath <- file.path(srcRootPath, path)
tgtPath <- path
tgtPath <- Arguments$getWritablePath(tgtPath)

pTags <- c("probes", "genes")
for (pTag in pTags) {
  filename <- sprintf("%s,%s,from%s.xdr", chipType, pTag, dataSet)
  srcPathname <- file.path(srcPath, filename)
  tgtPathname <- file.path(tgtPath, filename)
  
  pn <- downloadFile(srcPathname, tgtPathname, skip=!force, overwrite=force, verbose=verbose)
}
verbose && exit(verbose)

############################################################################
# HISTORY:
# 2010-09-17
# o Added gene and probe annotation data.
# 2010-09-14
# o Created.
############################################################################
