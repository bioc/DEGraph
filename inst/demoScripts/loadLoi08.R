## - - - - - - - - - - - - - - - - - - -
## Setup
## - - - - - - - - - - - - - - - - - - -

## exprData/
##   Loi08/
##     HG-U133A+B/
##       Loi08,tam,geneLevel,median.xdr
##       Loi08,untreated,geneLevel,median.xdr
##
## annotationData/
##   samples/
##     Loi08/
##       Loi08,metastasis,tam.xdr
##       Loi08,metastasis,untreated.xdr
##       Loi08,relapse,tam.xdr
##       Loi08,relapse,untreated.xdr

dataSet <- "Loi08"
chipType <- "HG-U133A+B"

## - - - - - - - - - - - - - - - - - - -
## Setup
## - - - - - - - - - - - - - - - - - - -
lev <- "gene"                                 ## gene-level expression data
treat <- c("tam", "untreated")[1]             ## Tamoxifen treated patients
outcome <- c("relapse", "metastasis")[2]       ## Metastasis/non metastasis/NA

## - - - - - - - - - - - - - - - - - - -
## Load gene expression data
## - - - - - - - - - - - - - - - - - - -

expPath <- "exprData"
expPath <- Arguments$getReadablePath(expPath)

path <- file.path(expPath, dataSet, chipType)
path <- Arguments$getReadablePath(path)

pattern <- sprintf("%s,%s,%s", dataSet, treat, lev)
filename <- list.files(path, pattern=pattern)
pathname <- file.path(path, filename)

exprData <- loadObject(pathname)

## - - - - - - - - - - - - - - - - - - - - - - - - - -
## Load class information
## - - - - - - - - - - - - - - - - - - - - - - - - - -
annPath <- "annotationData"
annPath <- Arguments$getWritablePath(annPath)

path <- file.path(annPath, "samples", dataSet)
path <- Arguments$getReadablePath(path)

pattern <- sprintf("%s,%s,%s", dataSet, outcome, treat)
filename <- list.files(path, pattern=pattern)
pathname <- file.path(path, filename)

classData <- loadObject(pathname)
## sanity check
stopifnot(names(classData)==colnames(exprData))

## - - - - - - - - - - - - - - - - - - - - - - - - - -
## Load gene annotation information
## - - - - - - - - - - - - - - - - - - - - - - - - - -
annPath <- "annotationData"
annPath <- Arguments$getWritablePath(annPath)

path <- file.path(annPath, "chipTypes", chipType)
path <- Arguments$getReadablePath(path)
filename <- sprintf("%s,genes,from%s.xdr", chipType, dataSet)
pathname <- file.path(path, filename)

annData <- loadObject(pathname)
## sanity check
stopifnot(all.equal(rownames(annData), rownames(exprData)))

## discard NA:s
ww <- which(!is.na(classData))
exprData <- exprData[, ww]
classData <- classData[ww]

verbose && cat(verbose, "exprData:")
verbose && str(verbose, exprData)

verbose && cat(verbose, "classData:")
verbose && str(verbose, classData)

verbose && cat(verbose, "annData:")
verbose && str(verbose, annData)

############################################################################
# HISTORY:
# 2010-09-17
# o Added gene annotation data.
# 2010-09-13
# o Created.
############################################################################
