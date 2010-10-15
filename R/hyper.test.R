## Copyright 2010 Laurent Jacob, Pierre Neuvial and Sandrine Dudoit.

## This file is part of DEGraph.

## DEGraph is free software: you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.

## DEGraph is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.

## You should have received a copy of the GNU General Public License
## along with DEGraph.  If not, see <http://www.gnu.org/licenses/>.

#########################################################################/**
## @RdocFunction hyper.test
##
## @title "Performs an hypergeometric test of enrichment of a set of
##   hypotheses in significant elements"
##
## \description{
##  @get "title".
## }
##
## @synopsis
##
## \arguments{
##   \item{p.values}{A named @numeric vector giving the p-values of all
##     tested elements.}
##   \item{testSet}{A @character vector giving the ids of the elements in the
##     tested set. Elements of 'testSet' must have a match in 'names(p.values)'.}
##   \item{thr}{A @numeric value between 0 and 1 giving the threshold on
##     p-values at which an element is declared to be significant.}
##   \item{universe}{An @integer value giving the number of elelments in the
##     considered universe.  Defaults to 'length(p.values)'.}
##   \item{verbose}{If @TRUE, extra information is output.}
## }
##
## \value{
##  A @list with class "htest" containing the following components:
##  \describe{
##    \item{statistic}{A @numeric value, the test statistic.}
##    \item{p.value}{A @numeric value, the corresponding p-value.}
##   }
## }
##
## @author
##
## \seealso{
##   @see "AN.test"
##   @see "BS.test"
##   @see "graph.T2.test"
## }
##
## @examples "../incl/tests.Rex"
##
##*/########################################################################

hyper.test <- function(p.values, testSet, thr=0.001, universe=length(p.values), verbose=FALSE) {
  ## - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ## Validate arguments
  ## - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ## Argument 'p.values'
  p.values <- Arguments$getNumerics(p.values)
  pnames <- names(p.values)
  pnames <- Arguments$getCharacters(pnames)
  if (length(pnames)==0) {
    throw("'names(p.values)' should not be NULL for a hypergeometric test to be performed.")
  }

  ## Argument 'testSet'
  testSet <- Arguments$getCharacters(testSet)
  if (length(testSet)==0) {
    throw("'testSet' should have at least one element for a hypergeometric test to be performed.")
  }
  ## Argument 'thr'
  thr <- Arguments$getNumeric(thr)

  ## Argument 'universe'
  universe <- Arguments$getNumeric(universe)

  ## Argument 'verbose'
  verbose <- Arguments$getVerbose(verbose)
  if (verbose) {
    cat <- R.utils::cat
    pushState(verbose)
    on.exit(popState(verbose))
  } 

  ## - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ## Hypergeometric testing
  ## - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  mm <- match(testSet, pnames)
  isNA <- is.na(mm)
  idxs <- which(isNA)
  if (length(idxs)) {
    verbose && cat(verbose, "Elements of 'testSet' not found in 'p.values':")
    verbose && str(verbose, pnames[idxs])
  }
  idxs <- which(!isNA)
  mm <- mm[idxs]
  p <- length(mm)
  if (p==0) {
    warning("No elements of 'testSet' found in 'p.values':  returning NA.")
    qHG <- NA
    pHG <- NA
  } else {
    qHG <- sum(p.values[mm] < thr) # How many significant elements in the set
    mHG <- p # How many elements in the set
    nHG <- universe - p # How many elements outside the set
    kHG <- sum(p.values < thr)# How many significant elements in total
    pHG <- 1-phyper(qHG, mHG, nHG, kHG)
  }
  res <- list(statistic=qHG,  p.value=pHG)
  class(res) <- "htest"
  res
}

############################################################################
## HISTORY
## 2010-10-08
## o Now validating argument 'verbose'.
## 2010-09-25
## o Arguments are now validated.
## o Now returning NA when no elements of 'testSet' is in 'p.values'.
## o Now throwing an error when 'names(p.values)' is NULL.
## 2010-09-23
## o Added an example.
## o Made more generic by replacing 'gene' by 'element', 'DE' by
##   'significant', and 'gene set' by 'set'.
## o Now returns an object of class "htest".
############################################################################
