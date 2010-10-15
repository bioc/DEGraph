%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  ./getKEGGPathways.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{getKEGGPathways}
\alias{getKEGGPathways}


\title{Builds a graph for each of the KEGG pathways}

\description{
 Builds a graph for each of the KEGG pathways.
}

\usage{getKEGGPathways(path=NULL, rootPath="networkData/ftp.genome.jp/pub/kegg/xml/kgml", organism="hsa", metaTag=c("non-metabolic", "metabolic"), pattern=NULL, verbose=FALSE)}

\arguments{
  \item{path}{A \code{\link[base]{character}} value, the local _full_ path of KGML data.}
  \item{rootPath}{A \code{\link[base]{character}} value, the local _root_ path of KGML data.}
  \item{organism}{A \code{\link[base]{character}} value specifying the organism whose
     pathways should be considered. Defaults to "hsa" (Homo Sapiens).}
  \item{metaTag}{A \code{\link[base]{character}} value, specifying the type of pathways to
    be considered ("metabolic" or "non-metabolic"). Defaults to "non-metabolic".}
  \item{pattern}{An optional \code{\link[base]{character}} value specifying a file name pattern to
    look for.}
  \item{verbose}{If \code{\link[base:logical]{TRUE}}, extra information is output.}
}

\value{
  A \code{\link[base]{list}} containing a \code{\link[=graph-class]{graph}} object for each KEGG pathway with at least one edge.
}

\details{If 'path' is supplied, KGML files in this directory are loaded.
  Otherwise, KGML files are assumed to be in
  <rootPath>/<metaTag>/"organisms"/<organism>, which mirrors the
 structure of the KEGG KGML file repository.
}

\author{Laurent Jacob, Pierre Neuvial and Sandrine Dudoit}

\seealso{
  \code{\link{parseKGML}}
  \code{\link{KEGGpathway2Graph}}
}

\examples{
library("KEGGgraph")

## example of KGML files
path <- system.file("extdata", package="KEGGgraph")
grList <- getKEGGPathways(path=path, verbose=TRUE)
print(grList)

graph <- grList[[1]]
plotKEGGgraph(graph)

\dontrun{
  ## Download all human KGML pathways locally
  pathname <- system.file("downloadScripts", "downloadKeggXmlFiles.R", package="DEGraph")
  source(pathname)

  ## Load some of them
  grList <- getKEGGPathways(pattern="040", verbose=TRUE)
  print(grList)

  graph <- grList[[1]]
  plotKEGGgraph(graph)
}
}

