## Libraries
library(R.utils)
library(R.menu) ## source("http://aroma-project.org/hbLite.R") hbLite("R.menu")
## library(aroma.core) ## for 'stext'.... TODO: remove this dependency

require(graph)
require(mvtnorm) ## for 'rmvnorm'
require(rrcov) ## for 'T2.test'
require(corpcor)

require(KEGGgraph)
require(KEGG.db)
require(Rgraphviz)
require(RBGL)

library(fields) # For image.plot called in plotValuedGraph
## library(limma) # For venn diagrams

require(lattice)
require(marray)
require(RColorBrewer)

## library(XML) ## ?

verbose <- Arguments$getVerbose(-8, timestamp=TRUE)

