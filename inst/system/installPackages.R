## Bioconductor packages
if (!requireNamespace("BiocManager", quietly=TRUE))
    install.packages("BiocManager")
BiocManager::install("Rgraphviz")
BiocManager::install("KEGGgraph")
BiocManager::install("marray")

## CRAN packages
install.packages("rrcov")
install.packages("corpcor")
install.packages("RBGL")
install.packages("fields")

## HB's packages
source("http://www.braju.com/R/hbLite.R")
hbLite("R.utils")

