## Bioconductor packages
source("http://bioconductor.org/biocLite.R")
biocLite("KEGGgraph")
biocLite("marray")

## CRAN packages
install.packages("rrcov")
install.packages("corpcor")
install.packages("RBGL")
install.packages("fields")

## HB's packages
source("http://www.braju.com/R/hbLite.R")
hbLite("R.utils")
## temporary:
installPackages("http://www.braju.com/R/repos/R.oo_1.7.4.tar.gz")

