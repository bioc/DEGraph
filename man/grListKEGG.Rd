\name{grListKEGG}
\docType{data}
\alias{grListKEGG}
\title{KEGG networks used in the DEGraph package vignette}
\description{ This data set gives KEGGgraph objects for two KEGG
  non-metabolic pathways ("Natural killer cell mediated cytotoxicity"
  and "Insulin signaling pathway").}
\usage{grListKEGG}  
\format{A list of two elements.}
\examples{
library("Rgraphviz")

data("Loi2008_DEGraphVignette")

grListKEGG
plot(grListKEGG[[1]])
}
\author{Laurent Jacob, Pierre Neuvial and Sandrine Dudoit}
\keyword{datasets}
