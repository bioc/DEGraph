\name{annLoi2008}
\docType{data}
\alias{annLoi2008}
\title{Annotation data used in the DEGraph package vignette}
\description{This data set gives NCBI, Hugo and alternative gene symbols
  along with the cytoband and description for the 227 genes used in the
  DEGraph package vignette. This comes from the 15737 gene, 255 patient
  dataset of Loi et al. (2008) which was used to study resistance to
  tamoxifen treatment in hormone-dependent breast cancer.  }
\usage{annLoi2008}
\format{A matrix of 227 lines and 5 columns.}
\examples{
data("Loi2008_DEGraphVignette")

dim(annLoi2008)
head(annLoi2008)
}
\source{
  Loi et al.,
  \emph{Predicting prognosis using molecular profiling in estrogen
  receptor-positive breast cancer treated with tamoxifen.}
  BMC Genomics, 9(1):239, 2008.}
\references{
  Loi et al.,
  \emph{Predicting prognosis using molecular profiling in estrogen
  receptor-positive breast cancer treated with tamoxifen.}
  BMC Genomics, 9(1):239, 2008.
}
\author{Laurent Jacob, Pierre Neuvial and Sandrine Dudoit}
\keyword{datasets}
