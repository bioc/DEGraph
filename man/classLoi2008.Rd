\name{classLoi2008}
\docType{data}
\alias{classLoi2008}
\title{Tamoxifen treatment resistance status data used in the DEGraph
  package vignette}

\description{ This data set gives resistance status data for the 255
  patients used in the DEGraph package vignette. This comes from the
  15737 gene, 255 patient dataset of Loi et al. (2008) which was used to
  study resistance to tamoxifen treatment in hormone-dependent breast
  cancer.  }
\usage{classLoi2008}
\format{A vector of 255 elements which are either 0 (resistance to
  treatment) or 1 (sensitivity to treatment).}
\examples{
data("Loi2008_DEGraphVignette")

dim(classLoi2008)
head(classLoi2008)
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
