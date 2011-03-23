\name{exprLoi2008}
\docType{data}
\alias{exprLoi2008}
\title{Gene expression data used in the DEGraph package vignette}
\description{
  This data set gives gene expression data for a subset of 227
  genes used in the DEGraph package vignette. This comes from the 15737
  gene, 255 patient dataset of Loi et al. (2008) which was used to study
  resistance to tamoxifen treatment in hormone-dependent breast cancer.
}
\usage{exprLoi2008}
\format{A matrix of 227 lines and 255 columns.}
\examples{
data("Loi2008_DEGraphVignette")

dim(exprLoi2008)
head(exprLoi2008)
}

\details{
  The original data set corresponds to data processed by RMA and median-centered
  as available from the GSE6532 GEO archive:
  http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE6532.
  
  These data were summarized from the probe set level to the gene level as 
  follows. The expression level of a gene was defined as the expression level 
  of the probe set with largest alignment score among all probe sets mapping 
  to this gene according to the annotation in GSE6532. When the largest 
  alignment score was achieved by several probe sets, the median expression 
  level of those probe sets was taken.
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
