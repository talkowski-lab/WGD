\name{WGD.matrix.merge}
\alias{WGD.matrix.merge}
\title{
Merge WGD coverage matrices
}
\description{
Standardized merger & processing of two or more WGD coverage matrices
}
\usage{
WGD.matrix.merge(matrices, innerjoin = T, allosomes = F, quiet = F)
}
\arguments{
  \item{matrices}{
  List of matrices to merge. Each matrix must be an R object imported
  previously with WGD.readmatrix.
}
  \item{innerjoin}{
  Option to restrict merger only to bins appearing in all matrices.
}
  \item{allosomes}{
  Option to include non-numeric contigs, such as chrX and chrY.
}
  \item{quiet}{
  Disable verbose output.
}
}
\details{
Does not have contengency for uneven binning between matrices to merge.
}
\value{
\item{mat}{Merged coverage matrix}
\item{res}{Matrix of coverage residuals}
\item{stat}{Binwise summary statistics for coverage matrix}
\item{rstat}{Binwise summary statistics for coverage residuals}
\item{sstat}{Per-sample summary statistics for coverage residuals}
}
\seealso{
\url{https://github.com/RCollins13/WGD}
}
\author{
Ryan L. Collins <rlcollins@g.harvard.edu>
}
\keyword{ WGD }
\keyword{ Matrix }
\keyword{ Merge }
