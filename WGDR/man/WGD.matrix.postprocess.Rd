\name{WGD.matrix.postprocess}
\alias{WGD.matrix.postprocess}
\title{
Collect WGD metadata
}
\description{
Collects binwise summary statistics, binwise residuals, and binwise
summary statistics for a WGD matrix.
}
\usage{
WGD.matrix.postprocess(mat, allosomes = F)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{mat}{
  Matrix to filter. Must be imported previously with WGD.readmatrix.
}
  \item{allosomes}{
  Option to include non-numeric contigs, such as chrX and chrY.
}
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
\keyword{ Summarize }
\keyword{ Postprocess }
