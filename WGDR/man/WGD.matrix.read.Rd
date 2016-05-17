\name{WGD.matrix.read}
\alias{WGD.matrix.read}
\title{
Import WGD coverage matrix
}
\description{
Standardized import & processing of WGD coverage matrix
}
\usage{
WGD.matrix.read(path, allosomes = F, quiet = F)
}
\arguments{
  \item{path}{
  Full path to input matrix of normalized coverage. Must be in standardized WGD
  matrix format; see documentation for makeMatrix.sh in main WGD docs for
  details. Input can be gzipped.
}
  \item{allosomes}{
  Option to include non-numeric contigs, such as chrX and chrY.
}
  \item{quiet}{
  Disable verbose output.
}
}
\details{
\code{path} can specify a gzipped matrix. Input will be decompressed as it is
imported into memory.
}
\value{
\item{mat}{Input coverage matrix}
\item{res}{Matrix of coverage residuals}
\item{stat}{Binwise summary statistics for coverage matrix}
\item{rstat}{Binwise summary statistics for coverage residuals}
\item{rstat}{Per-sample summary statistics for coverage residuals}
}
\seealso{
\url{https://github.com/RCollins13/WGD}
}
\author{
Ryan L. Collins <rlcollins@g.harvard.edu>
}
\keyword{ WGD }
\keyword{ Matrix }
\keyword{ Import }
