\name{WGD.matrix.tranche}
\alias{WGD.matrix.tranche}
\title{
Apply tranche filter to WGD matrix
}
\description{
Subsets a WGD matrix to retain the innermost N percent of bins, based on
interquartile range. Designed to protect against pericentromeric and
peritelomeric outliers highly enriched for reads due to alpha satelites
as well as N-masked or poorly mapping bins not correctly blacklisted during
coverage matrix generation.
}
\usage{
WGD.matrix.tranche(mat, tranche = 0.995, quiet = F)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{mat}{
  Matrix to filter. Must be imported previously with WGD.readmatrix.
}
  \item{tranche}{
  Fraction of bins to keep; e.g. 0.995 = keep middle-most 99.5\% of bins.
}
  \item{quiet}{
  Disable verbose output.
}
}
\value{
\item{mat}{Merged coverage matrix}
\item{res}{Matrix of coverage residuals}
\item{stat}{Binwise summary statistics for coverage matrix}
\item{rstat}{Binwise summary statistics for coverage residuals}
}
\seealso{
\url{https://github.com/RCollins13/WGD}
}
\author{
Ryan L. Collins <rlcollins@g.harvard.edu>
}
\keyword{ WGD }
\keyword{ Matrix }
\keyword{ Filter }
\keyword{ Tranche }
