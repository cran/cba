\name{lminter}
\alias{lminter}
\title{Interpolating Logical Matrices}
\description{
  Interpolate a logical matrix into a lower-resolution representation.
}
\usage{
lminter(x, block.size = 1, nbin = 0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{a logical matrix}
  \item{block.size}{the interpolation block size}
  \item{nbin}{the number of density bins}
}
\details{
  Partitions a binary matrix into square blocks of specified size 
  (length) and interpolates the number of \code{TRUE} values per block 
  using the specified number of bins. 
  
  Note that the effective number of bins is one greater than the specified 
  number because the zero bin is always included. Excess rows and columns 
  at the lower or right margins of the matrix are ignored.

  If the number of bins is null counts are mapped to zero and one 
  thresholding at half of the number of distinct count values including 
  zero. Thus, for even numbered block sizes there is a bias towards zero.
}
\value{
  An integer matrix of bin numbers.
}
\author{Christian Buchta}
\note{Package internal function.}
\seealso{\code{\link{lmplot}} for plotting logical matrices }
\examples{
\dontrun{
x <- matrix(sample(c(FALSE, TRUE), 4 ,rep=TRUE), ncol=2)
lminter(x, block.size=2, nbin=2)
}}
\keyword{cluster}
\keyword{hplot}

