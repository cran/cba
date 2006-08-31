\name{fitted.proximus}
\alias{fitted.proximus}
\title{Extract from a Proximus Object}
\description{
  Get the full storage representation of the approximated matrix 
  and the pattern labels of the original data samples from an object
  of class \code{proximus}.
}
\usage{
\method{fitted}{proximus}(object, drop = FALSE, ...)
}
\arguments{
  \item{object}{an object of class \code{proximus}.}
  \item{drop}{optionally drop patterns that do not meet the mining criteria.}
  \item{\dots}{further (unused) arguments.}
}
\details{
  If option \code{drop} is \code{TRUE} only patterns that satisfy the 
  minimum size and maximum radius constraint are extracted. 
}
\value{
  A list with the following components:

  \item{x}{the fitted data matrix.}
  \item{pl}{a factor of pattern (cluster) labels. The indexes of the 
    original data samples are provided as attribute \code{Index}.}
}
\author{Christian Buchta}
\seealso{
  \code{\link{proximus}} for pattern mining with the Proximus algorithm.
}
\examples{
### see proximus
}
\keyword{cluster}

