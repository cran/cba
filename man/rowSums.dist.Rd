\name{rowSums.dist}
\alias{rowSums.dist}
\title{Row Sums of Sparse Symmetric Matrices}
\description{
  Compute the row (column) sums for a sparse symmetric (distance) matrix.
}
\usage{
rowSums.dist(x, na.rm = FALSE)
}
\arguments{
  \item{x}{an object of class \code{dist}.}
  \item{na.rm}{logical, should missing values (including \code{NaN}) be
    omitted from the summation?}
}
\details{
  This function is more efficient than expanding an object of
  class \code{dist} to matrix and using \code{rowSums}.
}
\value{
  A numeric vector of row sums.
}
\author{Christian Buchta}
\seealso{\code{as.matrix}, \code{as.dist}, and \code{rowSums}.}
\examples{
## 
x <- matrix(runif(10*2),ncol=2)
d <- dist(x)
rowSums(as.matrix(d))
rowSums.dist(d)
}
\keyword{cluster}
