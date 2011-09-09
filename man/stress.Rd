\name{stress}
\alias{stress}
\alias{stress.dist}
\title{Conciseness of Presentation Measures}
\description{
  Compute different types of conciseness of presentation measures.
}
\usage{
stress(x, rows = NULL, cols = NULL, type = "moore")

stress.dist(x, rows = NULL, cols = NULL, bycol = FALSE,
            type = "moore")
}
\arguments{
  \item{x}{a matrix object.}
  \item{rows}{a subscript vector indexing the rows.}
  \item{cols}{a subscript vector indexing the columns.}
  \item{bycol}{logical for computation over the columns.}
  \item{type}{the type of neighborhood to use.}
}
\details{
  Function \code{stress} computes the sum of squared distances of each
  matrix entry from its adjacent entries. The following types of neighborhoods
  are available:
  
  \describe{
    \item{\code{moore}:}{comprises the eight adjacent entries (five at the
	margins and three at the corners).}
    \item{\code{neumann}:}{comprises the four adjacent entries (three at the
	margins and two at the corners).}
  }

  Function \code{stress.dist} computes the auto-distance matrix for each pair
  of rows (or columns) given one of the above stress measures. Note that the
  result depends on the ordering of the unused dimension.

  As the computation can be reduced to summing the edge distances between any
  two neighboring points, only half of the value of the proposed measures is
  reported.

  Row and/or column indexes (or labels) can be supplied to test specific
  orderings, as well as subsets of indexes (labels).

  Note that the matrix should be normalized so that the distance computation
  is meaningful.
}
\value{
  \code{stress} returns a scalar real, i.e. half of the global stress measure.

  \code{stress.dist} returns an object of class \code{\link{dist}}, i.e. a
  lower triangular matrix in column format.
}
%\references{}
\author{Christian Buchta}
%\note{}
\seealso{
  \code{\link[proxy:dist]{dist}} for general distance computation.
}
\examples{
##
x1 <- matrix(sample(c(FALSE,TRUE),25,rep=TRUE),ncol=5)
stress(x1)
stress(x1, type="neumann")
##
x2 <- cbind(rbind(matrix(1,4,4),matrix(0,4,4)), 
	    rbind(matrix(0,4,4),matrix(1,4,4)))
stress.dist(x2)
stress.dist(x2, bycol=TRUE)
stress.dist(x2, type="neumann")
}
\keyword{hplot}
\keyword{cluster}
