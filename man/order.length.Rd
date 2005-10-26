\name{order.length}
\alias{order.length}
\title{Conciseness of Presentation Measures}
\description{
  Compute the length of a Hamilton path through a distance matrix.
}
\usage{
order.length(dist, order)
}
\arguments{
  \item{dist}{an object of class \code{dist}.}
  \item{order}{an optional permutation of the row (column) indexes.}
}
\details{
  Ordering a distance matrix such that low distance values are placed 
  close to the diagonal may improve its presentation. The length of an
  order is the corresponding objective measure.

  The order corresponds to a path through a graph where each node is
  visited only once, i.e. a Hamilton path. The length of a path is defined 
  as the sum of the edge weights, i.e. distances. 

  If \code{order} is missing the identity order is used.
  
  If \code{order} is not unique \code{NA} is returned.

  If there are non-finite distance values \code{NA} is returned.
}
\value{
  A single real value.
}
\references{
  R. Sedgewick. (2002) Algorithms in C. Part 5. Graph Algorithms.
  3rd Edition, Addison-Wesley.
}
\author{Christian Buchta}
\examples{
d <- dist(matrix(runif(10),ncol=2))
order.length(d)
o <- sample(5,5) # random order
order.length(d, o)
}
\keyword{hplot}
\keyword{cluster}
