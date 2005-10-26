\name{order.greedy}
\alias{order.greedy}
\title{Hierarchical Greedy Ordering}
\description{
  Compute a hierarchical greedy ordering of a data matrix.
}
\usage{
order.greedy(dist)
}
\arguments{
  \item{dist}{an object of class \code{dist}.}
}
\details{
  A single cluster is constructed by merging in each step the leaf 
  closest to one of the two endpoints of the cluster. The algorithm
  starts with a random leaf and uses tie-breaking.

  Clearly, the algorithm is more an ordering than a cluster algorithm. 
  However, it constructs a binary merge tree so that the linear ordring 
  of its leaves could be further improved.
}
\value{
  A list with the following items:
  
  \item{merge}{a matrix containing the merge tree.}
  \item{order}{a vector containing the leaf ordering.}
  \item{height}{a vector containing the merge heights.}
}
\references{
  F. Murtagh (1985). Multidimensional Cluster Algorithms. Lectures in
  Computational Statistics, Physica Verlag, pp. 15.
}
\author{Christian Buchta}
\note{The merge heights may not be monotonic.}
\seealso{
  \code{\link{hclust}} for hierarchical clustering,
  \code{\link{order.optimal}} for optimal leaf ordering, and
  \code{\link{order.length}} for computing the objective value of a 
    leaf ordering.
}
\examples{
d <- dist(matrix(runif(20), ncol=2))
hc <- hclust(d)
co <- order.optimal(d, hc$merge)
md <- -as.dist(crossprod(as.matrix(d)))	# Murtagh's distances
hg <- order.greedy(md)
go <- order.optimal(md, hg$merge)
### compare images
op <- par(mfrow=c(2,2), pty="s")
image(d[[hc$order]], main="hclust")
image(d[[co$order]], main="hlcust + optimal")
image(d[[hg$order]], main="greedy")
image(d[[go$order]], main="greedy + optimal")
par(op)
# compare lenghts
order.length(d, hc$order)
order.length(d, co$order)
order.length(d, hg$order)
order.length(d, go$order)
}
\keyword{hplot}
\keyword{cluster}
