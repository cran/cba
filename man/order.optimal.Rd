\name{order.optimal}
\alias{order.optimal}
\title{Optimal Leaf Ordering of Binary Trees.}
\description{
  Find an optimal linear leaf ordering of a binary merge tree as
  produced by a hierarchical cluster algorithm.
}
\usage{
order.optimal(dist, merge)
}
\arguments{
  \item{dist}{an object of class \code{dist}.}
  \item{merge}{a binary merge tree (see \code{\link{hclust}}).}
}
\details{
  A binary tree has 2\^(n-1) internal nodes (subtrees) and the same
  number of leaf orderings. That is, at each internal node the left 
  and right subtree (or leaves) can be swapped, or, in terms of a
  dendrogram, be flipped.

  An objective measure of a leaf ordering is the sum of the distances
  along the path connecting the leaves in the given order. An ordering 
  with a minimal path length is defined to be an optimal ordering.
  
  This function provides an interface to the optimal leaf ordering
  algorithm (see references) for tree representations that are used by 
  hierarchical cluster algorithms such as \code{\link{hclust}}.

  Note that non-finite distance values are not allowed.
}
\value{
  A list with the following components:

  \item{merge}{a matrix containing the merge tree corresponding with the 
    optimal leaf order.}
  \item{order}{a vector containing the optimal leaf order.}
  \item{length}{the length of the ordering.}
}
\references{
  Z. Bar-Joseph, E. D. Demaine, D. K. Gifford, and T. Jaakkola.
  (2001) Fast Optimal Leaf Ordering for Hierarchical Clustering.
  Bioinformatics, Vol. 17 Suppl. 1, pp. 22-29. 
}
\author{Christian Buchta}
\note{The time complexity of the algorithm is o(n\^3).}
\seealso{
  \code{\link{hclust}} for hierarchical clustering and
  \code{\link{order.length}} for computing the objective value of a 
    leaf ordering.
}
\examples{
d <- dist(matrix(runif(30), ncol=2))
hc <- hclust(d)
co <- order.optimal(d, hc$merge)
### compare dendrograms
ho <- hc
ho$merge <- co$merge
ho$order <- co$order
op <- par(mfrow=c(2,2), pty="s")
plot(hc, main="hclust")
plot(ho, main="optimal")
# compare images
image(d[[hc$order]])
image(d[[co$order]])
par(op)
### compare lengths
order.length(d, hc$order)
order.length(d, co$order)
cat("compare: ",co$length,"\n")
}
\keyword{hplot}
\keyword{cluster}
