\name{cluster.dist}
\alias{cluster.dist}
\title{Clustering a Sparse Symmetric Distance Matrix}
\description{
  Compute a clustering on a sparse symmetric distance matrix using
  graph cutting.
}
\usage{
cluster.dist(x, beta)
}
\arguments{
  \item{x}{an object of class \code{dist}.}
  \item{beta}{the distance threshold.}
}
\details{
  This function computes a clustering on an object of class \code{dist}
  by cutting the graph induced by the threshold \code{beta} into all
  disconnected subgraphs (the clusters). Two nodes are connected by a
  link if their distance is below the specified threshold. Note that
  the threshold is not strict, i.e. \code{>=}.

  Note that distances of value \code{NA} and \code{NaN} are ignored.
  This is not strictly correct but avoids computing \eqn{2^k} possible
  solutions if there are k \code{NA} values.

  The time complexity is \eqn{O(n^2)} with n the number of rows/columns.
}
\value{
  A factor of cluster labels (indexed 1,2,\dots,k).
}
%\references{
% fixme
%}
\author{Christian Buchta}
\note{
  Fixme: can the time complexity be improved?
}
\seealso{
  \code{\link{dist}} and \code{\link{sdists}} for distance computation.}
\examples{
## 3 clusters (1 = connected)
x <- matrix(c(1,1,0,0,0,0,
	      1,1,0,0,0,0,
	      0,0,1,1,0,0,
	      0,0,1,1,0,0,
	      0,0,0,0,1,1,
	      0,0,0,0,1,1), ncol=6)
c <- cluster.dist(as.dist(!x), beta = 0) # invert and note that 0 >= 0
c
}
\keyword{cluster}
