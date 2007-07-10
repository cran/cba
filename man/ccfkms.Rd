\name{ccfkms}
\alias{ccfkms}
\title{Clustering with Conjugate Convex Functions}
\description{
  Partition a data set into convex sets using conjugate convex functions.
}
\usage{
ccfkms(x, n, p = NULL, par = 2, max.iter = 100, opt.std = FALSE,
       opt.retry = 0, debug = FALSE)
}
\arguments{
  \item{x}{a data matrix.}
  \item{n}{optional number of prototypes.}
  \item{p}{a matrix of initial prototypes.}
  \item{par}{type or parameter of conjugate convex function.}
  \item{max.iter}{maximum number of iterations.}
  \item{opt.std}{optionally standardize the data.}
  \item{opt.retry}{number of retries.}
  \item{debug}{optionally turn on debugging output.}
}
\details{
  Two types of conjugate convex functions are available: one that is based on 
  powers of the norm of the prototype vectors and another that is based on a 
  logarithmic transformation of the norm. Both are intended to obtain more 
  robust partitions.

  Using \code{par} = 2 is equivalent to performing ordinary k-means with
  Euclidean distances. \code{par} = 1 is equivalent to LVQ of Kohonen type 
  (the directions of the prototypes from the center of the data are used), 
  and \code{par} = 0 is equivalent to using 2*ln(cosh(|p|))/2.
  
  Internally the algorithm uses sparse data structures and avoids computations 
  with zero data values. Thus, the data must not be centered (the algorithm 
  does this internally with the option to further standardize the data). For 
  dense data this is slightly inefficient.
  
  If initial prototypes are omitted the number of prototypes must be specified.
  In this case the initial prototypes are drawn from the data (without 
  replacement).

  If the number of retries is greater than zero the best among that number 
  of trial solutions is returned. Note that the number of prototypes must be
  specified as the initial prototypes are sampled from the data.

  The debugging output shows the iteration number, the inverted information 
  and the variance of the current partition as a percentage of the total (if 
  each data point were a cluster), and the number of active prototypes (those 
  with at least one member, i.e. a data point that is not closer to any 
  other prototype).

  Note that the algorithm uses tie-breaking when it determines the cluster
  memberships of the samples.
}
\value{
  A list with the following components:

  \item{centers}{a matrix of cluster means (final prototypes).}
  \item{size}{a vector of cluster sizes.}
  \item{cl}{a factor of cluster labels (indexes).}
  \item{inv.inf}{the inverted information of the partition.}
  \item{par}{see above.}
  \item{opt.std}{see above.}
}
\references{
  Helmut Strasser and Klaus Poetzelberger. Data Compression by 
  Unsupervised Classification. \emph{SFB Report Series}, No. 10, 1997.
}
\author{Christian Buchta}
\note{
Support for data matrices \code{x} in sparse \code{dgTMatrix} and 
\code{dgCMatrix} format (see package \pkg{Matrix}) is experimental. 
Support for the \code{dgRMatrix} format is currently suspended
due to problems with package \pkg{Matrix}.
}
\seealso{
  \code{\link{kmeans}}, \code{cmeans}, \code{kkmeans} for similar or related
  clustering techniques.
}
\examples{
### extend proximus example
x <- rlbmat()
rownames(x) <- seq(dim(x)[1])
cm <- ccfkms(x, n=4, opt.retry=10)
pcm <- predict(cm, x)
\dontrun{
### using sparse data may be more time-efficient
### depending on the goodness of the implementation
### of subset, etc. in package Matrix.
require(Matrix)
#sx <- as(x, "dgRMatrix")    # currently broken
sx <- as(x, "dgCMatrix")
system.time(scm <- ccfkms(sx, n=4, opt.retry=50))
system.time(cm <- ccfkms(x, n=4, opt.retry=50))
}
}
\keyword{cluster}
