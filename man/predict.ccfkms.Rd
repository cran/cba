\name{predict.ccfkms}
\alias{predict.ccfkms}
\title{Clustering with Conjugate Convex Functions.}
\description{
  Classify the rows of a data matrix using conjugate convex functions.
}
\usage{
\method{predict}{ccfkms}(object, x, drop = 1, ...)
}
\arguments{
  \item{object}{an object of class \code{ccfkms}}
  \item{x}{a data matrix containing test or new samples}
  \item{drop}{cluster size threshold}
  \item{\dots}{other (unused) arguments}
}
\details{
  This is a wrapper to \code{\link{ccfkms}} which uses a single iteration 
  for classifying the data.

  In the case a \code{drop} value greater than zero is specified, all
  clusters with size equal or less than this value are removed from the
  classifier.
}
\value{
  A list with the following elements:

  \item{centers}{a matrix of cluster means.}
  \item{size}{a vector of cluster sizes.}
  \item{cl}{a factor of cluster labels (indexes).}
  \item{inv.inf}{the inverted information of the partition.}
}
\author{Christian Buchta}
%\note{}
\seealso{
  \code{\link{ccfkms}} for obtaining a classifier.}
\examples{
### see ccfkms
}
\keyword{cluster}
