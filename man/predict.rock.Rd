\name{predict.rock}
\alias{predict.rock}
\alias{fitted.rock}
\title{Rock Clustering}
\description{
  Classify the rows of a data matrix using the Rock classifier.
}
\usage{
\method{predict}{rock}(object, x, drop = 1, ...)

\method{fitted}{rock}(object, ...)
}
\arguments{
  \item{object}{an object of class \code{rock}}
  \item{x}{a data matrix containing test or new samples}
  \item{drop}{cluster size threshold}
  \item{\dots}{further (unused) arguments}
}
\details{
  Provides a wrapper function to the Rock classifier: cluster memberships 
  of the (row) samples are determined by majority vote using the size (of 
  the cluster) weighted number of links. 
  
  The classifier uses random tie-breaking and assigns samples with zero link 
  counts to \code{NA}. This allows for detection of possible outliers, or
  interestin patterns.

  In the case a \code{drop} value greater than zero is specified, all 
  clusters with size equal or less than this value are removed from the 
  classifier. Especially, \code{fitted} uses a threshold of one because
  for singleton clusters the neighborhood is empty.

  Note that for the training data the predicted memberships need not 
  necessarily be the same as those obtained from the cluster algorithm.
}
\value{
  A list with the following elements:

  \item{cl}{a factor of cluster memberships of the samples.}
  \item{size}{a vector of cluster sizes.}
}
\seealso{
  \code{\link{rockCluster}} for obtaining a rock classifier and
  \code{\link{gknn}} for k-nearest neighbor classification.
}
\examples{
### example from paper
data(Votes)
x <- as.dummy(Votes[-17])
rc <- rockCluster(x, n=2, theta=0.73, debug=TRUE)
rp <- predict(rc, x)
table(rp$cl)
}
\keyword{cluster}
%%
