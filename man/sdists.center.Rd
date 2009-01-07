\name{sdists.center}
\alias{sdists.center}
\title{Centroid Sequences}
\description{
  Find centroid sequences among a collection of sequences.
}
\usage{
sdists.center(x, d = NULL, method = "ow", weight = c(1, 0, 2),
              exclude = c(NA, NaN, Inf, -Inf), FUN = NULL, ...,
	      unique = FALSE)
}
\arguments{
  \item{x}{a list (of vectors) of a vector of character.}
  \item{d}{a matrix or an object of class \code{dist}.}
  \item{method}{argument to \code{sdists}.}
  \item{weight}{argument to \code{sdists}.}
  \item{exclude}{argument to \code{sdists}.}
  \item{FUN}{a function to rank distances.}
  \item{\dots}{additional arguments to \code{FUN}.}
  \item{unique}{a logical specifying whether to return a unique
    set of sequences.}
}
\details{
  This function provides a wrapper to computing the distances among
  the sequences in \code{x}, unless \code{d} is supplied, and the
  subsequent selection of a set of centroid sequences with minimum
  sum of distances to any other sequence.
}
\value{
  A subset of \code{x}.
}
\author{Christian Buchta}
\seealso{\code{\link{sdists}} for distance computation.}
\examples{
x <- c("ABCD", "AD", "BCD", "ACF", "CDF", "BC")
sdists.center(x)
}
\keyword{cluster}
