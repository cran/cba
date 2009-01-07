\name{gknn}
\alias{gknn}
\title{Generalized k-Nearest Neighbor Classification}
\description{
  Compute the k-nearest neighbor classification given a matrix of 
  cross-distances and a factor of class values. For each row the majority 
  class is found, where ties are broken at random (default). If there are ties 
  for the kth nearest neighbor, all candidates are included in the vote 
  (default).
}
\usage{
gknn(x, y, k = 1, l = 0, break.ties = TRUE, use.all = TRUE,
     prob = FALSE)
}
\arguments{
  \item{x}{a cross-distances matrix.}
  \item{y}{a factor of class values of the columns of \code{x}.}
  \item{k}{number of nearest neighbors to consider.}
  \item{l}{minimum number of votes for a definite decision.}
  \item{break.ties}{option to break ties.}
  \item{use.all}{option to consider all neighbors that are tied with the 
    kth neighbor.}
  \item{prob}{optionally return proportions of winning votes.}
}
\details{
  The rows of the cross-distances matrix are interpreted as referencing the 
  test samples and the columns as referencing the training samples.

  The options are fashioned after \code{knn} in package \pkg{class} but are
  extended for tie breaking of votes, e.g. if only definite (majority) votes 
  are of interest.

  Missing class values are not allowed because that would collide with a 
  missing classification result.

  Missing distance values are ignored but with the possible consequence 
  of missing classification results. Note that this depends on the 
  options settings, e.g. 
}
\value{
  Returns a factor of class values (of the rows of \code{x}) which may be 
  \code{NA} in the case of doubt (no definite decision), ties, or missing
  neighborhood information.

  The proportions of winning votes are returned as attribute \code{prob}
  (if option \code{prob} was used).
}
\author{Christian Buchta}
\seealso{\code{\link[proxy:dist]{dist}} for efficient computation of cross-distances.}
\examples{
\dontrun{
### extend Rock example
data(Votes)
x <- as.dummy(Votes[-17])
rc <- rockAll(x, n=2, m=100, theta=0.73, predict=FALSE, debug=TRUE)
gc <- gknn(dist(x, rc$y, method="binary"), rc$cl, k=3)
table(gc[rc$s], rc$cl)
}
}
\keyword{cluster}
\keyword{classif}
