\name{sdists.center.align}
\alias{sdists.center.align}
\title{Align Sequences to a Center}
\description{
  Find a global alignment of a collection of sequences using the
  center-star-tree heuristic. 
}
\usage{
sdists.center.align(x, center, method = "ow", weight = c(1, 0, 2),
                    exclude = c(NA, NaN, Inf, -Inf), break.ties = TRUE,
		    transitive = FALSE, to.data.frame = FALSE)
}
\arguments{
  \item{x}{a list (of vectors) or a vector of character.}
  \item{center}{a vector}
  \item{method}{argument to \code{sdists}.}
  \item{weight}{argument to \code{sdists}.}
  \item{exclude}{arguments to \code{sdists}.}
  \item{break.ties}{a logical specifying whether random tie-breaking
    should be perforemd. Otherwise the first alignment is used.}
  \item{transitive}{a logical specifiying whether the sequences in
    \code{x} should be aligned with each other, too.}
  \item{to.data.frame}{a logical specifying whether the result should
    be converted to \code{data.frame}.}
}
\details{
  Each component of \code{x} is aligned with \code{center} in turn
  such that the latter is aligned with all sequences processed so
  far.

  If \code{center} is missing \code{i\link{sdists.center}} is used
  to compute an initial center.
}
\value{
  Either a list of sequences with attributes \code{center} and
  \code{ties}, or a \code{data.frame} with the sequences in
  the columns.
}
\references{
  D. Gusfield (1997). \emph{Algorithms on Strings, Trees, and
  Sequences}.  Cambridge University Press, Chapter XX.
}
\author{Christian Buchta}
\note{ 
  The global alignment may depend on the order of \code{x}. 
}
\seealso{
\code{\link{sdists}} for computation of distances, 
\code{\link{sdists.center}} for computation of centroids.
}
\examples{
## continue example
x <- c("ABCD", "AD", "BCD", "ACF", "CDF", "BC")
sdists.center.align(x)
sdists.center.align(x, transitive = TRUE, to.data.frame = TRUE)
}
\keyword{cluster}
