\name{rlbmat}
\alias{rlbmat}
\title{Block Uniform Logical Matrix Deviates}
\description{
  Generate a uniform logical matrix deviate with a possibly overlapping 
  block structure. 
}
\usage{
rlbmat(npat = 4, rows = 20, cols = 12, over = 4, noise = 0.01, prob = 0.8, perfect = FALSE)
}
\arguments{
  \item{npat}{number of patterns}
  \item{rows}{number of rows per pattern}
  \item{cols}{number of columns per pattern}
  \item{over}{number of additional columns per pattern that overlap}
  \item{noise}{the probability of observing a one in the background (non-pattern) matrix}
  \item{prob}{the probability of observing \code{TRUE} in a pattern}
  \item{perfect}{option for overlap of the first and the last pattern}
}
\details{
  Implements a test case for \code{proximus}. 
}
\value{
  A logical matrix
}
\author{Christian Buchta}
%\note{}
\seealso{
  \code{\link{lmplot}} and 
  \code{\link{clmplot}} for plotting a logical matrix
}
\examples{
x <- rlbmat()
lmplot(x)
}
\keyword{cluster}
