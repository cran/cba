\name{plotting}
\alias{image.matrix}
\alias{image.dist}
\title{Matrix Image Plots}
\description{
  Provides methods for proper plotting objects of class \code{matrix} and
  \code{dist}. 
}
\usage{
\method{image}{matrix}(x, xlab = "", ylab = "", ...)

\method{image}{dist}(x, xlab = "", ylab = "", ...)
}
\arguments{
  \item{x}{an object of class \code{matrix} or \code{dist}.}
  \item{xlab, ylab}{labels for the x and y axis.}
  \item{\dots}{further arguments to \code{image}.}
}
\details{
  Plots an object of class \code{matrix} in its original row and column
  orientation. This means, in a plot the columns become the x-coordinates 
  and the reversed rows the y-coordinates.

  If \code{x} is of class \code{dist} it is converted to full-storage
  representation and dispatched to the matrix method.
}
\author{Christian Buchta}
\seealso{
  \code{\link{image}} for details.}
\examples{
x <- matrix(sample(c(FALSE, TRUE),100,rep=TRUE),ncol=10)
image(x)
image(x, col=c("white","black"))
axis(1)				# the columns
axis(2, at=c(1,3,5,7,9), labels=c(9,7,3,5,1))	# the reversed rows
}
\keyword{cluster}
\keyword{hplot}
