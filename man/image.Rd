\name{image}
\alias{implot}
\title{Matrix Image Plots}
\description{
  Implements a wrapper function to \code{image} for proper plotting of
  objects of class \code{matrix} and \code{dist}.
}
\usage{
implot(x, xlab = "", ylab = "", axes = FALSE, ticks = 10,
       las = 2, ...)
}
\arguments{
  \item{x}{an object of class \code{matrix} or \code{dist}.}
  \item{xlab, ylab}{labels for the x and y axis.}
  \item{axes}{logical, indicating whether \code{dimnames(x)} should be
	      drawn on the plot.}
  \item{ticks}{the number of tick-marks to use.}
  \item{las}{the axis style to use (see \code{par}).}
  \item{\dots}{further arguments to \code{image}.}
}
\details{
  Plots an object of class \code{matrix} in its original row and column
  orientation. This means, in a plot the columns become the x-coordinates 
  and the reversed rows the y-coordinates.

  If \code{x} is of class \code{dist} it is coerced to full-storage
  \code{matrix} representation.
}
\value{
  Returns the transformed \code{x} \emph{invisibly}.
}
\author{Christian Buchta}
\seealso{
  \code{\link{image}} and \code{\link{par}} for details.}
\examples{
x <- matrix(sample(c(FALSE, TRUE),100,rep=TRUE),ncol=10,
	    dimnames=list(1:10,LETTERS[1:10]))
implot(x)
implot(x, col=c("white","black"), axes = TRUE)
}
\keyword{cluster}
\keyword{hplot}
