\name{lmplot}
\alias{lmplot}
\title{Plotting Logical Matrices}
\description{
  Implements a wrapper function to image that produces a black and white or 
  gray-scale plot of a logical matrix.
}
\usage{
lmplot(x, block.size = 1, gray = FALSE, xlab = "", ylab = "", ...)
}
\arguments{
  \item{x}{a logical matrix}
  \item{block.size}{the interpolation block size}
  \item{gray}{optionally use a gray scale}
  \item{xlab}{title for the x axis}
  \item{ylab}{title for the y axis}
  \item{\dots}{further arguments to \code{image}}
}
\details{
  \code{TRUE} is represented by the color white and \code{FALSE} by the 
  color black. 
  
  A lower resolution can be obtained by specifying an (interpolation) block 
  size greater than one. Block densities can then be visualized by using 
  the gray scale option. The number of levels of the palette corresponds to 
  the block size but is capped to 8 levels (excluding white). Note that 
  the opacity (blackness) corresponds with density (as on photographic 
  film).
}
\author{Christian Buchta}
%\section{Warning}{}
\seealso{ 
  \code{\link{lminter}} for interpolating logical matrices and 
  \code{\link{image}} for further plotting options}
\examples{
###
x <- matrix(sample(c(FALSE, TRUE), 64, rep=TRUE), ncol=8)
lmplot(x)
### use lower resolution
lmplot(x, block.size=2)
### use gray scale
lmplot(x, block.size=2, gray=TRUE)
}
\keyword{cluster}
\keyword{hplot}

