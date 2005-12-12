\name{subset.dist}
\alias{subset.dist}
\alias{[[.dist}
\title{Subsetting Compact Symmetric Matrices}
\description{
  This function conveniently returns a subset (submatrix) of an object of class \code{dist}.
}
\usage{
\method{subset}{dist}(x, subset, ...)
\method{[[}{dist}(x, subset, ...)

x[[i]]
}
\arguments{
  \item{x}{an object of class \code{dist}}
  \item{subset, i}{logical, character, or index vector}
  \item{\dots}{other (unused) arguments}
}
\details{
  For consistency, subsetting a symmetric matrix is constrained to a 
  single index (vector). However,  non-unique indexing is allowed. In
  this case zero is returned for diagonal elements.
  
  If \code{x} has a \code{Labels} attribute a character vector can
  be used for subsetting.
  
  Note that objects which contain diagonal elements are not handled.
}
\value{
  An object of class \code{dist}.
}
\author{Christian Buchta}
\section{Warning}{This is experimental code which may change in future releases.}
\seealso{\code{\link{dist}} for details}
\examples{
x <- as.dist(matrix(1:36, ncol=6))
subset(x, c(2,4,6))
attr(x, "Labels") <- letters[1:6] # hack!
subset(x, c(2,4,6))
subset(x, letters[c(2,4,6)])
subset(x, c(2,2,4,6))   # contains one zero 
# this is more user friendly
x[[c(2,4,6)]]
x[[c("b","d","e")]]
}
\keyword{cluster}
