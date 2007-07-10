\name{subset.dist}
\alias{subset.dist}
\alias{[[.dist}
\alias{row.dist}
\alias{col.dist}
\title{Subsetting Sparse Symmetric Matrices}
\description{
  This function conveniently returns a subset (submatrix) of a sparse
  symmetric (distance) matrix.
}
\usage{
\method{subset}{dist}(x, subset, ...)
\method{[[}{dist}(x, subset, ...)

x[[i]]

row.dist(x)
col.dist(x)
}
\arguments{
  \item{x}{an object of class \code{dist}.}
  \item{subset, i}{logical, character, or index vector.}
  \item{\dots}{further (unused) arguments.}
}
\details{
  For consistency, subsetting a symmetric matrix is constrained to a 
  single index (vector). However,  non-unique indexing is allowed. In
  this case zero is returned for diagonal elements.
  
  If \code{x} has a \code{Labels} attribute a character vector can
  be used for subsetting.

  Note that an error is raised if \code{subset} contains \code{NA}, 
  or a non-finite value.

  \code{row.dist} and \code{col.dist} return the matrix indexes of
  the row or columns, respectively.
}
\value{
  An object of class \code{dist}.
}
\author{Christian Buchta}
\section{Warning}{This is experimental code which may change in future releases.}
\seealso{\code{\link{dist}} for details.}
\examples{
x <- as.dist(matrix(1:36, ncol=6))
subset(x, c(2,4,6))
subset(x, c(4,2,6))               # same set but different order
attr(x, "Labels") <- letters[1:6] # hack!
subset(x, c(2,4,6))
subset(x, letters[c(2,4,6)])
subset(x, c(2,2,4,6))		  # contains one zero 
# this is more user friendly
x[[c(2,4,6)]]
x[[c("b","d","e")]]
## finally, don't use NA, etc.
z <- try(x[[c(2,NA,6)]], silent=TRUE)
z
# this is arguable
x[[1]]				  # vector of length zero

# indexes
row.dist(x)
col.dist(x)

# convenience
dim(x)
dimnames(x) <- LETTERS[1:6]
dimnames(x)
}
\keyword{cluster}
