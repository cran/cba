\name{dists}
\alias{dists}
\alias{dapply}
\title{Matrix Distance Computation}
\description{
  This function computes and return the auto-distance matrix between the 
  rows of a matrix, as well as the cross-distance matrix between two matrices.
}
\usage{
dists(x, y = NULL, method = "minkowski", p = 2)

dapply(x, y, FUN, ...)
}
\arguments{
  \item{x}{a numeric matrix object.}
  \item{y}{\code{NULL}, or a numeric matrix object.}
  \item{method}{a mnemonic string referencing the distance measure.}
  \item{p}{Minkowski metric parameter.}
  \item{FUN}{a user supplied function.}
  \item{\dots}{further arguments to the user supplied function.}
}
\details{
  The interface is fashioned after \code{dist}: you have to specify a 
  \code{method} to use, i.e. a (not so) mnemonic name. 

  Methods that are also implemented in \code{\link{dist}} are: 
  \code{minkowski}, \code{maximum}, \code{canberra}, and \code{binary}. 
  See the documentation there. Note that for \code{binary} the 
  arguments \code{x} (and \code{y}) must be \code{logical}.

  Additional methods implemented are:
  
\describe{
  \item{\code{ebinary}:}{one minus extended Jaccard similarities on 
    real-valued vectors using Euclidean distances.}

  \item{\code{fbinary}:}{one minus fuzzy Jaccard similarities on 
    positive real-valued vectors as proposed by Kurt Hornik.}

  \item{\code{angular}:}{one minus cosine similarities on real-valued 
    vectors.}
}

  Missing values are allowed but are excluded from all computations 
  involving the rows within which they occur. However, rows (and columns)
  of NAs are not dropped as in \code{dist}.

  For compatibility the distance is zero instead of \code{NA} in the case
  two (near) zero vectors are involved in the computation of \code{binary}, 
  \code{ebinary}, and \code{angular}. Note that this is inconsistent with
  the coding of \code{NA} by \code{as.dummy}.

  Function \code{dapply} allows the user to apply arbitrary distance
  functions that take as arguments at least two vectors (i.e. rows of
  \code{x}, etc.) and return a scalar real value.
}
\value{
  Auto distances are returned as an object of class \code{dist} and 
  cross-distances as an object of class \code{matrix}. 
}
\author{Christian Buchta}
\section{Warning}{The interface is experimental and may change in the future.}
}
\seealso{\code{\link{dist}} for compatibility information.}
\examples{
### binary data
x <- matrix(sample(c(FALSE,TRUE),8,rep=TRUE), ncol=2)
dists(x, method="binary")
### for real-valued data
dists(x, method="ebinary")
### for positive real-valued data
dists(x, method="fbinary")
### cross distances
dists(x, x, method="binary")
### this is the same but less efficient
as.matrix(dists(x, method="binary"))
## test inheritance of names
rownames(x) <- LETTERS[1:4]
dists(x)
dists(x,x)
## custom distance function
f <- function(x, y) sum(x*y)
dapply(x, FUN=f)
dapply(x,x, FUN=f)
}
\keyword{cluster}

