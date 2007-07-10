\name{order}
\alias{order.dist}
\alias{order.matrix}
\alias{order.data.frame}
\title{Improving the Presentation of Matrix Objects}
\description{
  High-level functions that improve the presentation of a matrix or data
  frame by reordering their rows and columns.
}
\usage{
order.dist(x, index = FALSE)
order.matrix(x, type = "neumann", by = c("both","rows","cols"),
	     index = FALSE)
order.data.frame(x, type = "neumann", by = c("both","rows","cols"),
                 index = FALSE)
}
\arguments{
  \item{x}{an object of class \code{dist}, \code{matrix}, or
    \code{data.frame}.}
  \item{type}{the type of stress measure to use (see details).}
  \item{by}{option to order either by rows, or columns, or both.}
  \item{index}{option to return the order index(ex) instead of the
    reordered object.}
}
\details{
  These functions try to improve the presentation of an object of class
  \code{dist}, \code{matrix}, or\cr\code{data.frame} by reordering the rows
  and columns such that similar entries are grouped together.

  \code{order.dist} uses a simple heuristic to solve the TSP problem of
  finding an ordering of minimum length (see \code{order.length}) for an
  object of class \code{dist}. Note that the heuristic used is quick but
  more elaborate TSP algorithms will produce better orderings.

  \code{order.matrix} tries to minimize the stress measure of a matrix
  (see \code{stress} by using the same TSP heuristic as above, once for
  the column and once for the row ordering (while the other dimension is
  fixed) if \code{by = "both"}.

  \code{order.data.frame} uses attributes of type \code{numeric}
  and \code{logical} only, combines them into a normalized matrix and
  finds an ordering as above.
}
\value{
  Either the reordered object supplied, or a vector of subscripts (for
  \code{reorder.dist}), or a list with components \code{rows} and
  \code{columns} containing the order indexes (for \code{reorder.matrix} 
  and \code{reorder.data.frame}).
}
%\references{}
\author{Christian Buchta}
\note{
  This is experimental code that may be integrated in a separate packages
  in the future.
}
\seealso{\code{\link{dist}}, \code{stress}, \code{stress.dist}.}
\examples{
## not a hard problem
data(iris)
d <- dist(iris[1:4])
image(order.dist(d))
data(townships)
x <- order.data.frame(townships)
x
}
\keyword{cluster}
