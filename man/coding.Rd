\name{coding}
\alias{as.dummy}
\alias{as.dummy.matrix}
\alias{as.dummy.data.frame}
\title{Dummy Coding}
\description{
  Functions that convert R objects to a dummy coded matrix (or a list of 
  matrices). 
}
\usage{
as.dummy(x, ...)

\method{as.dummy}{matrix}(x, sep = " ", drop = FALSE, ...)
\method{as.dummy}{data.frame}(x, sep = " ", drop = FALSE, ...)
}
\arguments{
  \item{x}{an R object (see below).}
  \item{sep}{separator used for construction of colnames.}
  \item{drop}{drop factors with less than two levels.}
  \item{\dots}{other (unused) arguments.}
}
\details{
  The generic is applicable to factor and to R objects that can be 
  converted to factor, i.e. \code{logical}, \code{integer}, or 
  \code{character}. For numeric data a discretization method has to be 
  used.

  A factor is converted to as many logical variables as there are levels 
  where the value \code{TRUE} indicates the presence of a level.

  The colnames are made of the concatenation of a variable name and the 
  level, separated by \code{sep}. For \code{matrix} and \code{data.frame} 
  variable names are created if necessary.

  A value of \code{NA} is mapped to \code{FALSE} across all levels.
}
\value{
  A \code{matrix} with a \code{levels} attribute which contains a list of 
  the levels of the coded variables.
}
\author{Christian Buchta}
\section{Warning}{This is experimental code which may change in the future.}
\seealso{\code{\link{as.logical}}.}
\examples{
### 
x <- as.integer(sample(3,10,rep=TRUE))
as.dummy(x)
is.na(x) <- c(3,5)
as.dummy(x)
x <- as.data.frame(x)
as.dummy(x)
}
\keyword{cluster}
