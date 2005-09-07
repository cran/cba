\name{cut.ordered}
\alias{cut.ordered}
\title{Converting Ordered Factors}
\description{
  Reduce the levels of an ordered factor.
}
\usage{
\method{cut}{ordered}(x, breaks, ...)
}
\arguments{
  \item{x}{an ordered factor.}
  \item{breaks}{a logical, character, or index vector of cut points.}
  \item{\dots}{further (unused) arguments}
}
\details{
  If \code{breaks} is of class \code{logical} it must have the same length 
  as the number of levels of \code{x}.

  If \code{breaks} is of class \code{character} partial matching with
  the levels of \code{x} is attempted.

  Otherwise \code{breaks} is assumed to index the levels.
}
\value{
  An ordered factor.
}
\references{Functions Missing in R: A Never Ending Story ;-)}
\author{Christian Buchta}
%\note{}
%\section{Warning }{}
\seealso{
  \code{\link{cut}} for converting numeric vectors to factor.
}
\examples{
x <- ordered(sample(letters[1:3],10,rep=TRUE))
cut(x, c(FALSE,TRUE,FALSE))
cut(x, "b")
cut(x, 2)
}
\keyword{manip}
