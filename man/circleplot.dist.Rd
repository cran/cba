\name{circleplot.dist}
\alias{circleplot.dist}
\title{Plotting Distance Graphs}
\description{
  Function for visualizing distance graphs using a circular layout.
}
\usage{
circleplot.dist(x, cutoff = 0.5, col = 1, circle = FALSE, scale = 1.4)
}
\arguments{
  \item{x}{an object of class \code{dist}.}
  \item{cutoff}{a numeric value specifying the threshold for edge exclusion.}
  \item{col}{a number or string specifying the edge color to use.}
  \item{circle}{a logical value specifying if a circle connecting the
    nodes should be drawn.}
  \item{scale}{a numeric value specifying the plot range (the default
    accommodates node labels).}
}
\details{
  Plots the distance graph of \code{x} placing its nodes on a circle
  such that the number of crossing edges is approximately minimized.
  This is achieved by using \code{order.dist} for seriation.
}
%\value{}
%\references{}
\author{Christian Buchta}
%\note{}
\seealso{\code{\link{order.dist}}.}
\examples{
##
data(iris)
d <- dist(iris[,-5])[[1:26]]
circleplot.dist(d, col = 2, scale = 1)

dimnames(d) <- LETTERS[1:26]
circleplot.dist(d)
}
\keyword{hplot}
\keyword{cluster}
