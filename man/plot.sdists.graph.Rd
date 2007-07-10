\name{plot.sdists.graph}
\alias{plot.sdists.graph}
\title{Plotting Edit Transcripts and Sequence Alignments}
\description{
  Function for visualizing the optimal transformation (alignment) graph
  for two symbol sequences.
}
\usage{
\method{plot}{sdists.graph}(x, circle.col = 1, graph.col = 2,
     circle.scale = c("mean", "max", "last"), main = "", ...)
}
\arguments{
  \item{x}{an object of class \code{sdists.graph.}}
  \item{circle.col}{color to be used for circles.}
  \item{graph.col}{color to be used for the graph.}
  \item{circle.scale}{scaling to be used for circles.}
  \item{main}{plot title.}
  \item{\dots}{further unused arguments.}
}
\details{
  This function plots the dynamic programming table, the (back)pointers
  and the combined graph of optimal edit transcripts (alignments) computed
  with \code{sdist.trace}. The fist sequence is represented by the y-axis
  and the second by the x-axis.

  The circumference of a circle is proportional to the minimum cost
  (maximum weight) of edit (alignment) operations leading to a table cell.
  \code{circle.scale} specifies the type of normalization performed where
  \code{last} means the last table entry (containing the optimum value).
  
  The (back)pointers, defining possible (traceback) paths, are plotted
  as \code{dotted} edges. Note that a traceback starts in the last cell
  of the table and ends at the origin.
  
  For the edges of the graph that is the union of all optimal paths, two
  line types are used: \code{solid} for insert, delete, and replace
  operations, and \code{dashed} for a match. The line width indicates the
  number of times an edge is on a path, but note that the interpretation
  is device-specific (compare \code{\link{par}}).
}
%\value{
%}
\references{
  D. Gusfield (1997). \emph{Algorithms on Strings, Trees, and Sequences}. 
  Cambridge University Press, Chapter 11.
	  
  Inspired by: \url{http://home.uchicago.edu/~aabbott/}.
}
\author{Christian Buchta}
\note{
  Some issues with grid were fixed in R.2.4.x (Fixme ?).
}
\seealso{\code{\link{sdists.trace}}, \code{\link{plot}}}
\examples{
## continue example in sdists.trace
x1 <- "vintner"
y1 <- "writers"
b11 <- sdists.trace(x1, y1, weight=c(1,0,1), graph = TRUE)
b11
plot(b11)
}
\keyword{cluster}
