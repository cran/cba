\name{clmplot}
\alias{clmplot}
\title{Plotting Logical Matrices}
\description{
  A wrapper function to image that produces a level plot with the option to 
  color the rows (or columns) based on a clustering and/or classification 
  of the data, and the option to reorder the rows and columns for better
  presentation.
}
\usage{
clmplot(x, col, col.bycol = FALSE, order=FALSE, dists.method = "binary",
        hclust.method = "average", axes = FALSE, xlab = "", ylab = "", ...)
}
\arguments{
  \item{x}{an logical matrix.}
  \item{col}{an optional vector defining a coloring.}
  \item{col.bycol}{option to color by columns.}
  \item{order}{option to (pre)order the rows and columns.}
  \item{dists.method}{method to be used by \code{dists}.}
  \item{hclust.method}{method to be used by \code{hclust}.}
  \item{axes}{option to plot axes.}
  \item{xlab, ylab}{labels for the x and y axis.}
  \item{\dots}{further arguments to \code{image}.}
}
\details{
  For dummy coded data the level \code{FALSE} is assumed to carry no 
  information and is colored "white". Thus, the level \code{TRUE} can 
  be colored according to some classification of the rows (or columns) 
  of the data matrix. If no color specification is supplied the default 
  color \code{black} is used.

  If \code{col} is of type character it is assumed to contain color codes. 
  Otherwise, it must be a factor and the levels are assigned colors
  from \code{\link{heat.colors}}. 
 
  If \code{order} is \code{TRUE} the rows and columns are ordered
  by \code{hclust} where the distances are computed by \code{dists}. 

  Note that an axis is only plotted if the number of elements (rows or
  columns) is less than 100.
}
\value{
  A list with the following components:

  \item{rows}{the row order.}
  \item{cols}{the column order.}
}
\author{Christian Buchta}
\seealso{
  \code{\link{lmplot}} for plotting of logical data at reduced resolutions,
  \code{\link{heatmap}} for ordered plotting of real-valued data, and
    package \code{gclus} for ordering functions.
}
\examples{
x <- matrix(sample(c(FALSE,TRUE),100,rep=TRUE),ncol=10)
clmplot(x, order=TRUE, axes=TRUE)
clmplot(x, col=rep(c(1,2),each=5))
clmplot(x, col=rep(c("red","blue"),each=5))
clmplot(x, col=rep(c("red","blue"),each=5), col.bycol=TRUE)
\dontrun{
### continue example (see rockCluster)
col <- Votes$Class                                  # color by party
levels(col) <- c("red","blue")
op <- par(mfrow=c(1,2), pty="s")
clmplot(x, order=TRUE, col=as.character(col), main="Parties")
col <- rf$cl                                        # color by cluster
levels(col) <- c("blue","red","green", "black")     # map NA to black
clmplot(x, order=TRUE, col=as.character(col), main="Clusters")
par(op)
}
}
\keyword{cluster}
\keyword{hplot}

