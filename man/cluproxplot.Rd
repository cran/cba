\name{cluproxplot}
\alias{cluproxplot}
\title{Cluster Proximity Plot}
\description{
  Visualizes cluster quality using an 
  image plot of proximities (a dissimilarity matrix as, 
  e.g., a \code{dist} objects)
  with reordered rows and columns. Reordering is done by seriation
  algorithms which try to place large similarities close to the diagonal.
  Good clusters
  are visible as dark squares (high similarity) on the diagonal of the plot.
  
  The visualization was inspired by CLUSION (see Strehl and Ghosh, 2002).
  }
\usage{
cluproxplot(x, labels = NULL, method = NULL, args = NULL, 
  clusterLabels = TRUE, averages = TRUE, lines = TRUE, silhouettes = TRUE,
  main = "Cluster proximity plot", 
  col = gray.colors(64, 0 , 1), colorkey = TRUE, linesCol = "black", 
  newpage = TRUE, pop = TRUE, \ldots)
}
\arguments{
  \item{x}{ an object of class \code{dist} (distance)
      or a matrix. }
  \item{labels}{ \code{NULL} or an integer vector of the same length
      as rows/columns in \code{x} indicating the membership for each element in
      \code{x} as consecutive integers starting with one.
      The labels are used to reorder the matrix.}
  \item{method}{ a vector of character strings indicating the used
      seriation algorithms. The first element indicates the inter-cluster
      and the second element the intra-cluster seriation method. See 
      \code{\link{seriation}}} for available seriation heuristics.
  \item{args}{\code{"list"}; contains arguments passed on 
      to the seriation algorithms. }
  \item{clusterLabels}{ \code{"logical"}; display cluster labels in the plot.}
  \item{averages}{ \code{"logical"}; display in the lower triangle of the
      plot the average 
      pair-wise dissimilarity instead of the individual dissimilarities.}
  \item{lines}{ \code{"logical"}; draw lines to separate clusters.} 
  \item{silhouettes}{ \code{"logical"}; include a silhouette plot 
	(see Rousseeuw, 1987).} 
  \item{main}{ title for the plot.}
  \item{col}{ colors used for the image plot (default: 64 shades of gray). }
  \item{colorkey}{ place a color key under the plot. }
  \item{linesCol}{ color used for the lines to separate clusters. }
  \item{newpage}{\code{logical"}; start plot on a new page. }
  \item{pop}{\code{logical"};  should the viewports created be popped?  }
  \item{\dots}{ further arguments passed on to \code{image}.}
}
\details{
% To come\ldots
}
\value{
An invisible list of the folowing elements:
\item{order}{ \code{NULL} or integer vector giving the order 
used to plot \code{x}.}
\item{method}{ vector of character strings indicating the seriation methods 
	used for plotting  \code{x}.}
\item{k}{ \code{NULL} or integer scalar giving the number of clusters
generated.}
\item{description}{ a \code{data.frame}  
   containing information (label, size, average intra-cluster dissimilarity
   and the average silhouette) 
   for the clusters as displayed in the 
   plot (from top/left to bottom/right).}
}
\references{
Alexander Strehl and Joydeep Ghosh. Relationship-based clustering and 
visualization for high-dimensional data mining. 
INFORMS Journal on Computing, pages 208-230, Spring 2003.

Rousseeuw, P.J. Silhouettes: A graphical aid to the 
interpretation and validation of cluster analysis. J. Comput. Appl. Math., 
20, 53-65, 1987.

}
\seealso{
  \code{\link[stats]{dist}}
  (in package \code{stats});
   package \code{grid},
   \code{\link{seriation}}.
}
\author{Michael Hahsler (hahsler@ai.wu-wien.ac.at)}
\examples{
data("Votes")

### create dummy coding (with removed party affiliation)
x <- as.dummy(Votes[-17])

### calculate distance matrix
d <- dists(x, method = "binary")

### plot dissimilarity matrix unseriated
res <- cluproxplot(d, method = "No seriation", main = "No seriation")

### plot matrix seriated
res <- cluproxplot(d, main = "Seriation - (Murtagh, 1985)")

### cluster with pam
library("cluster")
l <- pam(d, 8, cluster.only = TRUE)
res <- cluproxplot(d, l, main = "PAM + Seriation (Murtagh)")

res <- cluproxplot(d, l, method = c("Optimal", "Optimal"),
	main = "PAM + Seriation (Optimal Leaf ordering)")

### the result 
names(res)
res$description
}
\keyword{hplot}
\keyword{cluster}
