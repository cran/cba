\name{cluproxplot}
\alias{cluproxplot}
\title{Cluster Proximity Plot}
\description{
Visualizes cluster quality using shading of a rearranged proximity  matrix (see
Ling, 1973). Objects belonging to the same cluster are displayed in consecutive
order. The placement of clusters and the within cluster order is done by
various seriation algorithms which try to place large similarities close to
the diagonal. Compact clusters are visible as dark squares (high similarity)
on the diagonal of the plot.

Additionally a Silhouette plot (Rousseeuw, 1987) is added.

The visualization was also inspired by CLUSION (see Strehl and Ghosh, 2002).
  }
\usage{
cluproxplot(x, labels = NULL, method = NULL, args = NULL, 
            plot = TRUE, plotOptions = NULL, \ldots)   
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
  \item{plot}{\code{logical"};  if \code{FALSE}, no plot is produced. The 
    returned object can be plotted later using the function \code{plot}
    which takes as the second argument a a list of plotting options (see
    \code{plotOptions} below).}
  \item{plotOptions}{\code{"list"}; options for plotting the matrix. The list
    can contain the following elements:
    \describe{
  \item{clusterLabels}{ \code{"logical"}; display cluster labels in the plot.}
  \item{averages}{ \code{"logical"}; display in the lower triangle of the
    plot the average 
    pair-wise dissimilarity between clusters 
    instead of the individual dissimilarities.}
  \item{lines}{ \code{"logical"}; draw lines to separate clusters.} 
  \item{silhouettes}{ \code{"logical"}; include a silhouette plot 
    (see Rousseeuw, 1987).} 
  \item{threshold}{ \code{"numeric"}; only plot distances below the 
    threshold.} 
  \item{main}{ title for the plot.}
  \item{col}{ colors used for the image plot (default: 100 shades of gray using
    the hcl colorspace with \code{hcl(h = 0, c = 0, l = seq(5, 95, len = 100))}). }
  \item{colorkey}{ place a color key under the plot. }
  \item{linesCol}{ color used for the lines to separate clusters. }
  \item{newpage}{\code{logical"}; start plot on a new page 
    (see package \code{grid}). }
  \item{pop}{\code{logical"};  should the viewports created be popped 
    (see package \code{grid})?  }
    } 
    }
  \item{\dots}{ further arguments; currently unused. }
}
\details{
% To come\ldots
}
\value{
An invisible object of class \code{"cluProxMatrix"} of the following elements:
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
Ling, R.F. A computer generated aid for cluster analysis. Comm. of the ACM, 
16(6), 355-361, 1973.

Rousseeuw, P.J. Silhouettes: A graphical aid to the 
interpretation and validation of cluster analysis. J. Comput. Appl. Math., 
20, 53-65, 1987.

Strehl, A. and Ghosh, J. Relationship-based clustering and 
visualization for high-dimensional data mining. 
INFORMS Journal on Computing, 208-230, 2003.

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
res <- cluproxplot(d, method = "No seriation", 
	plotOptions = list(main = "No seriation"))

### plot matrix seriated
res <- cluproxplot(d, plotOptions = list(main = "Seriation - (Murtagh, 1985)"))

### cluster with pam
library("cluster")
l <- pam(d, 8, cluster.only = TRUE)
res <- cluproxplot(d, l, plotOptions = list(main = "PAM + Seriation (Murtagh)"))

### now we use a different seriation algorithm (hclust + optimal leaf ordering)
### and just do the seriation and then use plot to produce the plot
res <- cluproxplot(d, l, method = c("Optimal", "Optimal"), plot = FALSE)
res


### use blue (hue is 260 with decreasing chroma and  increasing luminance 
### towards a distance of 1)
plot(res, plotOptions = list(main = "PAM + Seriation (Optimal Leaf ordering)", 
	col = hcl(h = 260, c = seq(75,0, length=5), l = seq(30,95, length=5))))

### the result contains more information, e.g., the order used for reordering
### the matrix
names(res)
res$order
}
\keyword{hplot}
\keyword{cluster}
