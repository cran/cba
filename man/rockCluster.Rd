\name{rockCluster}
\alias{rockCluster}
\title{Rock Clustering}
\description{
  Cluster a data matrix using the Rock algorithm.
}
\usage{
  rockCluster(x, n, beta = 1-theta, theta = 0.5, fun = "dists", funArgs = list(method="binary"), debug = FALSE)
}
\arguments{
  \item{x}{a data matrix}
  \item{n}{the number of desired clusters}
  \item{beta}{optional distance threshold}
  \item{theta}{neighborhood parameter in the range [0,1)}
  \item{fun}{distance function to use}
  \item{funArgs}{a \code{list} of named parameter arguments to \code{fun}}
  \item{debug}{turn on/off debugging output}
}
\details{
  The intended area of application is the clustering of binary (logical) 
  data. For instance in a preprocessing step in data mining. However, 
  arbitrary distance metrics could be used (see \code{\link{dists}}). 

  According to the reference (see below) the distance threshold and the 
  neighborhood parameter are coupled. Thus, higher values of the neighborhood 
  parameter \code{theta} pose a tighter constraint on the neighborhood. For 
  any two data points the latter is defined as the number of other data points 
  that are neighbors to both. Further, points only are neighbors (or linked) 
  if their distance is less than or equal \code{beta}.

  Note that for a tight neighborhood specification the algorithm may be 
  running out of clusters to merge, i.e. may terminate with more than the 
  desired number of clusters. 
  
  The \code{debug} option can help in determining the proper settings by 
  examining lines suffixed with a plus which indicates that non-singleton 
  clusters were merged. 
  
  Note that tie-breaking is not implemented, i.e. the first max encountered 
  is used. However, permuting the order of the data can help in determining 
  the dependence of a solution on ties.
}
\value{
  An object of class \code{rock}. A list with the following elements:

  \item{x}{the data matrix or a subset of it.}
  \item{cl}{a factor of cluster labels.}
  \item{size}{a vector of cluster sizes.}
  \item{beta}{see above.}
  \item{theta}{see above.}
}
\references{
  S. Guha, R. Rastogi, and K. Shim. ROCK: A Robust Clustering Algorithm for 
  Categorical Attributes. Information Science Vol. 25, No. 5, 2000.
}
\author{Christian Buchta}
\seealso{
  \code{\link{dists}} for common distance functions,
  \code{\link{predict}} for classifying new data samples, and
  \code{\link{fitted}} for classifying the clustered data samples.
}
\examples{
### example from paper
data(Votes)
x <- as.dummy(Votes[-17])
rc <- rockCluster(x, n=2, theta=0.73, debug=TRUE)
print(rc)
rf <- fitted(rc)
table(Votes$Class, rf$cl)
\dontrun{
### large example from paper
data("Mushroom")
x <- as.dummy(Mushroom[-1])
rc <- rockCluster(x[sample(dim(x)[1],1000),], n=10, theta=0.8)
print(rc)
rp <- predict(rc, x)
table(Mushroom$classes, rp$cl)
}
### real valued example
gdist <- function(x, y=NULL) 1-exp(-dists(x, y)^2)
x <- matrix(rnorm(200, sd=0.6)+rep(rep(c(1,-1),each=50),2), ncol=2)
rc <- rockCluster(x, n=2, theta=0.75, fun=gdist, funArgs=NULL)
print(rc)
}
\keyword{cluster}
%