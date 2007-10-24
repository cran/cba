\name{seriation}
\alias{seriation}
\title{Seriation heuristics for symmetric matrices}
\description{
Seriation heuristics for symmetric matrices used for \code{cluproxplot}. 
Seriation tries to move large values towards the main diagonal by reordering
columns and rows at the same time.
}
\usage{
seriation(x, method = NULL, args = NULL)
}
\arguments{
  \item{x}{ a symmetric matrix to be seriated.}
  \item{method}{ a character string with the name of the seriation method. }
  \item{args}{\code{"list"}; contains arguments for the seriation heuristic.}
}
\details{
    Currently the following methods are implemented:
    
\describe{
    \item{Murtagh (1985), Algorithm B}{a simple heuristic (default method).}
    \item{Hierarchical Clustering}{ uses hierarchical clustering on the matrix; 
      the method can be passed as the element \code{method} in \code{args} 
      (defaults to \code{"average"}). The returned ordering 
      is the ordering of the leaves in the 
      resulting dendrogram.}
    \item{Gruvaeus and Wainer (1972)}{ uses hierarchical clustering and 
      reorders the leaves using a heuristic. }
    \item{Optimal Leaf Ordering}{ uses hierarchical clustering and 
	finds the optimal permutation of the leaves in the dendrogram. }
    }	

}
\value{
The order as a vector of \code{integer}. The attribute \code{"method"} contains
the method as a character string.
}
\references{
F. Murtagh (1985). Multidimensional Cluster Algorithms. \emph{Lectures
in Computational Statistics}, Physica Verlag, pp. 15.

Gruvaeus, G. and Wainer, H. (1972), "Two Additions to
Hierarchical Cluster Analysis", \emph{British
Journal of Mathematical and Statistical Psychology}, 25, 200-206.

Hurley, Catherine B. Clustering Visualizations of Multidimensional Data. 
\emph{Journal of Computational and Graphical Statistics}, 13(4),
pages 788-806, 2004.

Z. Bar-Joseph, E. D. Demaine, D. K. Gifford, and T. Jaakkola.
(2001). Fast Optimal Leaf Ordering for Hierarchical Clustering.
\emph{Bioinformatics}, Vol. 17 Suppl. 1, pp. 22-29.
}
\author{Michael Hahsler (hahsler@ai.wu-wien.ac.at)}
\examples{
data("Votes")

### create dummy coding (with removed party affiliation)
x <- as.dummy(Votes[-17])

### calculate distance matrix
d <- dists(x, method = "binary")


### use standard method
order <- seriation(d)
order

### use hclust w/complete linkage + optimal leaf ordering
order2 <- seriation(d, method = "Optimal")
order2


### use hclust w/single linkage + optimal leaf ordering
order3 <- seriation(d, method = "Optimal", args = list(method = "single"))
order3

}
\keyword{optimize}
