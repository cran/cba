\name{sdists}
\alias{sdists}
\title{Sequence Distance Computation}
\description{
  This function computes and returns the auto-distance matrix between 
  the vectors of a list or between the character strings of a vector 
  treating them as sequences of symbols, as well as the cross-distance 
  matrix between two such lists or vectors.
}
\usage{
sdists(x, y = NULL, method = "ow", weight = c(1, 1, 0, 2), 
       exclude = c(NA, NaN, Inf, -Inf), pairwise = FALSE)
}
\arguments{
  \item{x,y}{a list (of vectors) or a vector of character.}
  \item{method}{a mnemonic string referencing a distance measure.}
  \item{weight}{vector or matrix of parameter values.}
  \item{exclude}{argument to factor.}
  \item{pairwise}{compute distances for the parallel pairs of \code{x}
    and \code{y} only.}
}
\details{
  This function provides a common interface to different methods for 
  computation of distances between sequences, such as the edit a.k.a. 
  Levenshtein distance. Conversely, in the context of sequence alignment 
  the similarity of the maximizing alignment is computed.
  
  Note that negative similarities are returned as distances. So be careful
  to use a proper weighting (scoring) scheme.

  The following methods are currently implemented:

\describe{
    \item{\code{ow}:}{operation-weight edit distance. Weights have to be
	specified for deletion, insertion, match, and replacement. Other
	weights for initial operations can be specified as
	\code{weight[5:6]}.}
    \item{\code{aw}:}{alphabet-weight sequential alignment similarity.
	A matrix of weights (scores) for all possible symbol replacements 
	needs to be specified with the convention that the first row/column 
	defines the replacement with the empty (space) symbol. The colnames 
	of this matrix are used as the levels argument for the encoding as 
	\code{factor}. Consequently, unspecified symbols are mapped to
	\code{NA}.}
    \item{\code{awl}:}{alphabet-weight local sequential alignment similarity.
	The weight matrix must be as described above. However, note that
	zero acts as threshold for a 'restart' of the search for a local
	alignment and at the same time indicates that the solution is the 
	empty substring. Thus, you normally would use non-negative scores
	for matches and non-positive weights otherwise.}
}

  Missing (and non-finite) values should be avoided, i.e. either be removed
  or recoded (and appropriately weighted). By default they are excluded
  when coercing to factor and therefore mapped to \code{NA}. The result 
  is then defined to be \code{NA} as we cannot determine a match!

  The time complexity is O(n*m) for two sequences of length n and m.

  Note that in the case of auto-distances the weight matrix must be
  (exactly) symmetric. Otherwise, for asymmetric weights \code{y}
  must not be \code{NULL}. For instance, \code{x} may be supplied 
  twice (see the examples).
}
\value{
  Auto distances are returned as an object of class \code{dist} and 
  cross-distances as an object of class \code{matrix}. 
}
\references{
    D. Gusfield (1997). \emph{Algorithms on Strings, Trees, and Sequences}. 
    Cambridge University Press, Chapter 11.
}
\author{Christian Buchta}
\section{Warning}{The interface is experimental and may change in the future} 
\seealso{
\code{\link{sdists.trace}} for computation of edit transcripts and sequence alignments,
\code{\link[proxy:dist]{dist}} for computation of common distances,
\code{\link{agrep}} for searches for approximate matches.
}
\examples{
### numeric data
sdists(list(c(2,2,3),c(2,4,3)))			# 2
sdists(list(c(2,2,3),c(2,4,3)),weight=c(1,1,0,1)) # 1

### character data
w <- matrix(-1,nrow=8,ncol=8)			# weight/score matrix for
diag(w) <- 0					# longest common subsequence
colnames(w) <- c("",letters[1:7])
x <- sapply(rbinom(3,64,0.5),function(n,x)
    paste(sample(x,n,rep=TRUE),collapse=""),
    colnames(w)[-1])
x
sdists(x,method="aw",weight=w)
sdists(x,x,method="aw",weight=w)		# check
## pairwise
sdists(x,rev(x),method="aw",weight=w,pairwise = TRUE)
diag(w) <- seq(0,7)
sdists(x,method="aw", weight=w)			# global alignment
sdists(x,method="awl",weight=w)			# local alignment

## empty strings
sdists("", "FOO")
sdists("", list(c("F","O","O")))
sdists("", list(""))				# space symbol
sdists("", "abc", method="aw", weight=w)
sdists("", list(""), method="aw", weight=w)

### asymmetric weights
w[] <- matrix(-sample(0:5,64,TRUE),ncol=8)
diag(w) <- seq(0,7)
sdists(x,x,method="aw", weight=w)
sdists(x,x,method="awl",weight=w)

### missing values
sdists(list(c(2,2,3),c(2,NA,3)),exclude=NULL)	# 2 (include anything)
sdists(list(c(2,2,3),c(2,NA,3)),exclude=NA)	# NA

}
\keyword{cluster}
