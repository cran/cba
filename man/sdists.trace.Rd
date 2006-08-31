\name{sdists.trace}
\alias{sdists.trace}
\title{Edit Transcripts and Sequence Alignments}
\description{
  This function computes and returns the set of all optimal but equivalent
  edit transcripts that transforms one sequences into another at minimum
  cost, as well as the corresponding aligned sequences, or, alternatively 
  a combined edit graph. 
}
\usage{
sdists.trace(x, y, method = "ow", weight = c(1, 0, 2),
	     exclude = c(NA,NaN,Inf,-Inf), graph = FALSE)
}
\arguments{
  \item{x,y}{a numeric or string vector.}
  \item{method}{a mnemonic string referencing a distance measure.}
  \item{weight}{vector or matrix of parameter values.}
  \item{exclude}{argument to factor.}
  \item{graph}{option to compute the combined edit graph.}
}
\details{
  Function \code{sdists.trace} complements the distance computation between
  sequences by \code{\link{sdists}}. So, please, see the details of
  \code{method}, \code{weight}, and \code{exclude} there. However, note the
  following differences: 1) you can supply only two sequences, either as
  vectors of numeric symbol codes, factors, or as strings, i.e. scalar 
  vectors of type \code{character}. 2) you can supply a weight matrix with
  the rownames and colnames representing the symbol sets of the first and
  second sequence. For instance, this allows you to align a sequence with
  the profile of a multiple alignment. 3) if \code{method = "ow"} the
  space symbol \code{""} is included in the factor levels so that you can
  conveniently replace \code{NA} in the aligned sequences.

  A transcript uses the character codes \code{I}, \code{D}, \code{R}, and
  \code{M}, for insert, delete, replace, and match operations, which
  transform the first into the second sequence. Thus, conceptually a symbol
  has to be inserted into the first, deleted from the second, replaced in the
  first sequence, or matched in both, to obtain the second sequence. However,
  in the aligned sequences you will see \code{NA}, where an insert or delete
  would take place, indicating space. 

  In the case of a local alignment different symbols are used for the 
  prefix and/or suffix of the alignment: \code{i}, \code{d}, and \code{?}
  for insert, delete, and replace or match operations. However, note that
  their sole purpose is to obtain a common representation of the two
  sequences. Finally, only alignments of maximal length are reported.

  The time complexity of finding a transcript is O(n+m) for two sequences 
  of length n and m, respectively O(n*m) for the local alignment problem.
  However, note that the runtime for generating all transcripts can be 
  O((n*m)\^3) in the worst case.
}
\value{
  A list with components each a list of two factors, the aligned sequences.
  The names of the components are the edit transcripts, and the attribute
  \code{value} contains the minimum cost, i.e. the distance (or negative
  similarity).
  
  If \code{graph = TRUE} a vector or edit transcripts is returned with
  attributes \code{value}, \code{table}, \code{pointer}, and \code{graph}.
  The second contains the values of the dynamic programming table and the
  third a list of vectors \code{x0, y0, x1, y1} representing the
  (back)pointers. Similarly, the fourth attribute is a list of vectors
  \code{x0, y0, x1, y1, weight} representing the edge set of all optimal
  paths. That is, each tuple contains the \code{from} and \code{to}
  coordinates as used by \code{segments}, each representing a pair of
  indexes into the first and second sequence, and the number of times an
  edge occurs on a path. Note that the origin of the coordinate system
  (0,0) corresponds to the element of \code{table} indexed by ("",""),
  where \code{""} indicates the space symbol. Thus, if used as subscripts
  the coordinates have to be offset by one.
}
\references{
    D. Gusfield (1997) Algorithms on Strings, Trees, and Sequences. 
    Cambridge University Press, Chapter 11.
}
\author{Christian Buchta}
\section{Warning}{The interface is experimental and may change in the future} 
\seealso{
\code{\link{sdists}} for computation of distances between sequences,
\code{\link{segments}} for plotting of edge sets,
\code{\link{plot.sdists.graph}} for visualizing alignments.
}
\examples{
### from the book
x1 <- "vintner"
y1 <- "writers"
b1 <- sdists.trace(x1, y1, weight=c(1,0,1))
b1
## longest common subsequence ?
sdists.trace("a","b", weight=c(0,-1,0))
## from the book
w2 <- matrix(-2,ncol=13,nrow=13)
w2[1,] <- w2[,1] <- -1
diag(w2) <- c(0,rep(2,12))
x2 <- "pqraxabcstvq"
y2 <- "xyaxbacsll"
colnames(w2) <- c("",unique(strsplit(paste(x2, y2, sep = ""),"")[[1]]))
b2 <- sdists.trace(x2, y2, method="awl", weight=w2)
b2
## alignment with different symbol sets
x3 <- "121314"
y3 <- "ABACAD"
w3 <- matrix(-1,nrow=5,ncol=5)
diag(w3) <- 0
rownames(w3) <- c("","1","2","3","4")
colnames(w3) <- c("","A","B","C","D")
b3 <- sdists.trace(x3, y3, method="aw", weight=w3)
b3
### test graphics computations
b11 <- sdists.trace(x1, y1, weight=c(1,0,1), graph = TRUE)
b11
}
\keyword{cluster}
