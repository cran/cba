\name{proximus}
\alias{proximus}
\title{Proximus}
\description{
  Cluster the rows of a logical matrix using the Proximus algorithm. The 
  compression rate of the algorithm can be influenced by the choice of 
  the maximum cluster radius and the minimum cluster size. 
}
\usage{
proximus(x, max.radius = 2, min.size = 1, min.retry = 10, max.iter = 16,
         debug = FALSE)
}
\arguments{
  \item{x}{a logical matrix.}
  \item{max.radius}{the maximum number of bits a member in a row set may 
    deviate from its dominant pattern.}
  \item{min.size}{the minimum split size of a row set.}
  \item{min.retry}{number of retries to split a pure rank-one approximation 
    (translates into a resampling rate).}
  \item{max.iter}{the maximum number of iterations for finding a local 
    rank-one approximation.}
  \item{debug}{optional debugging output.}
}
\details{
  The intended area of application is the compression of high-dimensional 
  binary data into representative patterns. For instance, purchase incidence 
  (market basket date) or term-document matrices may be preprocessed by 
  Proximus for later association rule mining.

  The algorithm is of a recursive partitioning type. 
  Specifically, at each step a binary split is attempted using a local 
  rank-one approximation of the current submatrix (row set). That is a 
  specialization of principal components to binary data which represents a 
  matrix as the outer product of two binary vectors. The node expansion stops 
  if a submatrix is pure, i.e. the column (presence set) vector indicates 
  all the rows and the Hamming distances from the row (dominant attribute 
  set) pattern vector, or the size of the row set, are equal or less than 
  the specified threshold. In the case the rank-one approximation does not 
  result in a split but the radius constraint is violated, the matrix is 
  split using a random row and the radius constraint.

  The debug option can be used to gain some insight into how the algorithm 
  proceeds: a right angle bracket indicates a split and the return to
  a recursion level is indicated by a left one. Leafs in the recursion tree 
  are indicated by an asterisk and retries by a plus sign. The number of 
  retries is bounded by the size of the current set divided by 
  \code{min.retry}. 
  Double angle brackets indicate a random split (see above). The numbers 
  between square brackets indicate the current set size, the size of the 
  presence (sub)set, and its radius. The adjoining numbers indicate the 
  depth of the recursion and the count of retries. Finally, a count of 
  the leaf nodes found so far is shown to the right of an asterisk.
}
\value{
  An object of class \code{proximus} with the following components:

  \item{nr}{the number of rows of the data matrix.}
  \item{nc}{the number of columns of the data matrix.}
  \item{a}{a list containing the approximations (patterns).}
  \item{a$x}{a vector of row (presence set) indexes.}
  \item{a$y}{a vector of column (dominant attribute set) indexes.}
  \item{a$n}{the number of ones in the approximated submatrix.}
  \item{a$c}{the absolute error reduction by the approximation.}
  \item{max.radius}{see arguments.}
  \item{min.size}{see arguments.}
  \item{rownames}{rownames of the data matrix.}
  \item{colnames}{colnames of the data matrix.}
}
\references{
  M. Koyutürk, A. Graham, and N. Ramakrishnan. Compression, Clustering, and 
  Pattern Discovery in Very High-Dimensional Discrete-Attribute Data Sets. 
  IEEE Transactions On Knowledge and Data Engineering, Vol. 17, No. 4, 
  (April) 2005.
}
\author{Christian Buchta}
\note{
  The size of a set need not be equal or greater than the user defined 
  threshold.
}
\section{Warning}{Deep recursions may exhaust your computer.}
\seealso{
  \code{\link{summary.proximus}} for summaries,
  \code{\link{fitted}} for obtaining the approximated matrix and the 
    pattern labels of the samples, and
  \code{\link{lmplot}} for plotting logical matrices.
}
\examples{
x <- matrix(sample(c(FALSE, TRUE), 200, rep=TRUE), ncol=10)
pr <- proximus(x, max.radius=8)
summary(pr)
### example from paper
x <- rlbmat()
pr <- proximus(x, max.radius=8, debug=TRUE)
op <- par(mfrow=c(1,2), pty="s")
lmplot(x, main="Data")
box()
lmplot(fitted(pr)$x, main="Approximation")
box()
par(op)
}
\keyword{cluster}

