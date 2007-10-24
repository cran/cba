\name{summary.proximus}
\alias{summary.proximus}
\alias{print.summary.proximus}
\title{Summarizing Proximus Objects}
\description{
  \code{summary} method for an object of class \code{proximus}.
}
\usage{
\method{summary}{proximus}(object, ...)
}
\arguments{
  \item{object}{an object of class \code{proximus}.}
  \item{\dots}{further (unused) arguments.}
}
\value{
  An object of class \code{summary.proximus} with the following elements:

  \item{nr}{the number of rows of the approximated matrix.}
  \item{nc}{the number of columns of the approximated matrix.}
  \item{error}{the relative error of the total approximation.}
  \item{fnorm}{the Frobenius norm of the total approximation.}
  \item{jsim}{the Jaccard similarity of the total approximation.}
  \item{valid}{the number of patterns that satisfy the mining constraints.}
  \item{pattern}{a \code{data.frame} of pattern summaries:}
  \item{pattern$Size}{the absolute size of the presence set.}
  \item{pattern$Length}{the number of \code{TRUE} values of the dominant 
    pattern.}
  \item{pattern$Radius}{the Hamming radius of the presence set.}
  \item{pattern$Error}{the relative error of the presence set.}
  \item{pattern$Fnorm}{the Frobenius norm of the presence set.}
  \item{pattern$Jsim}{the Jaccard similarity of the presence set.}
  \item{pattern$Valid}{a logical indicating if the constraints are satisfied.}
}
\author{Christian Buchta}
\section{Warning}{The function may change in future releases}
\seealso{
  \code{\link{proximus}} for obtaining a result object.
}
\examples{
### see proximus
}
\keyword{cluster}

