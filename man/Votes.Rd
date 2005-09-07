\name{Votes}
\alias{Votes}
\docType{data}
\title{Congressional Votes 1984 Data Set}
\description{
  This data set includes votes for each of the U.S. House of
  Representatives Congressmen on the 16 key votes identified by the
  CQA.  The CQA lists nine different types of votes: voted for, paired
  for, and announced for (these three simplified to yea), voted
  against, paired against, and announced against (these three
  simplified to nay), voted present, voted present to avoid conflict
  of interest, and did not vote or otherwise make a position known
  (these three simplified to an unknown disposition).
}
\usage{data(Votes)}
\format{
  A data frame with 435 observations on the following 17 variables.
  \describe{
    \item{\code{handicapped-infants}}{a factor with levels \code{n} \code{y}}
    \item{\code{water-project-cost-sharing}}{a factor with levels \code{n} \code{y}}
    \item{\code{adoption-of-the-budget-resolution}}{a factor with levels \code{n} \code{y}}
    \item{\code{physician-fee-freeze}}{a factor with levels \code{n} \code{y}}
    \item{\code{el-salvador-aid}}{a factor with levels \code{n} \code{y}}
    \item{\code{religious-groups-in-schools}}{a factor with levels \code{n} \code{y}}
    \item{\code{anti-satellite-test-ban}}{a factor with levels \code{n} \code{y}}
    \item{\code{aid-to-nicaraguan-contras}}{a factor with levels \code{n} \code{y}}
    \item{\code{mx-missile}}{a factor with levels \code{n} \code{y}}
    \item{\code{immigration}}{a factor with levels \code{n} \code{y}}
    \item{\code{synfuels-corporation-cutback}}{a factor with levels \code{n} \code{y}}
    \item{\code{education-spending}}{a factor with levels \code{n} \code{y}}
    \item{\code{superfund-right-to-sue}}{a factor with levels \code{n} \code{y}}
    \item{\code{crime}}{a factor with levels \code{n} \code{y}}
    \item{\code{duty-free-exports}}{a factor with levels \code{n} \code{y}}
    \item{\code{export-administration-act-south-africa}}{a factor with levels \code{n} \code{y}}
    \item{\code{Class}}{a factor with levels \code{democrat} \code{republican}}
  }
}
\details{
  The records are drawn from:
  
  \emph{Congressional Quarterly Almanac}, 98th Congress, 
  2nd session 1984, Volume XL: Congressional Quarterly Inc. 
  Washington, D.C., 1985.

  It is important to recognize that \code{NA} in this database does 
  not mean that the value of the attribute is unknown.  It 
  means simply, that the value is not "yea" or "nay" (see above).
}
\source{
\url{http://www.ics.uci.edu/~mlearn/MLRepository.html}
}
\references{
  Blake, C.L. \& Merz, C.J. (1998).
  UCI Repository of Machine Learning Databases.
  Irvine, CA: University of California, Department of Information and
  Computer Science.	  
}
\examples{
data(Votes)
summary(Votes)
## maybe str(Votes) ; plot(Votes) ...
}
\keyword{datasets}
