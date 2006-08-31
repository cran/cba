\name{townships}
\alias{townships}
\docType{data}
\title{Bertin's Characteristics and Townships Data Set}
\description{
  This data set was used to illustrate that the conciseness of presentation
  can be improved by reordering the rows and columns.
}
\usage{data(townships)}
\format{
  A data frame with 16 observations on the following 10 variables.
  \describe{
    \item{\code{Township}}{a factor with levels \code{A} \code{B} \code{C} \code{D} \code{E} \code{F} \code{G} \code{H} \code{I} \code{J} \code{K} \code{L} \code{M} \code{N} \code{O} \code{P}}
    \item{\code{High.School}}{a logical vector}
    \item{\code{Agricultural.Coop.}}{a logical vector}
    \item{\code{Railway.Station}}{a logical vector}
    \item{\code{One.Room.School}}{a logical vector}
    \item{\code{Veterinary}}{a logical vector}
    \item{\code{No.Doctor}}{a logical vector}
    \item{\code{No.Water.Supply}}{a logical vector}
    \item{\code{Police.Station}}{a logical vector}
    \item{\code{Land.Reallocation}}{a logical vector}
  }
}
\details{
  \code{townships} is a data set with 16 \code{logical} variables
  indicating the presence (\code{TRUE}) or absence (\code{FALSE})
  of characteristics of townships.
}
%\source{}
\references{
  Bertin, J. (1981) \emph{Graphics and Graphic Information Processing}.
  Berlin, Walter de Gruyter.
}
\examples{
## see order.data.frame
}
\keyword{datasets}
