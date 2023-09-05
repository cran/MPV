\name{table.b25}
\alias{table.b25}
\title{Table B25 Golf Data}
\usage{data(table.b25)}
\description{
The \code{table.b25} data frame has 50 observations on 6 variables.
}
\format{
  This data frame contains the following columns:
  \describe{
    \item{Player}{character vector}
    \item{Average.Score}{numeric vector}
    \item{SG..Off.the.Tee}{numeric vector}
    \item{SG..Approach.to.Green}{numeric vector}
    \item{SG..Around.the.Green}{numeric vector}
    \item{SG..Putting}{numeric vector}
  }
}
\source{
 Montgomery, D.C., Peck, E.A., and Vining, C.G. (2021)
Introduction to Linear Regression Analysis. 6th Edition, John Wiley and Sons.
}
\examples{
pairs(table.b25[,-1])
}
\keyword{datasets}
