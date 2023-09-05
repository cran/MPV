\name{table.b23}
\alias{table.b23}
\title{Table B23}
\usage{data(table.b23)}
\description{
The \code{table.b23} data frame has 59 observations on 8 variables.
}
\format{
  This data frame contains the following columns:
  \describe{
    \item{Player}{character vector}
    \item{Per}{numeric vector}
    \item{Lane.Agility.Time..Seconds.}{numeric vector}
    \item{Shuttle.Run..Seconds.}{numeric vector}
    \item{Three.Quarter.Sprint..Seconds.}{numeric vector}
    \item{Standing.Vertical.Leap..Inches.}{numeric vector}
    \item{Max.Vertical.Leap..Inches.}{numeric vector}
    \item{Position}{character vector}
  }
}
\source{
 Montgomery, D.C., Peck, E.A., and Vining, C.G. (2021)
Introduction to Linear Regression Analysis. 6th Edition, John Wiley and Sons.
}
\examples{
pairs(table.b23[,-c(1, 8)])
}
\keyword{datasets}
