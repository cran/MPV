\name{p8.3}
\alias{p8.3}
\title{Data Set for Problem 8-3}
\usage{data(p8.3)}
\description{
The \code{p8.3} data frame has 25 observations on delivery
times taken by a vending machine route driver.
}
\format{
  This data frame contains the following columns:
  \describe{
    \item{y}{delivery time (in minutes)}
    \item{x1}{number of cases of product stocked}
    \item{x2}{distance walked by route driver}
  }
}
\source{
 Montgomery, D.C., Peck, E.A., and Vining, C.G. (2001)
Introduction to Linear Regression Analysis. 3rd Edition, John Wiley and Sons.
}
\examples{
data(p8.3)
pairs(p8.3)
}
\keyword{datasets}
