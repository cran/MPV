\name{table.b24}
\alias{table.b24}
\title{Table B24 - Rental Data}
\usage{data(table.b24)}
\description{
The \code{table.b24} data frame has 51 observations on 6 variables.
}
\format{
  This data frame contains the following columns:
  \describe{
    \item{City}{character vector}
    \item{Population}{numeric vector}
    \item{X95th.Percentile.Income}{numeric vector}
    \item{Median.Sale.Price}{numeric vector}
    \item{Median.Price.sqft}{numeric vector}
    \item{Rental.Price}{numeric vector}
  }
}
\source{
 Montgomery, D.C., Peck, E.A., and Vining, C.G. (2021)
Introduction to Linear Regression Analysis. 6th Edition, John Wiley and Sons.
}
\examples{
pairs(table.b24[,-1])
}
\keyword{datasets}
