\name{table.b22}
\alias{table.b22}
\title{Table B22 - Baseball Data}
\usage{data(table.b22)}
\description{
The \code{table.b22} data frame has 30 observations on 12 variables.
}
\format{
  This data frame contains the following columns:
  \describe{
    \item{Team}{character vector}
    \item{Wins}{numeric vector}
    \item{Batter.Age}{numeric vector}
    \item{Runs}{numeric vector}
    \item{HRs}{numeric vector}
    \item{SLG}{numeric vector}
    \item{Pitcher.Age}{numeric vector}
    \item{ERA}{numeric vector}
    \item{SO}{numeric vector}
    \item{HRA}{numeric vector}
    \item{RA.G}{numeric vector}
    \item{Errors}{numeric vector}
  }
}
\source{
 Montgomery, D.C., Peck, E.A., and Vining, C.G. (2021)
Introduction to Linear Regression Analysis. 6th Edition, John Wiley and Sons.
}
\examples{
pairs(table.b22[,-1])
}
\keyword{datasets}
