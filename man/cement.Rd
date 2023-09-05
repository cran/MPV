\name{cement}
\alias{cement}
\alias{table.b21}
\title{Table B21 - Cement Data }
\usage{data(cement)}
\description{
The \code{cement} data frame has 13 rows and 5 columns.
}
\format{
  This data frame contains the following columns:
  \describe{
    \item{y}{a numeric vector}
    \item{x1}{a numeric vector}
    \item{x2}{a numeric vector}
    \item{x3}{a numeric vector}
    \item{x4}{a numeric vector}
  }
}
\source{
 Montgomery, D.C., Peck, E.A., and Vining, C.G. (2001)
Introduction to Linear Regression Analysis. 3rd Edition, John Wiley and Sons.
}
\examples{
data(cement)
pairs(cement)
}
\keyword{datasets}
