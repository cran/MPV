\name{table.b2}
\alias{table.b2}
\title{Table B2}
\usage{data(table.b2)}
\description{
The \code{table.b2} data frame contains 29 observations on 6 variables 
related to a solar thermal energy test.
}
\format{
  This data frame contains the following columns:
  \describe{
    \item{y}{a numeric vector measuring total heat flux (kwatts)}
    \item{x1}{a numeric vector measuring insulation (watts/m^2)}
    \item{x2}{a numeric vector measuring position of focal point in east direction (inches)}
    \item{x3}{a numeric vector measuring position of focal point in south direction (inches)}
    \item{x4}{a numeric vector measuring position of focal point in north direction (inches)}
    \item{x5}{a numeric vector representing time of day}
  }
}
\source{
 Montgomery, D.C., Peck, E.A., and Vining, C.G. (2001)
Introduction to Linear Regression Analysis. 3rd Edition, John Wiley and Sons.
}
\examples{
data(table.b2)
pairs(table.b2)
}
\keyword{datasets}
