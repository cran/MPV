\name{p9.10}
\alias{p9.10}
\title{Data Set for Problem 9-10}
\usage{data(p9.10)}
\description{
The \code{p9.10} data frame has 31 observations
on the rut depth of asphalt pavements prepared under
different conditions.
}
\format{
  This data frame contains the following columns:
  \describe{
    \item{y}{change in rut depth/million wheel passes (log scale)}
    \item{x1}{viscosity (log scale)}
    \item{x2}{percentage of asphalt in surface course}
    \item{x3}{percentage of asphalt in base course}
    \item{x4}{indicator}
    \item{x5}{percentage of fines in surface course}
    \item{x6}{percentage of voids in surface course}
  }
}
\source{
 Montgomery, D.C., Peck, E.A., and Vining, C.G. (2001)
Introduction to Linear Regression Analysis. 3rd Edition, John Wiley and Sons.
}
\references{
Gorman and Toman (1966) 
}
\examples{
data(p9.10)
pairs(p9.10)
}
\keyword{datasets}
