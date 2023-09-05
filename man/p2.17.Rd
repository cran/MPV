\name{p2.17}
\alias{p2.17}
\title{Data Set for Problem 2-17}
\usage{data(p2.17)}
\description{
The \code{p2.17} data frame has 17 observations on the
boiling point of water (in Fahrenheit degrees)
for various barometric pressures (in inches of mercury).
}
\format{
  This data frame contains the following columns:
  \describe{
    \item{BoilingPoint}{numeric vector}
    \item{BarometricPressure}{numeric vector}
  }
}
\source{
 Montgomery, D.C., Peck, E.A., and Vining, C.G. (2021)
Introduction to Linear Regression Analysis. 6th Edition, John Wiley and Sons.
}
\references{
Atkinson, A.C. (1985) Plots, Transformations and Regression, Clarendon Press,
Oxford.
}
\examples{
data(p2.17)
attach(p2.17)
plot(BoilingPoint ~ BarometricPressure, pch=16)
detach(p2.17)
}
\keyword{datasets}
