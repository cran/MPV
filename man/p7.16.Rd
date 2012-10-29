\name{p7.16}
\alias{p7.16}
\title{Data Set for Problem 7-16}
\usage{data(p7.16)}
\description{
The \code{p7.16} data frame has 26 observations on the
observed mole fraction solubility of a solute at a
constant temperature.
}
\format{
  This data frame contains the following columns:
  \describe{
    \item{y}{negative logarithm of the mole fraction solubility}
    \item{x1}{dispersion partial solubility}
    \item{x2}{dipolar partial solubility}
    \item{x3}{hydrogen bonding Hansen partial solubility}
  }
}
\source{
 Montgomery, D.C., Peck, E.A., and Vining, C.G. (2001)
Introduction to Linear Regression Analysis. 3rd Edition, John Wiley and Sons.
}
\references{
(1991) Journal of Pharmaceutical Sciences 80, 971-977.
}
\examples{
data(p7.16)
pairs(p7.16)
}
\keyword{datasets}
