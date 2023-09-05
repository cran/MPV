\name{p7.18}
\alias{p7.18}
\title{Data Set for Problem 7-18}
\usage{data(p7.18)}
\description{
The \code{p7.18} data frame has 26 observations on the
observed mole fraction solubility of a solute at a
constant temperature. This data set is the same as 
\code{p7.16} which was for problem 7-16 in the third
edition of MPV.
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
pairs(p7.18)
}
\keyword{datasets}
