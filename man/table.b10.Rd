\name{table.b10}
\alias{table.b10}
\title{Table B10}
\usage{data(table.b10)}
\description{
The \code{table.b10} data frame has 40 observations
on kinematic viscosity of a certain solvent system.
}
\format{
  This data frame contains the following columns:
  \describe{
    \item{x1}{Ratio of 2-methoxyethanol to 1,2-dimethoxyethane}
    \item{x2}{Temperature (in degrees Celsius)}
    \item{y}{Kinematic viscosity (.000001 m2/s}
  }
}
\source{
 Montgomery, D.C., Peck, E.A., and Vining, C.G. (2001)
Introduction to Linear Regression Analysis. 3rd Edition, John Wiley and Sons.
}
\references{
Viscosimetric Studies on 2-Methoxyethanol + 1, 2-Dimethoxyethane
Binary Mixtures from -10 to 80C. Canadian Journal of Chemical
Engineering, 75, 494-501.
}
\examples{
data(table.b10)
attach(table.b10)
y.lm <- lm(y ~ x1 + x2)
summary(y.lm)
detach(table.b10)
}
\keyword{datasets}
