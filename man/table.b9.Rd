\name{table.b9}
\alias{table.b9}
\non_function{}
\title{Data Set for Table B9}
\usage{data(table.b9)}
\description{
The \code{table.b9} data frame has 62 observations on an
experimental pressure drop.
}
\format{
  This data frame contains the following columns:
  \describe{
    \item{x1}{Superficial fluid velocity of the gas (cm/s)}
    \item{x2}{Kinematic viscosity}
    \item{x3}{Mesh opening (cm)}
    \item{x4}{Dimensionless number relating superficial fluid velocity
of the gas to the superficial fluid velocity of the liquid}
    \item{y}{Dimensionless factor for the pressure drop through a
bubble cap}
  }
}
\source{
 Montgomery, D.C., Peck, E.A., and Vining, C.G. (2001)
Introduction to Linear Regression Analysis. 3rd Edition, John Wiley and Sons.
}
\references{
Liu, C.H., Kan, M., and Chen, B.H. A Correlation of Two-Phase Pressure
Drops in Screen-Plate Bubble Column. Canadian Journal of Chemical 
Engineering, 71, 460-463.
}
\examples{
data(table.b9)
attach(table.b9)
# Partial Solution to Problem 3.13:
y.lm <- lm(y ~ x1 + x2 + x3 + x4)
summary(y.lm)
detach(table.b9)
}
\keyword{datasets}
