\name{table.b14}
\alias{table.b14}
\non_function{}
\title{Table B14}
\usage{data(table.b14)}
\description{
The \code{table.b14} data frame has 25 observations on the transient
points of an electronic inverter.
}
\format{
  This data frame contains the following columns:
  \describe{
    \item{x1}{width of the NMOS Device}
    \item{x2}{length of the NMOS Device}
    \item{x3}{width of the PMOS Device}
    \item{x4}{length of the PMOS Device}
    \item{x5}{a numeric vector}
    \item{y}{transient point of PMOS-NMOS Inverters}
  }
}
\source{
 Montgomery, D.C., Peck, E.A., and Vining, C.G. (2001)
Introduction to Linear Regression Analysis. 3rd Edition, John Wiley and Sons.
}
\examples{
data(table.b14)
y.lm <- lm(y ~ x1 + x2 + x3 + x4, data=table.b14)
plot(y.lm, which=1)
}
\keyword{datasets}
