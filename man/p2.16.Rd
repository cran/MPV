\name{p2.16}
\alias{p2.16}
\non_function{}
\title{Data Set for Problem 2-16}
\usage{data(p2.16)}
\description{
The \code{p2.16} data frame has 33 observations on the
pressure in a tank the volume of liquid.  
}
\format{
  This data frame contains the following columns:
  \describe{
    \item{volume}{volume of liquid}
    \item{pressure}{pressure in the tank}
  }
}
\source{
 Montgomery, D.C., Peck, E.A., and Vining, C.G. (2001)
Introduction to Linear Regression Analysis. 3rd Edition, John Wiley and Sons.
}
\reference{
Carroll and Spiegelman (1986) The Effects of Ignoring Small Measurement
Errors in Precision Instrument Calibration. Journal of Quality Technology,
18, 170-173.  
}
\examples{
data(p2.16)
attach(p2.16)
plot(pressure ~ volume, pch=16)
pressure.lm <- lm(pressure ~ volume)
plot(pressure.lm, which=1)
summary(pressure.lm)
detach(p2.16)
}
\keyword{datasets}
