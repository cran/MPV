\name{p12.11}
\alias{p12.11}
\title{Data Set for Problem 12-11}
\usage{data(p12.11)}
\description{
The \code{p12.11} data frame has 44 observations on the fraction
of active chlorine in a chemical product as a function of time
after manufacturing.
}
\format{
  This data frame contains the following columns:
  \describe{
    \item{xi}{time}
    \item{yi}{available chlorine}
  }
}
\source{
 Montgomery, D.C., Peck, E.A., and Vining, C.G. (2001)
Introduction to Linear Regression Analysis. 3rd Edition, John Wiley and Sons.
}
\examples{
data(p12.11)
plot(p12.11)
lines(lowess(p12.11))
}
\keyword{datasets}
