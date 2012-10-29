\name{p7.11}
\alias{p7.11}
\title{Data Set for Problem 7-11}
\usage{data(p7.11)}
\description{
The \code{p7.11} data frame has 11 observations on production cost
versus production lot size.
}
\format{
  This data frame contains the following columns:
  \describe{
    \item{x}{production lot size}
    \item{y}{average production cost per unit}
  }
}
\source{
 Montgomery, D.C., Peck, E.A., and Vining, C.G. (2001)
Introduction to Linear Regression Analysis. 3rd Edition, John Wiley and Sons.
}
\examples{
data(p7.11)
plot(y ~ x, data=p7.11)
}
\keyword{datasets}
