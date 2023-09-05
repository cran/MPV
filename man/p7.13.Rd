\name{p7.13}
\alias{p7.13}
\title{Data Set for Problem 7-13}
\usage{data(p7.13)}
\description{
The \code{p7.13} data frame has 11 observations on production cost
versus production lot size. (This data set was for problem 7-11 in 
the third edition of MPV). 
}
\format{
  This data frame contains the following columns:
  \describe{
    \item{x}{production lot size}
    \item{y}{average production cost per unit}
  }
}
\source{
 Montgomery, D.C., Peck, E.A., and Vining, C.G. (2021)
Introduction to Linear Regression Analysis. 6th Edition, John Wiley and Sons.
}
\examples{
plot(y ~ x, data=p7.13)
}
\keyword{datasets}
