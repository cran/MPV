\name{p7.1}
\alias{p7.1}
\non_function{}
\title{Data Set for Problem 7-1}
\usage{data(p7.1)}
\description{
The \code{p7.1} data frame has 10 observations on a predictor variable.
}
\format{
  This data frame contains the following columns:
  \describe{
    \item{x}{a numeric vector}
  }
}
\source{
 Montgomery, D.C., Peck, E.A., and Vining, C.G. (2001)
Introduction to Linear Regression Analysis. 3rd Edition, John Wiley and Sons.
}
\examples{
data(p7.1)
attach(p7.1)
x2 <- x^2
detach(p7.1)
}
\keyword{datasets}
