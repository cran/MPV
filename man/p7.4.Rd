\name{p7.4}
\alias{p7.4}
\non_function{}
\title{Data Set for Problem 7-4}
\usage{data(p7.4)}
\description{
The \code{p7.4} data frame has 12 observations on two variables.
}
\format{
  This data frame contains the following columns:
  \describe{
    \item{x}{a numeric vector}
    \item{y}{a numeric vector}
  }
}
\source{
 Montgomery, D.C., Peck, E.A., and Vining, C.G. (2001)
Introduction to Linear Regression Analysis. 3rd Edition, John Wiley and Sons.
}
\examples{
data(p7.4)
y.lm <- lm(y ~ x + I(x^2), data = p7.4)
summary(y.lm)
}
\keyword{datasets}
