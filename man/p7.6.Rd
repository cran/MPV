\name{p7.6}
\alias{p7.6}
\non_function{}
\title{Data Set for Problem 7-6}
\usage{data(p7.6)}
\description{
The \code{p7.6} data frame has 12 observations on softdrink
carbonation.
}
\format{
  This data frame contains the following columns:
  \describe{
    \item{y}{carbonation}
    \item{x1}{temperature}
    \item{x2}{pressure}
  }
}
\source{
 Montgomery, D.C., Peck, E.A., and Vining, C.G. (2001)
Introduction to Linear Regression Analysis. 3rd Edition, John Wiley and Sons.
}
\examples{
data(p7.6)
y.lm <- lm(y ~ x1 + I(x1^2) + x2 + I(x2^2) + I(x1*x2), data=p7.6)
summary(y.lm)
}
\keyword{datasets}
