\name{p7.19}
\alias{p7.19}
\non_function{}
\title{Data Set for Problem 7-19}
\usage{data(p7.19)}
\description{
The \code{p7.19} data frame has 10 observations on the concentration
of green liquor and paper machine speed from a kraft paper
machine.
}
\format{
  This data frame contains the following columns:
  \describe{
    \item{y}{green liquor (g/l)}
    \item{x}{paper machine speed (ft/min)}
  }
}
\source{
 Montgomery, D.C., Peck, E.A., and Vining, C.G. (2001)
Introduction to Linear Regression Analysis. 3rd Edition, John Wiley and Sons.
}
\references{
(1986) Tappi Journal.
}
\examples{
data(p7.19)
y.lm <- lm(y ~ x + I(x^2), data=p7.19)
summary(y.lm)
}
\keyword{datasets}
