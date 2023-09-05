\name{p7.20}
\alias{p7.20}
\title{Data Set for Problem 7-20}
\usage{data(p7.20)}
\description{
The \code{p7.20} data frame has 10 observations on the concentration
of green liquor and paper machine speed from a kraft paper
machine.This data set is the same as \code{p7.19} which was used
in problem 7.19 of the third edition of MPV.
}
\format{
  This data frame contains the following columns:
  \describe{
    \item{y}{green liquor (g/l)}
    \item{x}{paper machine speed (ft/min)}
  }
}
\source{
 Montgomery, D.C., Peck, E.A., and Vining, C.G. (2021)
Introduction to Linear Regression Analysis. 6th Edition, John Wiley and Sons.
}
\references{
(1986) Tappi Journal.
}
\examples{
data(p7.20)
y.lm <- lm(y ~ x + I(x^2), data=p7.20)
summary(y.lm)
}
\keyword{datasets}
