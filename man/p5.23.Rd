\name{p5.23}
\alias{p5.23}
\title{Data Set for Problem 5-23}
\usage{data(p5.23)}
\description{
The \code{p5.23} data frame has 18 observations on
3 variables.
}
\format{
  This data frame contains the following columns:
  \describe{
    \item{Batch}{a character vector}
    \item{Pressure}{a numeric vector}
    \item{Strength}{a numeric vector}
  }
}
\source{
 Montgomery, D.C., Peck, E.A., and Vining, C.G. (2021)
Introduction to Linear Regression Analysis. 6th Edition, John Wiley and Sons.
}
\examples{
 anova(lm(Strength ~ Pressure*Batch, data = p5.23))
}
\keyword{datasets}
