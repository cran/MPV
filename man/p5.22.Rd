\name{p5.22}
\alias{p5.22}
\title{Data Set for Problem 5-22}
\usage{data(p5.22)}
\description{
The \code{p5.22} data frame has 18 observations on
2 variables.
}
\format{
  This data frame contains the following columns:
  \describe{
    \item{Temp}{a numeric vector}
    \item{Density}{a numeric vector}
  }
}
\source{
 Montgomery, D.C., Peck, E.A., and Vining, C.G. (2021)
Introduction to Linear Regression Analysis. 6th Edition, John Wiley and Sons.
}
\examples{
 anova(lm(Density ~ Temp, data = p5.22))
}
\keyword{datasets}
