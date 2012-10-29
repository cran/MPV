\name{p7.15}
\alias{p7.15}
\title{Data Set for Problem 7-15}
\usage{data(p7.15)}
\description{
The \code{p7.15} data frame has 6 observations 
on vapor pressure of water at various temperatures.
}
\format{
  This data frame contains the following columns:
  \describe{
    \item{y}{vapor pressure (mm Hg)}
    \item{x}{temperature (degrees Celsius)}
  }
}
\source{
 Montgomery, D.C., Peck, E.A., and Vining, C.G. (2001)
Introduction to Linear Regression Analysis. 3rd Edition, John Wiley and Sons.
}
\examples{
data(p7.15)
y.lm <- lm(y ~ x, data=p7.15)
plot(y ~ x, data=p7.15)
abline(coef(y.lm))
plot(y.lm, which=1)
}
\keyword{datasets}
