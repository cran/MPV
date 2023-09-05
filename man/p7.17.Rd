\name{p7.17}
\alias{p7.17}
\title{Data Set for Problem 7-17}
\usage{data(p7.17)}
\description{
The \code{p7.17} data frame has 6 observations 
on vapor pressure of water at various temperatures. This
data set is the same as \code{p7.15} which was used
for exercise 7-15 in the third edition of MPV. 
}
\format{
  This data frame contains the following columns:
  \describe{
    \item{y}{vapor pressure (mm Hg)}
    \item{x}{temperature (degrees Celsius)}
  }
}
\source{
 Montgomery, D.C., Peck, E.A., and Vining, C.G. (2021)
Introduction to Linear Regression Analysis. 6th Edition, John Wiley and Sons.
}
\examples{
y.lm <- lm(y ~ x, data=p7.17)
plot(y ~ x, data=p7.17)
abline(coef(y.lm))
plot(y.lm, which=1)
}
\keyword{datasets}
