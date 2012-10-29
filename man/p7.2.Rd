\name{p7.2}
\alias{p7.2}
\title{Data Set for Problem 7-2}
\usage{data(p7.2)}
\description{
The \code{p7.2} data frame has 10 observations on solid-fuel
rocket propellant weight loss.
}
\format{
  This data frame contains the following columns:
  \describe{
    \item{x}{months since production}
    \item{y}{weight loss (kg)}
  }
}
\source{
 Montgomery, D.C., Peck, E.A., and Vining, C.G. (2001)
Introduction to Linear Regression Analysis. 3rd Edition, John Wiley and Sons.
}
\examples{
data(p7.2)
y.lm <- lm(y ~ x + I(x^2), data=p7.2)
summary(y.lm)
plot(y ~ x, data=p7.2)
}
\keyword{datasets}
