\name{p4.19}
\alias{p4.19}
\non_function{}
\title{Data Set for Problem 4-19}
\usage{data(p4.19)}
\description{
The \code{p4.19} data frame has 14 observations on
a designed experiment studying the relationship
between abrasion index for a tire tread compound
and three factors.
}
\format{
  This data frame contains the following columns:
  \describe{
    \item{x1}{hydrated silica level}
    \item{x2}{silane coupling agent level}
    \item{x3}{sulfur level}
    \item{y}{abrasion index for a tire tread compound}
  }
}
\source{
 Montgomery, D.C., Peck, E.A., and Vining, C.G. (2001)
Introduction to Linear Regression Analysis. 3rd Edition, John Wiley and Sons.
}
\references{
Derringer and Suich (1980) Simultaneous Optimization of Several
Response Variables. Journal of Quality Technology.
}
\examples{
data(p4.19)
attach(p4.19)
y.lm <- lm(y ~ x1 + x2 + x3)
summary(y.lm)
plot(y.lm, which=1)
y.lm <- lm(y ~ x1)
detach(p4.19)
}
\keyword{datasets}
