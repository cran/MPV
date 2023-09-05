\name{p5.24}
\alias{p5.24}
\title{Data Set for Problem 5-24}
\usage{data(p5.24)}
\description{
The \code{p5.24} data frame has 13 observations on
7 variables.
}
\format{
  This data frame contains the following columns:
  \describe{
    \item{Location}{a character vector}
    \item{x1}{a numeric vector}
    \item{x2}{a numeric vector}
    \item{x3}{a numeric vector}
    \item{x4}{a numeric vector}
    \item{x5}{a numeric vector}
    \item{y}{a numeric vector}
  }
}
\source{
 Montgomery, D.C., Peck, E.A., and Vining, C.G. (2021)
Introduction to Linear Regression Analysis. 6th Edition, John Wiley and Sons.
}
\references{
French, R.J. and Schultz, J.E. "Water Use Efficiency of Wheat in a Mediterranean-type Environment, I The Relation between Yield, Water Use, and Climate," Australian Journal of Agricultural Research, 35, 743-764, 1984.
}
\examples{
 lm(y ~ x1 + x2 + x3 + x4 + x5, data = p5.24)
}
\keyword{datasets}
