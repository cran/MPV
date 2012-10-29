\name{p11.15}
\alias{p11.15}
\title{Data set for Problem 11-15}
\usage{data(p11.15)}
\description{
The \code{p11.15} data frame has 9 rows and 2 columns.
}
\format{
  This data frame contains the following columns:
  \describe{
    \item{x}{a numeric vector}
    \item{y}{a numeric vector}
  }
}
\source{
 Montgomery, D.C., Peck, E.A., and Vining, C.G. (2001)
Introduction to Linear Regression Analysis. 3rd Edition, John Wiley and Sons.
}
\references{
Ryan (1997), Stefanski (1991)
}
\examples{
data(p11.15)
plot(p11.15)
attach(p11.15)
lines(lowess(x,y))
detach(p11.15)
}
\keyword{datasets}
