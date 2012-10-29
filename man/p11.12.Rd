\name{p11.12}
\alias{p11.12}
\title{Data For Problem 11-12}
\usage{data(p11.12)}
\description{
The \code{p11.12} data frame has 19 observations on satellite cost.
}
\format{
  This data frame contains the following columns:
  \describe{
    \item{cost}{first-unit satellite cost}
    \item{x}{weight of the electronics suite}
  }
}
\source{
  Montgomery, D.C., Peck, E.A., and Vining, C.G. (2001)
Introduction to Linear Regression Analysis. 3rd Edition, John Wiley and Sons.
}
\references{
Simpson and Montgomery (1998)
}
\examples{
data(p11.12)
attach(p11.12)
plot(cost~x)
detach(p11.12)
}
\keyword{datasets}
