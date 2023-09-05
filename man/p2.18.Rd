\name{p2.18}
\alias{p2.18}
\title{Data Set for Problem 2-18}
\usage{data(p2.18)}
\description{
The \code{p2.18} data frame has 21 observations on the
advertising expenses (in millions of US dollars) and retain
impressions (in millions per week)
for various companies.
}
\format{
  This data frame contains the following columns:
  \describe{
    \item{Firm}{character vector}
    \item{Amount.Spent}{numeric vector}
    \item{Returned.Impressions}{numeric vector}
  }
}
\source{
 Montgomery, D.C., Peck, E.A., and Vining, C.G. (2021)
Introduction to Linear Regression Analysis. 6th Edition, John Wiley and Sons.
}
\examples{
data(p2.18)
attach(p2.18)
plot(Returned.Impressions ~ Amount.Spent, pch=16)
detach(p2.18)
}
\keyword{datasets}
