\name{p4.20}
\alias{p4.20}
\non_function{}
\title{Data Set for Problem 4-20
}
\usage{data(p4.20)}
\description{
The \code{p4.20} data frame has 26 observations
on a designed experiment to determine the influence
of five factors on the whiteness of rayon.
}
\format{
  This data frame contains the following columns:
  \describe{
    \item{acidtemp}{acid bath temperature}
    \item{acidconc}{cascade acid concentration}
    \item{watertemp}{water temperature}
    \item{sulfconc}{sulfide concentration}
    \item{amtbl}{amount of chlorine bleach}
    \item{y}{a measure of the whiteness of rayon}
  }
}
\source{
 Montgomery, D.C., Peck, E.A., and Vining, C.G. (2001)
Introduction to Linear Regression Analysis. 3rd Edition, John Wiley and Sons.
}
\references{
Myers and Montgomery (1995) Response Surface Methodology, pp. 267-268.
}
\examples{
data(p4.20)
y.lm <- lm(y ~ acidtemp, data=p4.20)
summary(y.lm)
}
\keyword{datasets}
