\name{p2.7}
\alias{p2.7}
\title{Data Set for Problem 2-7}
\usage{data(p2.7)}
\description{
The \code{p2.7} data frame has 20 observations on the
purity of oxygen produced by a fractionation process.  It
is thought that oxygen purity is related to the percentage
of hydrocarbons in the main condensor of the processing
unit.  
}
\format{
  This data frame contains the following columns:
  \describe{
    \item{purity}{oxygen purity (percentage)}
    \item{hydro}{hydrocarbon (percentage)}
  }
}
\source{
 Montgomery, D.C., Peck, E.A., and Vining, C.G. (2001)
Introduction to Linear Regression Analysis. 3rd Edition, John Wiley and Sons.
}

\examples{
data(p2.7)
attach(p2.7)
purity.lm <- lm(purity ~ hydro)
summary(purity.lm)
# confidence interval for mean purity at 1% hydrocarbon:
predict(purity.lm,newdata=data.frame(hydro = 1.00),interval="confidence")
detach(p2.7)
}
\keyword{datasets}
