\name{p2.10}
\alias{p2.10}
\non_function{}
\title{Data Set for Problem 2-10}
\usage{data(p2.10)}
\description{
The \code{p2.10} data frame has 26 observations on weight and
systolic blood pressure for randomly selected males in the 25-30
age group.  
}
\format{
  This data frame contains the following columns:
  \describe{
    \item{weight}{in pounds}
    \item{sysbp}{systolic blood pressure}
  }
}
\source{
 Montgomery, D.C., Peck, E.A., and Vining, C.G. (2001)
Introduction to Linear Regression Analysis. 3rd Edition, John Wiley and Sons.
}
\examples{
data(p2.10)
attach(p2.10)
cor.test(weight, sysbp, method="pearson")  # tests rho=0
                                           # and computes 95% CI for rho
                                           # using Fisher's Z-transform
}
\keyword{datasets}
