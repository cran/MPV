\name{p8.11}
\alias{p8.11}
\non_function{}
\title{Data Set for Problem 8-11}
\usage{data(p8.11)}
\description{
The \code{p8.11} data frame has 25 observations on the tensile
strength of synthetic fibre used for men's shirts. 
}
\format{
  This data frame contains the following columns:
  \describe{
    \item{y}{tensile strength}
    \item{percent}{percentage of cotton}
  }
}
\source{
 Montgomery, D.C., Peck, E.A., and Vining, C.G. (2001)
Introduction to Linear Regression Analysis. 3rd Edition, John Wiley and Sons.
}
\references{
Montgomery (2001) 
}
\examples{
data(p8.11)
y.lm <- lm(y ~ percent, data=p8.11)
model.matrix(y.lm)
}
\keyword{datasets}
