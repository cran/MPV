\name{table.b4}
\alias{table.b4}
\non_function{}
\title{Table B4}
\usage{data(table.b4)}
\description{
The \code{table.b4} data frame has 24 observations on property
valuation.
}
\format{
  This data frame contains the following columns:
  \describe{
    \item{y}{sale price of the house (in thousands of dollars)}
    \item{x1}{taxes (in thousands of dollars)}
    \item{x2}{number of baths}
    \item{x3}{lot size (in thousands of square feet)}
    \item{x4}{living space (in thousands of square feet)}
    \item{x5}{number of garage stalls}
    \item{x6}{number of rooms}
    \item{x7}{number of bedrooms}
    \item{x8}{age of the home (in years)}
    \item{x9}{number of fireplaces}
  }
}
\source{
 Montgomery, D.C., Peck, E.A., and Vining, C.G. (2001)
Introduction to Linear Regression Analysis. 3rd Edition, John Wiley and Sons.
}
\references{
Narula, S.C. and Wellington (1980) Prediction, Linear Regression
and Minimum Sum of Relative Errors. Technometrics, 19, 1977. 
}
\examples{
data(table.b4)
attach(table.b4)
y.lm <- lm(y ~ x1 + x2 + x3 + x4 + x5 + x6 + x7 + x8 + x9)
summary(y.lm)
detach(table.b4)
}
\keyword{datasets}
