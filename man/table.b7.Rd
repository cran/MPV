\name{table.b7}
\alias{table.b7}
\non_function{}
\title{Data Set for Table B7}
\usage{data(table.b7)}
\description{
The \code{table.b7} data frame has 16 observations on 
oil extraction from peanuts.
}
\format{
  This data frame contains the following columns:
  \describe{
    \item{x1}{CO2 pressure (bar)}
    \item{x2}{CO2 temperature (in degrees Celsius)}
    \item{x3}{peanut moisture (percent by weight)}
    \item{x4}{CO2 flow rate (L/min)}
    \item{x5}{peanut particle size (mm)}
    \item{y}{total oil yield}
  }
}
\source{
 Montgomery, D.C., Peck, E.A., and Vining, C.G. (2001)
Introduction to Linear Regression Analysis. 3rd Edition, John Wiley and Sons.
}
\references{
Kilgo, M.B.  An Application of Fractional Experimental Designs. 
Quality Engineering, 1, 19-23.
}
\examples{
data(table.b7)
attach(table.b7)
# partial solution to Problem 3.11:
peanuts.lm <- lm(y ~ x1 + x2 + x3 + x4 + x5)
summary(peanuts.lm)
detach(table.b7)
}
\keyword{datasets}
