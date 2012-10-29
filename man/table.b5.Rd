\name{table.b5}
\alias{table.b5}
\title{Data Set for Table B5}
\usage{data(table.b5)}
\description{
The \code{table.b5} data frame has 27 observations on liquefaction.
}
\format{
  This data frame contains the following columns:
  \describe{
    \item{y}{CO2}
    \item{x1}{Space time (in min)}
    \item{x2}{Temperature (in degrees Celsius)}
    \item{x3}{Percent solvation}
    \item{x4}{Oil yield (g/100g MAF)}
    \item{x5}{Coal total}
    \item{x6}{Solvent total}
    \item{x7}{Hydrogen consumption}
  }
}
\source{
 Montgomery, D.C., Peck, E.A., and Vining, C.G. (2001)
Introduction to Linear Regression Analysis. 3rd Edition, John Wiley and Sons.
}
\references{
(1978) Belle Ayr Liquefaction Runs with Solvent. Industrial Chemical 
Process Design Development, 17, 3.
}
\examples{
data(table.b5)
attach(table.b5)
y.lm <- lm(y ~ x6 + x7)
summary(y.lm)
detach(table.b5)
}
\keyword{datasets}
