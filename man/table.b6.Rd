\name{table.b6}
\alias{table.b6}
\non_function{}
\title{Data Set for Table B6}
\usage{data(table.b6)}
\description{
The \code{table.b6} data frame has 28 observations on 
a tube-flow reactor.
}
\format{
  This data frame contains the following columns:
  \describe{
    \item{y}{Nb0Cl3 concentration (g-mol/l)}
    \item{x1}{COCl2 concentration (g-mol/l)}
    \item{x2}{Space time (s)}
    \item{x3}{Molar density (g-mol/l)}
    \item{x4}{Mole fraction CO2}
  }
}
\source{
 Montgomery, D.C., Peck, E.A., and Vining, C.G. (2001)
Introduction to Linear Regression Analysis. 3rd Edition, John Wiley and Sons.
}
\references{
(1972) Kinetics of Chlorination of Niobium oxychloride by
Phosgene in a Tube-Flow Reactor.  Industrial and Engineering
Chemistry, Process Design Development, 11(2).
}
\examples{
data(table.b6)
# Partial Solution to Problem 3.9
attach(table.b6)
y.lm <- lm(y ~ x1 + x4)
summary(y.lm)
detach(table.b6)
}
\keyword{datasets}
