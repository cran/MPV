\name{p2.15}
\alias{p2.15}
\non_function{}
\title{Data Set for Problem 2-15}
\usage{data(p2.15)}
\description{
The \code{p2.15} data frame has 8 observations on the impact
of temperature on the viscosity of toluene-tetralin blends.
This particular data set deals with blends with a 0.4 molar
fraction of toluene.
}
\format{
  This data frame contains the following columns:
  \describe{
    \item{temp}{temperature (in degrees Celsius)}
    \item{visc}{viscosity (mPa s)}
  }
}
\source{
 Montgomery, D.C., Peck, E.A., and Vining, C.G. (2001)
Introduction to Linear Regression Analysis. 3rd Edition, John Wiley and Sons.
}
\reference{
Byers and Williams (1987) Viscosities of Binary and Ternary Mixtures
of Polynomatic Hydrocarbons. Journal of Chemical and Engineering
Data, 32, 349-354.
}
\examples{
data(p2.15)
attach(p2.15)
plot(visc ~ temp, pch=16)
visc.lm <- lm(visc ~ temp)
plot(visc.lm, which=1)
detach(p2.15)
}
\keyword{datasets}
