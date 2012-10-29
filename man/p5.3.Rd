\name{p5.3}
\alias{p5.3}
\title{Data Set for Problem 5-3}
\usage{data(p5.3)}
\description{
The \code{p5.3} data frame has 12 observations on the
number of bacteria surviving in a canned food product and the
number of minutes of exposure to 300 degree Fahrenheit heat.
}
\format{
  This data frame contains the following columns:
  \describe{
    \item{bact}{number of surviving bacteria}
    \item{min}{number of minutes of exposure}
  }
}
\source{
 Montgomery, D.C., Peck, E.A., and Vining, C.G. (2001)
Introduction to Linear Regression Analysis. 3rd Edition, John Wiley and Sons.
}
\examples{
data(p5.3)
plot(bact~min, data=p5.3)
}
\keyword{datasets}
