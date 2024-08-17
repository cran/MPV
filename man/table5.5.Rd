\name{table5.5}
\alias{table5.5}
\title{Table 5.5}
\usage{data(table5.5)}
\description{
The \code{table5.5} data frame has 25 observations on 
wind velocity (mph) and corresponding DC output from a windmill turbine.
This is the Windmill Data of Example 5.2.  
}
\format{
  This data frame contains the following columns:
  \describe{
    \item{v}{numeric vector of velocities}
    \item{DC}{numeric vector of DC output values}
  }
}
\source{
 Montgomery, D.C., Peck, E.A., and Vining, C.G. (2001)
Introduction to Linear Regression Analysis. 3rd Edition, John Wiley and Sons.
}
\examples{
plot(DC ~ v, data = table5.5)
}
\keyword{datasets}
