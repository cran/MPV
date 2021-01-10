\name{Wpgtemp}
\alias{Wpgtemp}
\title{Winnipeg Maximum Temperatures}
\usage{data(Wpgtemp)}
\description{
The \code{Wpgtemp} data frame has 7671 observations on 
daily maximum temperatures at the Winnipeg International Airport for the years 1960
through 1980.  
}
\format{
  This data frame contains the following columns:
  \describe{
    \item{temperature}{
A numeric vector containing the temperatures in degrees Celsius}
    \item{day}{A numeric vector denoting the observation date in numbers
of days after December 31, 1959}
  }
}
\source{Environment Canada}
\examples{
summary(Wpgtemp)
}
\keyword{datasets}
