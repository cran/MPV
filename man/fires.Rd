\name{fires}
\alias{fires}
\title{Micro-fires recorded in a lab setting}
\description{
Rate of spread measurements (inches/s) in each direction: East, West, North 
and South for each of 31 experimental runs at given slopes, measured over 
the given time period of each (measured in seconds).
}
\usage{fires}
\format{
  A data frame with 31 observations on the following 7 variables.
  \describe{
\item{\code{Run}}{numeric}
\item{\code{Slope}}{numeric: vertical rise divided by horizontal run, inclined
from East to West}
\item{\code{ROS_E}}{numeric: rate of spread measured in easterly direction}
\item{\code{ROS_W}}{numeric: rate of spread measured in westerly direction}
\item{\code{ROS_S}}{numeric: rate of spread measured in southerly direction}
\item{\code{ROS_N}}{numeric: rate of spread measured in northerly direction}
\item{\code{Time}}{numeric}
  }
}
\source{Braun, W.J. and Woolford, D.G. (2013) Assessing a stochastic fire spread
simulator.  Journal of Environmental Informatics. 22:1-12.  
}
\keyword{datasets}
