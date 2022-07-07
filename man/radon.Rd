\name{radon}
\alias{radon}
\docType{data}
\title{Radon Release}
\description{
Percentage of radon from water released in showers with orifices
of various diameters.  Four replicates were obtained, but it should
be noted that the temperatures for the replicates (in degrees Celsius) are
21, 30, 38, and 46, respectively.  This information should really be
accounted for in any serious analysis of the data.  
}

\usage{data("radon")}
\format{
  A data frame with 15 observations on the following 2 variables.
  \describe{
    \item{\code{diameter}}{shower orifice diameter in mm}
    \item{\code{rep 1}}{percentage radon released in first run}
    \item{\code{rep 2}}{percentage radon released in second run}
    \item{\code{rep 3}}{percentage radon released in third run}
    \item{\code{rep 4}}{percentage radon released in fourth run}
  }
}
\source{Hazin, C.A. and Eichholz, G.G. (1992) Influence of Water
Temperature and Shower Head Orifice Size on the Release of Radon
During Showering, Environment International, 18, 363-369. 
}
\keyword{datasets}
