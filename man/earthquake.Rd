\name{earthquake}
\alias{earthquake}
\title{Earthquakes Data}
\usage{earthquake}
\description{
The earthquake data frame contains measurements of latitude,
longitude, focal depth and magnitude for all earthquakes
having magnitude greater than 5.8 between 1964 and 1985.
}
\format{
  This data frame contains 2178 observations on the following columns:
  \describe{
    \item{depth}{numeric vector of focal depths.}
    \item{latitude}{latitudinal coordinate.}
    \item{longitude}{longitudinal coordinate.}
    \item{magnitude}{numeric vector of magnitudes.}
  }
}

\source{
Jeffrey S. Simonoff (1996), Smoothing Methods in Statistics,  
Springer-Verlag, New York.  
}
\examples{
summary(earthquake)
}
\keyword{datasets}

