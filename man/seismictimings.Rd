\name{seismictimings}
\alias{seismictimings}
\title{Seismic Timing Data}
\usage{seismictimings}
\description{
The \code{seismictimings} data frame has 504 rows and 3 columns.
Thickness of a layer of Alberta substratum as measured by
several transects of geophones.
}
\format{
  This data frame contains the following columns:
  \describe{
    \item{x}{longitudinal coordinate of geophone.}
    \item{y}{latitudinal coordinate of geophone.}
    \item{z}{time for signal to pass through substratum.}
  }
}
\examples{
plot(y ~ x, data = seismictimings)
}
\keyword{datasets}
