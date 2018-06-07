\name{windWin80}
\alias{windWin80}
\title{Winnipeg Wind Speed}
\usage{data(windWin80)}
\description{
The \code{windWin80} data frame has 366 observations on midnight and noon windspeed
at the Winnipeg International Airport for the year 1980.  
}
\format{
  This data frame contains the following columns:
  \describe{
    \item{h0}{a numeric vector containing the wind speeds at midnight.}
    \item{h12}{a numeric vector containing the wind spees at the following noon.}
  }
}
\examples{
data(windWin80)
ts.plot(windWin80$h12^2)
}
\keyword{datasets}
