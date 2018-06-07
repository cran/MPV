\name{tarimage}
\alias{tarimage}
\title{ target image }
\usage{data(tarimage)}
\description{
The \code{tarimage} is a list.
Most of the values are 0, but there are small regions of 1's.  
}
\format{
  This list contains the following elements:
  \describe{
    \item{x}{a numeric vector having 101 elements.}
    \item{y}{a numeric vector having 101 elements.}
    \item{xy}{a numeric matrix having 101 rows and columns}
  }
}
\examples{
with(tarimage, image(x, y, xy))
}
\keyword{datasets}
