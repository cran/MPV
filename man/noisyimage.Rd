\name{noisyimage}
\alias{noisyimage}
\title{ noisy image }
\usage{data(noisyimage)}
\description{
The \code{noisyimage} is a list.  The third component is 
noisy version of the third component of \link{tarimage}.  
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
with(noisyimage, image(x, y, xy))
}
\keyword{datasets}
