\name{rectangles}
\alias{rectangles}
\title{Length Measurements on Rectangular Objects}
\description{
Observations of heights, widths and diagonal lengths of several
rectangular objects, such as books, photographs, and so on were measured. 
Only the data in MPV versions 1.62 and later can be trusted; there were
errors in the third column in previous versions.
}
\usage{rectangles}
\format{
  A data frame with 51 observations on the following 4 variables.
  \describe{
\item{\code{h}}{numeric, heights in centimeters}
\item{\code{w}}{numeric, widths in centimeters}
\item{\code{d}}{numeric, diagonal lengths in centimeters}
\item{\code{index}}{numeric, sum of squares of heights and widths}
  }
}
\examples{
x <- sqrt(rectangles$index)
y <- rectangles$d
y.lp <- locpoly(x, y, bandwidth=dpill(x,y), degree=1)
plot(y ~ x)  
lines(y.lp, col=2, lty=2)
abline(0,1) # y = x + measurement error
plot(y.lp$y - y.lp$x, type="l", col=2)
}
\keyword{datasets}
