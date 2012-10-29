\name{p5.5}
\alias{p5.5}
\title{Data Set for Problem 5-5}
\usage{data(p5.5)}
\description{
The \code{p5.5} data frame has 14 observations on the average
number of defects per 10000 bottles due to stones in the bottle
wall and the number of weeks since the last furnace overhaul.
}
\format{
  This data frame contains the following columns:
  \describe{
    \item{defects}{a numeric vector}
    \item{weeks}{a numeric vector}
  }
}
\source{
 Montgomery, D.C., Peck, E.A., and Vining, C.G. (2001)
Introduction to Linear Regression Analysis. 3rd Edition, John Wiley and Sons.
}
\examples{
data(p5.5)
defects.lm <- lm(defects~weeks, data=p5.5)
plot(defects.lm, which=1)
}
\keyword{datasets}
