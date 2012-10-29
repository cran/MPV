\name{p2.13}
\alias{p2.13}
\title{Data Set for Problem 2-13}
\usage{data(p2.13)}
\description{
The \code{p2.13} data frame has 16 observations on the number
of days the ozone levels exceeded 0.2 ppm in the
South Coast Air Basin of California for the years 1976 through
1991.  It is believed that these levels are related to temperature.
}
\format{
  This data frame contains the following columns:
  \describe{
    \item{days}{number of days ozone levels exceeded 0.2 ppm}
    \item{index}{a seasonal
meteorological index giving the seasonal average 850 millibar
temperature.}
  }
}
\source{
 Montgomery, D.C., Peck, E.A., and Vining, C.G. (2001)
Introduction to Linear Regression Analysis. 3rd Edition, John Wiley and Sons.
}
\references{
Davidson, A. (1993) Update on Ozone Trends in California's South Coast Air 
Basin. \emph{Air Waste}, 43, 226-227.
}
\examples{
data(p2.13)
attach(p2.13)
plot(days~index, ylim=c(-20,130))
ozone.lm <- lm(days ~ index)
summary(ozone.lm)
# plots of confidence and prediction intervals:
ozone.conf <- predict(ozone.lm, interval="confidence")
lines(sort(index), ozone.conf[order(index),2], col="red")
lines(sort(index), ozone.conf[order(index),3], col="red")
ozone.pred <- predict(ozone.lm, interval="prediction")
lines(sort(index), ozone.pred[order(index),2], col="blue")
lines(sort(index), ozone.pred[order(index),3], col="blue")
detach(p2.13)
}
\keyword{datasets}
