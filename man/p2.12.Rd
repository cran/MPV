\name{p2.12}
\alias{p2.12}
\title{Data Set for Problem 2-12}
\usage{data(p2.12)}
\description{
The \code{p2.12} data frame has 12 observations on 
the number of pounds of steam used per month at a plant and
the average monthly ambient temperature.
}
\format{
  This data frame contains the following columns:
  \describe{
    \item{temp}{ambient temperature (in degrees F)}
    \item{usage}{usage (in thousands of pounds)}
  }
}
\source{
 Montgomery, D.C., Peck, E.A., and Vining, C.G. (2001)
Introduction to Linear Regression Analysis. 3rd Edition, John Wiley and Sons.
}
\examples{
data(p2.12)
attach(p2.12)
usage.lm <- lm(usage ~ temp)
summary(usage.lm)
predict(usage.lm, newdata=data.frame(temp=58), interval="prediction")
detach(p2.12)
}
\keyword{datasets}
