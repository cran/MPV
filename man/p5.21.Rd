\name{p5.21}
\alias{p5.21}
\title{Data Set for Problem 5-21}
\usage{data(p5.21)}
\description{
The \code{p5.21} data frame has 4 observations on
2 variables (replicated 4 times).
}
\format{
  This data frame contains the following columns:
  \describe{
    \item{Mix.Rate}{a numeric vector}
    \item{y1}{a numeric vector}
    \item{y2}{a numeric vector}
    \item{y3}{a numeric vector}
    \item{y4}{a numeric vector}
  }
}
\source{
 Montgomery, D.C., Peck, E.A., and Vining, C.G. (2021)
Introduction to Linear Regression Analysis. 6th Edition, John Wiley and Sons.
}
\examples{
 cementStrength <- reshape(p5.21, idvar = "Mix.Rate", varying=list(2:5), 
    direction="long", v.names=c("TensileStrength"))
 rownames(cementStrength) <- NULL
 anova(lm(TensileStrength ~ Mix.Rate*time, data = cementStrength))
}
\keyword{datasets}
