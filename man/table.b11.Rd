\name{table.b11}
\alias{table.b11}
\non_function{}
\title{Table B11}
\usage{data(table.b11)}
\description{
The \code{table.b11} data frame has 38 observations on the
quality of Pinot Noir wine.
}
\format{
  This data frame contains the following columns:
  \describe{
    \item{Clarity}{a numeric vector}
    \item{Aroma}{a numeric vector}
    \item{Body}{a numeric vector}
    \item{Flavor}{a numeric vector}
    \item{Oakiness}{a numeric vector}
    \item{Quality}{a numeric vector}
    \item{Region}{a numeric vector}
  }
}
\source{
 Montgomery, D.C., Peck, E.A., and Vining, C.G. (2001)
Introduction to Linear Regression Analysis. 3rd Edition, John Wiley and Sons.
}
\examples{
data(table.b11)
attach(table.b11)
Quality.lm <- lm(Quality ~ Clarity + Aroma + Body + Flavor + Oakiness + 
factor(Region))
summary(Quality.lm)
detach(table.b11)
}
\keyword{datasets}
