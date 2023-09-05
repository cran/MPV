\name{p5.13}
\alias{p5.13}
\title{Data Set for Problem 5-13}
\usage{data(p5.13)}
\description{
The \code{p5.13} data frame has 8 observations on 4 variables, with
responses replicated 3 times.  
}
\format{
  This data frame contains the following columns:
  \describe{
    \item{x1}{numeric}
    \item{x2}{numeric}
    \item{x3}{numeric}
    \item{x4}{numeric}
    \item{y.1}{response 1}
    \item{y.2}{response 2}
    \item{y.3}{response 3}
  }
}
\source{
 Montgomery, D.C., Peck, E.A., and Vining, C.G. (2021)
Introduction to Linear Regression Analysis. 6th Edition, John Wiley and Sons.
}
\references{
Schubert, K., M. W., Kerber, S. R., Schmidt, and Jones, S.E. (1992) 
"The catapult problem; enhanced engineering modeling using experimental 
design," Quality Engineering, 4, 463-473.
}
\examples{
y.lm <- lm(I((y.1+y.2+y.3)/3) ~ x1 + x2 + x3 + x4, data = p5.13)
plot(y.lm, which=1)
}
\keyword{datasets}
