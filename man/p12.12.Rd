\name{p12.12}
\alias{p12.12}
\non_function{}
\title{Data Set for Problem 12-12}
\usage{data(p12.12)}
\description{
The \code{p12.12} data frame has 18 observations on an
chemical experiment. A nonlinear model relating concentration to 
reaction time and temperature with an additive error is proposed to
fit these data.  
}
\format{
  This data frame contains the following columns:
  \describe{
    \item{x1}{reaction time (in minutes)}
    \item{x2}{temperature (in degrees Celsius)}
    \item{y}{concentration (in grams/liter)}
  }
}
\source{
 Montgomery, D.C., Peck, E.A., and Vining, C.G. (2001)
Introduction to Linear Regression Analysis. 3rd Edition, John Wiley and Sons.
}
\examples{
data(p12.12)
attach(p12.12)
# fitting the linearized model 
logy.lm <- lm(I(log(y))~I(log(x1))+I(log(x2)))
summary(logy.lm)
plot(logy.lm, which=1)  # checking the residuals
# fitting the nonlinear model
library(nls) # nonlinear regression library
 y.nls <- nls(y ~ theta1*I(x1^theta2)*I(x2^theta3), start=list(theta1=.95, 
theta2=.76, theta3=.21))
 summary(y.nls)
 plot(resid(y.nls)~fitted(y.nls)) # checking the residuals 
}
\keyword{datasets}
