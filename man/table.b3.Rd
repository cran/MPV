\name{table.b3}
\alias{table.b3}
\non_function{}
\title{Table B3}
\usage{data(table.b3)}
\description{
The \code{table.b3} data frame has observations on gasoline
mileage performance for 32 different automobiles.
}
\format{
  This data frame contains the following columns:
  \describe{
    \item{y}{Miles/gallon}
    \item{x1}{Displacement (cubic in)}
    \item{x2}{Horsepower (ft-lb)}
    \item{x3}{Torque (ft-lb)}
    \item{x4}{Compression ratio}
    \item{x5}{Rear axle ratio}
    \item{x6}{Carburetor (barrels)}
    \item{x7}{No. of transmission speeds}
    \item{x8}{Overall length (in)}
    \item{x9}{Width (in)}
    \item{x10}{Weight (lb)}
    \item{x11}{Type of transmission (1=automatic, 0=manual)}
  }
}
\source{
 Montgomery, D.C., Peck, E.A., and Vining, C.G. (2001)
Introduction to Linear Regression Analysis. 3rd Edition, John Wiley and Sons.
}
\references{
Motor Trend, 1975
}
\examples{
data(table.b3)
attach(table.b3)
y.lm <- lm(y ~ x1 + x6)
summary(y.lm)
# testing for the significance of the regression:
y.null <- lm(y ~ 1)
anova(y.null, y.lm)
# 95% CI for mean gas mileage:
predict(y.lm, newdata=data.frame(x1=275, x6=2), interval="confidence")
# 95% PI for gas mileage:
predict(y.lm, newdata=data.frame(x1=275, x6=2), interval="prediction")
detach(table.b3)
}
\keyword{datasets}
