\name{table.b8}
\alias{table.b8}
\non_function{}
\title{Table B8}
\usage{data(table.b8)}
\description{
The \code{table.b8} data frame has 36 observations on Clathrate
formation.
}
\format{
  This data frame contains the following columns:
  \describe{
    \item{x1}{Amount of surfactant (mass percentage)}
    \item{x2}{Time (min)}
    \item{y}{Clathrate formation (mass percentage)}
  }
}
\source{
 Montgomery, D.C., Peck, E.A., and Vining, C.G. (2001)
Introduction to Linear Regression Analysis. 3rd Edition, John Wiley and Sons.
}
\references{
Tanii, T., Minemoto, M., Nakazawa, K., and Ando, Y. Study on
a Cool Storage System Using HCFC-14 lb Clathrate. Canadian
Journal of Chemical Engineering, 75, 353-360.
}
\examples{
data(table.b8)
attach(table.b8)
clathrate.lm <- lm(y ~ x1 + x2)
summary(clathrate.lm)
detach(table.b8)
}
\keyword{datasets}
