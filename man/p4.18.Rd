\name{p4.18}
\alias{p4.18}
\non_function{}
\title{Data Set for Problem 4-18}
\usage{data(p4.18)}
\description{
The \code{p4.18} data frame has 13 observations on an
experiment to produce a synthetic analogue to jojoba oil.
}
\format{
  This data frame contains the following columns:
  \describe{
    \item{x1}{reaction temperature}
    \item{x2}{initial amount of catalyst}
    \item{x3}{pressure}
    \item{y}{yield}
  }
}
\source{
 Montgomery, D.C., Peck, E.A., and Vining, C.G. (2001)
Introduction to Linear Regression Analysis. 3rd Edition, John Wiley and Sons.
}
\reference{
Coteron, Sanchez, Matinez, and Aracil (1993) 
Optimization of the Synthesis of an Analogue of Jojoba Oil Using
a Fully Central Composite Design. Canadian Journal of Chemical
Engineering.
}

\examples{
data(p4.18)
y.lm <- lm(y ~ x1 + x2 + x3, data=p4.18)
summary(y.lm)
y.lm <- lm(y ~ x1, data=p4.18)
}
\keyword{datasets}
