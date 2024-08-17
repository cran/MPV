\name{table5.2}
\alias{table5.2}
\title{Table 5.2}
\usage{data(table5.2)}
\description{
The \code{table5.2} data frame has 53 observations on 
energy usage (KWH) and corresponding demand (KW) at a sample of 
residences.
This is the Electric Utility Data of Example 5.1.  
}
\format{
  This data frame contains the following columns:
  \describe{
    \item{Customer}{a numeric vector of customer IDs}
    \item{x}{a numeric vector of energy usage values}
    \item{y}{a numeric vector of demand values}
  }
}
\source{
 Montgomery, D.C., Peck, E.A., and Vining, C.G. (2001)
Introduction to Linear Regression Analysis. 3rd Edition, John Wiley and Sons.
}
\examples{
plot(y ~ x, xlab = "Usage", ylab = "Demand", data = table5.2)
anova(lm(y ~ x, data = table5.2)) # Note the typo in Table 5.3 for SS Regression
}
\keyword{datasets}
