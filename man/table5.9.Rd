\name{table5.9}
\alias{table5.9}
\title{Table 5.9}
\usage{data(table5.9)}
\description{
The \code{table5.9} data frame has 30 observations on 
wind income (dollars) and corresponding advertising expense.
This is the Restaurant Food Sales Data of Example 5.5.  
}
\format{
  This data frame contains the following columns:
  \describe{
    \item{y}{numeric vector of incomes}
    \item{x}{numeric vector of advertising expenses}
  }
}
\source{
 Montgomery, D.C., Peck, E.A., and Vining, C.G. (2001)
Introduction to Linear Regression Analysis. 3rd Edition, John Wiley and Sons.
}
\examples{
plot(y ~ x, xlab = "expense", ylab = "income", data = table5.9)
# carrying out the calculations in the example to obtain the regression
# weights:
indices <- rep(1:10, c(3, 2, 1, 5, 5, 1, 6, 2, 1, 4))
xbar <- sapply(split(table5.9$x, indices), mean)
yvarhat <- sapply(split(table5.9$y, indices), var)
xbar <- xbar[!is.na(yvarhat)]
yvarhat <- yvarhat[!is.na(yvarhat)]
eg55.lm <- lm(yvarhat ~ xbar)
wts <- 1/predict(eg55.lm, newdata = data.frame(xbar = table5.9$x))
# the values are different from those of the textbook; there seems
# to be some problem with either the calculations or the recorded values
}
\keyword{datasets}
