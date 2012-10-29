\name{table.b1}
\alias{table.b1}
\title{Table B1}
\usage{data(table.b1)}
\description{
The \code{table.b1} data frame has 28 observations on National
Football League 1976 Team Performance.
}
\format{
  This data frame contains the following columns:
  \describe{
    \item{y}{Games won in a 14 game season}
    \item{x1}{Rushing yards}
    \item{x2}{Passing yards}
    \item{x3}{Punting average (yards/punt)}
    \item{x4}{Field Goal Percentage (FGs made/FGs attempted)}
    \item{x5}{Turnover differential (turnovers acquired - turnovers lost)}
    \item{x6}{Penalty yards}
    \item{x7}{Percent rushing (rushing plays/total plays)}
    \item{x8}{Opponents' rushing yards}
    \item{x9}{Opponents' passing yards}
  }
}
\source{
 Montgomery, D.C., Peck, E.A., and Vining, C.G. (2001)
Introduction to Linear Regression Analysis. 3rd Edition, John Wiley and Sons.
}
\examples{
data(table.b1)
attach(table.b1)
y.lm <- lm(y ~ x2 + x7 + x8)
summary(y.lm)
# over-all F-test:
y.null <- lm(y ~ 1)
anova(y.null, y.lm)
# partial F-test for x7:
y7.lm <- lm(y ~ x2 + x8)
anova(y7.lm, y.lm)
detach(table.b1)
}
\keyword{datasets}
