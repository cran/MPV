\name{Uplot}
\alias{Uplot}
\title{
Plot of Multipliers in Regression ANOVA Plot
}
\description{
This function graphically displays the coefficient multipliers
used in the Regression Plot for the given predictor.
}
\usage{
Uplot(X.qr, Xcolumn = 1)
}
\arguments{
  \item{X.qr}{The design matrix or the QR decomposition of the design matrix.
}
  \item{Xcolumn}{The column(s) of the design matrix under study; this
can be either integer valued or a character string.}
}
\value{A bar plot is displayed.}
\author{
W. John Braun
}
\examples{
# Jojoba oil data set
X <- p4.18[,-4]
Uplot(X, 1:4)
# NFL data set; see Qyplot result first
X <- table.b1[,-1]
Uplot(X, c(2,3,9))
# In this example, x8 is the only predictor in
# the true model:
X <- pathoeg[,-10]
y <- pathoeg[,10]
pathoeg.F <- Qyplot(X, y, plotIt=FALSE)
Uplot(X, "x8")
Uplot(X, 9) # same as above
Uplot(pathoeg.F$QR, 9) # same as above
X <- table.b1[,-1]
Uplot(X, c("x2", "x3", "x9"))
}
\keyword{graphics}
