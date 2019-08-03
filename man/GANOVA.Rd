\name{GANOVA}
\alias{GANOVA}
\title{
Graphical ANOVA Plot
}
\description{
Graphical analysis of one-way ANOVA data.  It 
allows visualization of the usual F-test.
}
\usage{
GANOVA(dataset, var.equal=TRUE, type="QQ", center=TRUE, shift=0)
}
\arguments{
  \item{dataset}{A data frame, whose first column must be the factor variable
and whose second column must be the response variable.
}
  \item{var.equal}{Logical: if TRUE, within-sample variances are assumed to be equal}
  \item{type}{"QQ" or "hist"}
  \item{center}{if TRUE, center and scale the means to match the scale of the errors}
  \item{shift}{on the histogram, lift the points representing the means above the horizontal axis by this amount.}
}
\value{A QQ-plot or a histogram and rugplot}
\author{
W. John Braun and Sarah MacQueen
}
\source{Braun, W.J. 2013. Naive Analysis of Variance. Journal of Statistics Education.}

\keyword{graphics}
