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
GANOVA(dataset, var.equal=TRUE, type="QQ")
}
\arguments{
  \item{dataset}{A data frame, whose first column must be the factor variable
and whose second column must be the response variable.
}
  \item{var.equal}{Logical: if TRUE, within-sample variances are assumed to be equal}
  \item{type}{"QQ" or "hist"}
}
\value{A QQ-plot or a histogram and rugplot}
\author{
W. John Braun
}
\source{Braun, W.J. 2013. Naive Analysis of Variance. Journal of Statistics Education.}

\keyword{graphics}
