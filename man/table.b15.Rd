\name{table.b15}
\alias{table.b15}
\title{Table B15 - Air Pollution and Mortality Data}
\usage{data(table.b15)}
\description{
The \code{table.b15} data frame has 60 observations on the mortality, environment, and demographic variables for a sample of American cities.  
}
\format{
  This data frame contains the following columns:
  \describe{
    \item{City}{character vector}
    \item{Mort}{numeric vector, age-adjusted mortality from all causes per 100000}
    \item{Precip}{numeric vector, precipitation in inches}
    \item{Educ}{numeric vector, median number of school years completed}
    \item{Nonwhite}{numeric vector, percentage of 1960 population that is nonwhite}
    \item{Nox}{numeric vector, relative pollution potential of nitrous oxides}
    \item{SO2}{numeric vector, relative pollution potential of sulfur dioxide}
  }
}
\source{
 Montgomery, D.C., Peck, E.A., and Vining, C.G. (2021)
Introduction to Linear Regression Analysis. 6th Edition, John Wiley and Sons.
}
\references{
McDonald, G. C. and Ayers, J.A. [1978], "Some applications of Chernuff faces: 
A technique for graphically representing multivariate data", in 
Graphical Representation of Multivariate Data, Academic Press, New York.
}
\examples{
data(table.b15)
pairs(table.b15[,-1])
}
\keyword{datasets}
