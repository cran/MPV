\name{rftest}
\alias{rftest}
\title{
Pseudorandom Number Testing via Random Forest
}
\description{
Given a sequence of pseudorandom numbers, this function constructs
a random forest prediction model for successive values, based on
previous values up to a given lag.  The ability of the random forest
model to predict future values is inversely related to the quality
of the sequence as an approximation to locally random numbers.
}
\usage{
rftest(u, m=5) 
}
\arguments{
  \item{u}{numeric, a vector of pseudorandom numbers to test}
  \item{m}{numeric, number of lags to test}
}
\value{
  Side effect is a two way layout of graphs showing effectiveness
of prediction on a training and a testing subset of data.  Good 
predictions indicate a poor quality sequence. 
}
\author{
W. John Braun
}
\examples{
    x <- runif(200)
    rftest(x, m = 4)
}
\keyword{graphics}
